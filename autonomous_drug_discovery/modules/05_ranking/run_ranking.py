"""
Module 05: Ranking — Multi-criteria final ranker.

Wired into the orchestrator (`orchestrator.py run` + standalone `rank` subcommand).
Combines docking score, ADMET flags, and (optional) AiZynthFinder retrosynthetic
feasibility into a single composite score and a ranked scorecard for
medicinal-chemistry expert review.

Composite weights: 0.5 docking + 0.3 ADMET + 0.2 synthesis (see DEFAULT_WEIGHTS).

Synthesis scoring: when `--aizynth_config` is supplied (or the AIZYNTH_CONFIG
environment variable is set), the top-N provisional candidates — ranked on
docking + ADMET alone — are scored with AiZynthFinder, a retrosynthetic
tree-search planner that runs in the separate `aizynth_env`. A molecule with a
route to purchasable building blocks scores high (fewer steps = better); one
with no route scores low. AiZynth is slow (~30-60 s/molecule), so it is only
applied to the top-N contenders; the rest keep a neutral 0.5. Without a config,
every candidate keeps 0.5 and synthesis has no effect on the ranking.
See scripts/aizynth_score.py (the in-env worker) and envs/env_aizynth.yml.

Input contract:
  - docking_results.csv  (from Stage 4) — required
  - screening_report.json (from Stage 3) — optional, for ADMET enrichment

Output contract:
  - ranked_candidates.json + run_metadata.json
"""

import os
import sys
import csv
import json
import argparse
import subprocess
import tempfile
import traceback
from datetime import datetime, timezone
from pathlib import Path

MODULE_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = MODULE_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from telemetry import TelemetryDB

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

# Lead-like heavy-atom window for the size-desirability factor. Vina affinity
# scales with molecule size, so a raw docking score rewards fragments. The
# docking sub-score is multiplied by a trapezoidal factor that is 1.0 inside
# [SIZE_PLATEAU_LO, SIZE_PLATEAU_HI] and ramps to 0 outside [SIZE_MIN, SIZE_MAX]
# — a 12-heavy-atom biphenyl that slipped past screening earns almost no
# docking credit even if its raw affinity looks strong.
SIZE_MIN, SIZE_PLATEAU_LO, SIZE_PLATEAU_HI, SIZE_MAX = 12.0, 20.0, 35.0, 50.0


DEFAULT_WEIGHTS = {
    "docking": 0.5,      # binding affinity (kcal/mol, more negative = better)
    "admet": 0.3,        # composite ADMET safety score
    "synthesis": 0.2,    # retrosynthetic feasibility (AiZynthFinder, top-N only)
}

# Synthesis scoring is expensive (~30-60 s/molecule), so AiZynthFinder is run
# only on the strongest candidates by docking + ADMET; the rest keep PLACEHOLDER.
DEFAULT_SYNTHESIS_TOP_N = 15
PLACEHOLDER_SYNTH = 0.5
# scripts/ sits at the repo root, one level above the application package.
SYNTH_WORKER = PROJECT_ROOT.parent / "scripts" / "aizynth_score.py"


def _load_docking_results(csv_path: Path) -> list[dict]:
    with open(csv_path) as f:
        return list(csv.DictReader(f))


def _load_screening_report(json_path: Path) -> dict[str, dict]:
    """Return SMILES → screening record (incl. ADMET annotations).

    Stage 3 writes survivors under the `survivors` key (each a flat dict of
    `smiles` + `desc_*` / `admet_*` properties). Older drafts used `molecules`;
    both are accepted so a stale report doesn't silently drop ADMET data.
    """
    if not json_path.exists():
        return {}
    with open(json_path) as f:
        report = json.load(f)
    entries = report.get("survivors") or report.get("molecules") or []
    by_smiles = {}
    for entry in entries:
        smi = entry.get("smiles")
        if smi:
            by_smiles[smi] = entry
    return by_smiles


def _admet_score(screening_record: dict) -> float:
    """Composite ADMET score in [0, 1]; higher = safer.

    Reads the flat `admet_hERG` / `admet_AMES` / `admet_DILI` toxicity
    probabilities that Stage 3 stamps onto each survivor (higher prob = more
    toxic, so safety = 1 - prob). Returns 0.5 (neutral) when ADMET data is
    missing — unscreened molecules aren't penalised, but they don't earn
    safety credit they haven't proven either.
    """
    if not screening_record:
        return 0.5
    risks = []
    for key in ("admet_hERG", "admet_AMES", "admet_DILI"):
        v = screening_record.get(key)
        if isinstance(v, (int, float)):
            risks.append(1.0 - float(v))
    return sum(risks) / len(risks) if risks else 0.5


def _heavy_atom_count(smiles: str) -> int | None:
    """Heavy-atom count from SMILES, or None if RDKit is unavailable / parse fails."""
    if not HAS_RDKIT or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol is not None else None


def _size_desirability(heavy_atoms: int | None) -> float:
    """Trapezoidal lead-like size factor in [0, 1].

    1.0 inside the plateau, linearly ramping to 0 at SIZE_MIN / SIZE_MAX.
    Unknown size (RDKit unavailable) returns 1.0 — no penalty without evidence.
    """
    if heavy_atoms is None:
        return 1.0
    ha = float(heavy_atoms)
    if ha <= SIZE_MIN or ha >= SIZE_MAX:
        return 0.0
    if ha < SIZE_PLATEAU_LO:
        return (ha - SIZE_MIN) / (SIZE_PLATEAU_LO - SIZE_MIN)
    if ha > SIZE_PLATEAU_HI:
        return (SIZE_MAX - ha) / (SIZE_MAX - SIZE_PLATEAU_HI)
    return 1.0


def _run_aizynth(candidates: list[dict], aizynth_config: str,
                 conda_bin: str = "conda", env_name: str = "aizynth_env") -> dict[str, dict]:
    """Score `candidates` with AiZynthFinder via the `aizynth_env` worker.

    `candidates` is a list of {ligand_id, smiles}. Returns {smiles: synth_dict}
    where synth_dict carries `solved`, `n_steps`, `top_score`, `synthesis_score`.

    Returns {} (with a logged reason) if the worker, config, or env is missing
    or the run fails — the caller then falls back to the neutral placeholder, so
    a missing AiZynthFinder setup degrades ranking gracefully instead of crashing.
    """
    if not SYNTH_WORKER.exists():
        print(f"[Ranking] AiZynth worker not found at {SYNTH_WORKER}; "
              "skipping synthesis scoring.")
        return {}
    cfg = Path(aizynth_config).expanduser()
    if not cfg.exists():
        print(f"[Ranking] AiZynth config not found at {cfg}; skipping synthesis scoring.")
        return {}

    with tempfile.TemporaryDirectory() as tmp:
        in_path = Path(tmp) / "candidates.json"
        out_path = Path(tmp) / "synthesis.json"
        in_path.write_text(json.dumps(candidates))
        cmd = [
            conda_bin, "run", "--no-capture-output", "-n", env_name, "python",
            str(SYNTH_WORKER),
            "--config", str(cfg),
            "--input", str(in_path),
            "--output", str(out_path),
        ]
        print(f"[Ranking] Scoring {len(candidates)} candidates with AiZynthFinder "
              f"(env={env_name}) — this is slow, ~30-60 s/molecule ...")
        try:
            subprocess.check_call(cmd)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"[Ranking] AiZynthFinder scoring failed ({e}); "
                  "using placeholder synthesis scores.")
            return {}
        if not out_path.exists():
            print("[Ranking] AiZynthFinder produced no output; "
                  "using placeholder synthesis scores.")
            return {}
        return json.loads(out_path.read_text())


def _normalise_docking(affinity: float, best: float, worst: float) -> float:
    if worst == best:
        return 0.5
    return (worst - affinity) / (worst - best)


def rank_candidates(
    docking_csv: Path,
    screening_json: Path | None,
    weights: dict[str, float] = DEFAULT_WEIGHTS,
    aizynth_config: str | None = None,
    synthesis_top_n: int = DEFAULT_SYNTHESIS_TOP_N,
    conda_bin: str = "conda",
) -> list[dict]:
    """Rank docked candidates by a docking + ADMET + synthesis composite.

    When `aizynth_config` is given, the top `synthesis_top_n` candidates by a
    provisional docking + ADMET score are scored with AiZynthFinder; the rest
    keep a neutral placeholder. Synthesis therefore only re-orders the leading
    contenders — exactly where it matters and where the cost is affordable.
    """
    rows = _load_docking_results(docking_csv)
    if not rows:
        return []

    screening = _load_screening_report(screening_json) if screening_json else {}

    affinities = [float(r["affinity"]) for r in rows if r.get("affinity") not in (None, "")]
    best, worst = min(affinities), max(affinities)

    # Pass 1 — docking + ADMET for every candidate.
    # The docking sub-score is the raw-affinity normalisation scaled by a
    # lead-like size factor, so a fragment that docks well per-atom cannot
    # out-rank a properly sized lead with comparable affinity.
    partials = []
    for r in rows:
        smiles = r.get("smiles", "")
        affinity = float(r["affinity"])
        heavy_atoms = _heavy_atom_count(smiles)
        size_factor = _size_desirability(heavy_atoms)
        dock_raw = _normalise_docking(affinity, best, worst)
        ligand_eff = (round(-affinity / heavy_atoms, 4)
                      if heavy_atoms else None)
        partials.append({
            "ligand_id": r.get("ligand_id"),
            "smiles": smiles,
            "affinity": affinity,
            "heavy_atoms": heavy_atoms,
            "ligand_efficiency": ligand_eff,
            "size_factor": round(size_factor, 4),
            "dock": dock_raw * size_factor,
            "admet": _admet_score(screening.get(smiles, {})),
        })

    # Pass 2 — AiZynthFinder synthesis scoring on the top-N provisional set.
    synth_map: dict[str, dict] = {}
    synthesis_enabled = bool(aizynth_config)
    if synthesis_enabled:
        prov_denom = weights["docking"] + weights["admet"] or 1.0

        def _provisional(p):
            return (weights["docking"] * p["dock"] + weights["admet"] * p["admet"]) / prov_denom

        top_candidates: list[dict] = []
        seen_smiles: set[str] = set()
        for p in sorted(partials, key=_provisional, reverse=True):
            if p["smiles"] and p["smiles"] not in seen_smiles:
                seen_smiles.add(p["smiles"])
                top_candidates.append({"ligand_id": p["ligand_id"], "smiles": p["smiles"]})
            if len(top_candidates) >= synthesis_top_n:
                break
        if top_candidates:
            synth_map = _run_aizynth(top_candidates, aizynth_config, conda_bin)

    # Pass 3 — combine into the final composite score.
    ranked = []
    for p in partials:
        synth_rec = synth_map.get(p["smiles"])
        if synth_rec:
            synth = synth_rec["synthesis_score"]
            if synth_rec.get("error"):
                status = f"AiZynthFinder: search error ({synth_rec['error']})"
            elif synth_rec.get("solved"):
                status = f"AiZynthFinder: route found ({synth_rec['n_steps']} steps)"
            else:
                status = "AiZynthFinder: no route to purchasable building blocks"
        elif synthesis_enabled:
            synth = PLACEHOLDER_SYNTH
            status = "not evaluated (outside top-N synthesis triage)"
        else:
            synth = PLACEHOLDER_SYNTH
            status = "placeholder (AiZynthFinder not run)"

        combined = (
            weights["docking"] * p["dock"]
            + weights["admet"] * p["admet"]
            + weights["synthesis"] * synth
        )

        ranked.append({
            "ligand_id": p["ligand_id"],
            "smiles": p["smiles"],
            "docking_affinity": p["affinity"],
            "heavy_atoms": p["heavy_atoms"],
            "ligand_efficiency": p["ligand_efficiency"],
            "size_factor": p["size_factor"],
            "scores": {"docking": p["dock"], "admet": p["admet"], "synthesis": synth},
            "combined_score": combined,
            "synthesis_status": status,
        })

    ranked.sort(key=lambda x: x["combined_score"], reverse=True)
    for i, entry in enumerate(ranked, start=1):
        entry["final_rank"] = i
    return ranked


def main():
    parser = argparse.ArgumentParser(description="Module 05: Multi-criteria final ranking")
    parser.add_argument("--docking_csv", required=True, help="Path to docking_results.csv")
    parser.add_argument("--screening_json", default=None, help="Path to screening_report.json (optional)")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--db_path", default=None, help="Path to telemetry database")
    parser.add_argument("--campaign_id", default=None, help="Campaign ID for telemetry")
    parser.add_argument("--aizynth_config", default=os.environ.get("AIZYNTH_CONFIG"),
                        help="AiZynthFinder config.yml — enables synthesis scoring of "
                             "the top-N candidates (default: $AIZYNTH_CONFIG)")
    parser.add_argument("--synthesis_top_n", type=int, default=DEFAULT_SYNTHESIS_TOP_N,
                        help="How many top candidates to score with AiZynthFinder")
    args = parser.parse_args()

    docking_csv = Path(args.docking_csv)
    screening_json = Path(args.screening_json) if args.screening_json else None
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    conda_bin = os.environ.get("CONDA_EXE", "conda")
    synthesis_on = bool(args.aizynth_config)

    timestamp = datetime.now(timezone.utc).isoformat()
    db = None
    run_id = None
    if args.db_path and args.campaign_id:
        db = TelemetryDB(args.db_path)
        run_id = db.start_run(
            campaign_id=args.campaign_id,
            module_name="05_ranking",
            input_path=str(docking_csv),
            parameters={
                "weights": DEFAULT_WEIGHTS,
                "aizynth_config": args.aizynth_config,
                "synthesis_top_n": args.synthesis_top_n if synthesis_on else None,
            },
        )

    try:
        ranked = rank_candidates(
            docking_csv, screening_json,
            aizynth_config=args.aizynth_config,
            synthesis_top_n=args.synthesis_top_n,
            conda_bin=conda_bin,
        )
        output_path = out_dir / "ranked_candidates.json"
        with open(output_path, "w") as f:
            json.dump({
                "timestamp": timestamp,
                "weights": DEFAULT_WEIGHTS,
                "n_candidates": len(ranked),
                "candidates": ranked,
            }, f, indent=2)

        with open(out_dir / "run_metadata.json", "w") as f:
            json.dump({
                "module": "05_ranking",
                "timestamp": timestamp,
                "input_docking_csv": str(docking_csv),
                "input_screening_json": str(screening_json) if screening_json else None,
                "output": str(output_path),
                "n_candidates": len(ranked),
                "status": "success",
                "synthesis_scoring": (
                    f"AiZynthFinder on top {args.synthesis_top_n} candidates"
                    if synthesis_on else
                    "disabled — synthesis weight contributes a constant 0.5 "
                    "(pass --aizynth_config to enable)"
                ),
            }, f, indent=2)

        print(f"[Ranking] Ranked {len(ranked)} candidates → {output_path}")
        if db and run_id:
            db.complete_run(run_id, "success", str(output_path))

    except Exception as e:
        error_msg = traceback.format_exc()
        print(f"[Ranking] FAILED: {e}")
        if db and run_id:
            db.complete_run(run_id, "failed", error_trace=error_msg)
        sys.exit(1)
    finally:
        if db:
            db.close()


if __name__ == "__main__":
    main()
