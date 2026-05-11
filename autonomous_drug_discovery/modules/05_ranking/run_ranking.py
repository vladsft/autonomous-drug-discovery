"""
Module 05: Ranking — Multi-criteria final ranker.

Wired into the orchestrator (`orchestrator.py run` + standalone `rank` subcommand).
Combines docking score, ADMET flags, and (optional) AiZynthFinder retrosynthetic
feasibility into a single composite score and a ranked scorecard for
medicinal-chemistry expert review.

Composite weights: 0.5 docking + 0.3 ADMET + 0.2 synthesis (see DEFAULT_WEIGHTS).

Synthesis scoring: the `_synthesis_score()` slot currently returns a constant
0.5 (neutral) — AiZynthFinder integration via the `aizynth_env` is a deferred
add-on (it requires a per-molecule subprocess + ~30-60 s/molecule even on GPU,
which is too slow to run in the default pipeline). To integrate, replace
`_synthesis_score()` with a subprocess call into the aizynth_env.

Input contract:
  - docking_results.csv  (from Stage 4) — required
  - screening_report.json (from Stage 3) — optional, for ADMET enrichment

Output contract:
  - ranked_candidates.json + run_metadata.json
"""

import sys
import csv
import json
import argparse
import traceback
from datetime import datetime, timezone
from pathlib import Path

MODULE_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = MODULE_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from telemetry import TelemetryDB


DEFAULT_WEIGHTS = {
    "docking": 0.5,      # binding affinity (kcal/mol, more negative = better)
    "admet": 0.3,        # composite ADMET safety score
    "synthesis": 0.2,    # retrosynthetic feasibility (AiZynth, future)
}


def _load_docking_results(csv_path: Path) -> list[dict]:
    with open(csv_path) as f:
        return list(csv.DictReader(f))


def _load_screening_report(json_path: Path) -> dict[str, dict]:
    """Return SMILES → screening record (incl. ADMET annotations)."""
    if not json_path.exists():
        return {}
    with open(json_path) as f:
        report = json.load(f)
    by_smiles = {}
    for entry in report.get("molecules", []):
        smi = entry.get("smiles")
        if smi:
            by_smiles[smi] = entry
    return by_smiles


def _admet_score(screening_record: dict) -> float:
    """Composite ADMET score in [0, 1]; higher = safer.

    Uses ADMET-AI flags when present. Returns 0.5 (neutral) when ADMET data
    is missing — i.e. unscreened molecules don't get penalised, but they
    also don't earn safety credit they haven't proven.
    """
    admet = screening_record.get("admet") if screening_record else None
    if not admet:
        return 0.5
    risks = []
    for key in ("hERG", "AMES", "DILI"):
        v = admet.get(key)
        if isinstance(v, (int, float)):
            risks.append(1.0 - float(v))
    return sum(risks) / len(risks) if risks else 0.5


def _synthesis_score(smiles: str) -> float:
    """Retrosynthetic feasibility — PLACEHOLDER.

    Real implementation will call AiZynthFinder here on the top-N survivors:
        from aizynthfinder.aizynthfinder import AiZynthFinder
        finder = AiZynthFinder(configdict)
        finder.target_smiles = smiles
        finder.tree_search()
        return _normalise(finder.routes)

    Until then, returns a neutral 0.5 so synthesis weight has no effect on
    ranking. Do not interpret current ranks as synthesis-aware.
    """
    return 0.5


def _normalise_docking(affinity: float, best: float, worst: float) -> float:
    if worst == best:
        return 0.5
    return (worst - affinity) / (worst - best)


def rank_candidates(
    docking_csv: Path,
    screening_json: Path | None,
    weights: dict[str, float] = DEFAULT_WEIGHTS,
) -> list[dict]:
    rows = _load_docking_results(docking_csv)
    if not rows:
        return []

    screening = _load_screening_report(screening_json) if screening_json else {}

    affinities = [float(r["affinity"]) for r in rows if r.get("affinity") not in (None, "")]
    best, worst = min(affinities), max(affinities)

    ranked = []
    for r in rows:
        smiles = r.get("smiles", "")
        affinity = float(r["affinity"])
        screening_record = screening.get(smiles, {})

        dock = _normalise_docking(affinity, best, worst)
        admet = _admet_score(screening_record)
        synth = _synthesis_score(smiles)

        combined = (
            weights["docking"] * dock
            + weights["admet"] * admet
            + weights["synthesis"] * synth
        )

        ranked.append({
            "ligand_id": r.get("ligand_id"),
            "smiles": smiles,
            "docking_affinity": affinity,
            "scores": {"docking": dock, "admet": admet, "synthesis": synth},
            "combined_score": combined,
            "synthesis_status": "placeholder (AiZynthFinder not yet integrated)",
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
    args = parser.parse_args()

    docking_csv = Path(args.docking_csv)
    screening_json = Path(args.screening_json) if args.screening_json else None
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now(timezone.utc).isoformat()
    db = None
    run_id = None
    if args.db_path and args.campaign_id:
        db = TelemetryDB(args.db_path)
        run_id = db.start_run(
            campaign_id=args.campaign_id,
            module_name="05_ranking",
            input_path=str(docking_csv),
            parameters={"weights": DEFAULT_WEIGHTS},
        )

    try:
        ranked = rank_candidates(docking_csv, screening_json)
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
                "notes": "AiZynthFinder synthesis scoring not yet integrated; synthesis weight contributes a constant.",
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
