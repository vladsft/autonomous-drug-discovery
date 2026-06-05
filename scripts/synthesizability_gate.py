#!/usr/bin/env python3
"""Synthesizability gate — Stage 2.5, between generation/screening and docking.

The pipeline's blind spot (measured 2026-05-30): 3D generators (TargetDiff/
Pocket2Mol) produce pocket-fitting molecules with no synthesizability prior, so
a large fraction cannot be made from purchasable building blocks — 0/28 top
kinase candidates had a retrosynthetic route. Docking + ranking then faithfully
score molecules that can't exist in a flask.

This gate runs AiZynthFinder retrosynthesis on a molecule set BEFORE the
expensive docking stage and keeps only the makeable ones, so docking and the
chemist's attention are spent on actionable matter. It is the upstream
counterpart to the post-hoc synthesis score already folded into Stage 5 ranking.

Flow:
  SDF in -> extract (molecule_id, smiles) -> aizynth_score.py (in aizynth_env)
         -> keep solved (or synthesis_score >= --min-score) -> SDF out + report.

Usage:
  python scripts/synthesizability_gate.py \
      --input  data/<campaign>/screened/screened_molecules.sdf \
      --config ~/aizynthfinder_data/config.yml \
      --output data/<campaign>/gated/makeable.sdf \
      --report data/<campaign>/gated/gate_report.json \
      [--require-solved | --min-score 0.30]

Runs the orchestrator side in `base`; shells out to `aizynth_env` for the
search exactly as Stage 5 does, so the ML stacks never co-resolve.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
AIZYNTH_SCORER = REPO_ROOT / "scripts" / "aizynth_score.py"


def _read_sdf(path: Path):
    """Yield (molecule_id, smiles, rdkit_mol) for each record."""
    from rdkit import Chem
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    for mol in supplier:
        if mol is None:
            continue
        mid = mol.GetProp("molecule_id") if mol.HasProp("molecule_id") else (
            mol.GetProp("_Name") if mol.HasProp("_Name") else None)
        smi = mol.GetProp("smiles") if mol.HasProp("smiles") else Chem.MolToSmiles(mol)
        yield mid, smi, mol


def run_aizynth(items: list[dict], config: Path) -> dict:
    """Call aizynth_score.py in aizynth_env; return {smiles: {solved,n_steps,...}}."""
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fin:
        json.dump(items, fin)
        in_path = fin.name
    out_path = in_path.replace(".json", "_out.json")
    cmd = ["conda", "run", "--no-capture-output", "-n", "aizynth_env", "python",
           str(AIZYNTH_SCORER), "--config", str(config),
           "--input", in_path, "--output", out_path]
    subprocess.run(cmd, check=True)
    return json.load(open(out_path))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="SDF of candidate molecules")
    ap.add_argument("--config", required=True, help="AiZynthFinder config.yml")
    ap.add_argument("--output", required=True, help="SDF to write makeable survivors")
    ap.add_argument("--report", required=True, help="JSON gate report path")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--require-solved", action="store_true",
                   help="Keep only molecules with a complete route (default)")
    g.add_argument("--min-score", type=float, default=None,
                   help="Instead, keep molecules with synthesis_score >= this")
    ap.add_argument("--proxy", choices=["rascore", "none"], default="rascore",
                    help="Fast pre-filter before AiZynth (default: rascore if "
                         "installed, else auto-fallback to full AiZynth)")
    ap.add_argument("--proxy-threshold", type=float, default=0.5,
                    help="Keep molecules with proxy P(solvable) >= this for the "
                         "AiZynth confirm step")
    ap.add_argument("--no-confirm", action="store_true",
                    help="Trust the proxy alone; skip the AiZynth route search "
                         "(fastest, but no real route is produced)")
    args = ap.parse_args()

    from rdkit import Chem

    records = list(_read_sdf(Path(args.input)))
    if not records:
        print(f"[gate] no molecules in {args.input}", file=sys.stderr)
        return 1
    items = [{"ligand_id": mid or f"mol_{i:04d}", "smiles": smi}
             for i, (mid, smi, _) in enumerate(records)]

    # Fast tier: proxy pre-filter. Reduces the expensive AiZynth set to the
    # molecules the proxy thinks are makeable. Falls back gracefully when the
    # proxy isn't installed (proxy_scores empty → AiZynth on everything).
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    import synth_proxy
    proxy_scores: dict[str, float] = {}
    if args.proxy == "rascore" and synth_proxy.available():
        proxy_scores = synth_proxy.score([smi for _, smi, _ in records])
        kept_smi = {s for s, v in proxy_scores.items() if v >= args.proxy_threshold}
        items = [it for it in items if it["smiles"] in kept_smi]
        print(f"[gate] proxy (RAScore) kept {len(items)}/{len(records)} "
              f"(P>= {args.proxy_threshold}) for the AiZynth confirm step")
    else:
        print("[gate] no fast proxy available — running full AiZynth on all "
              "molecules (install RAScore to enable the fast path)")

    if args.no_confirm and proxy_scores:
        # Trust the proxy: synthesize a pseudo-result so downstream is uniform.
        scores = {it["smiles"]: {"solved": True, "n_steps": -1,
                                 "synthesis_score": round(proxy_scores[it["smiles"]], 3)}
                  for it in items}
    elif items:
        print(f"[gate] {len(items)} molecules → AiZynth retrosynthesis…")
        scores = run_aizynth(items, Path(args.config))
    else:
        scores = {}

    def passes(smi: str) -> bool:
        s = scores.get(smi, {})
        if args.min_score is not None:
            return s.get("synthesis_score", 0.0) >= args.min_score
        return bool(s.get("solved", False))

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(out_path))
    kept = []
    for mid, smi, mol in records:
        s = scores.get(smi, {})
        if passes(smi):
            mol.SetProp("synthesis_solved", str(s.get("solved", False)))
            mol.SetProp("synthesis_steps", str(s.get("n_steps", -1)))
            mol.SetProp("synthesis_score", str(s.get("synthesis_score", 0.0)))
            writer.write(mol)
            kept.append({"molecule_id": mid, "smiles": smi, **s})
    writer.close()

    report = {
        "input": str(args.input),
        "criterion": "min_score>=%.2f" % args.min_score if args.min_score is not None
                     else "require_solved",
        "total_in": len(records),
        "total_kept": len(kept),
        "pass_rate": round(100.0 * len(kept) / len(records), 1),
        "kept": kept,
    }
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    json.dump(report, open(args.report, "w"), indent=2)
    print(f"[gate] kept {len(kept)}/{len(records)} makeable "
          f"({report['pass_rate']}%) → {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
