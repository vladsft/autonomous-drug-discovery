"""Score synthetic accessibility with AiZynthFinder. Runs inside `aizynth_env`.

Stage 5 ranking (modules/05_ranking/run_ranking.py) invokes this as a
subprocess — `conda run -n aizynth_env python scripts/aizynth_score.py` — so
AiZynthFinder's ML stack never has to co-resolve with the orchestrator env.

Contract:
  --input   JSON list of {"ligand_id": ..., "smiles": ...}
  --config  AiZynthFinder config.yml (from `download_public_data`)
  --output  JSON written as {smiles: {solved, n_steps, top_score,
            synthesis_score}}

For each SMILES it runs a retrosynthetic tree search: can this molecule be
made from purchasable building blocks, and in how many steps?

synthesis_score is in [0, 1], higher = easier to synthesise:
  - route found:  max(0.30, 1.0 - 0.12 * (n_steps - 1))   fewer steps is better
  - no route:     0.15                                     a soft penalty, not 0
"""

import argparse
import json
import sys
from pathlib import Path


def synthesis_score(solved: bool, n_steps: int) -> float:
    """Map an AiZynthFinder result to a [0, 1] feasibility score."""
    if not solved:
        return 0.15
    steps = max(int(n_steps), 1)
    return round(max(0.30, min(1.0, 1.0 - 0.12 * (steps - 1))), 4)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True, help="AiZynthFinder config.yml")
    ap.add_argument("--input", required=True, help="JSON list of {ligand_id, smiles}")
    ap.add_argument("--output", required=True, help="Path for the results JSON")
    args = ap.parse_args()

    cfg = Path(args.config).expanduser()
    if not cfg.exists():
        print(f"[aizynth_score] config not found: {cfg}", file=sys.stderr)
        return 1

    try:
        from aizynthfinder.aizynthfinder import AiZynthFinder
    except ImportError as e:
        print(f"[aizynth_score] aizynthfinder not importable: {e}", file=sys.stderr)
        return 1

    candidates = json.loads(Path(args.input).read_text())

    finder = AiZynthFinder(configfile=str(cfg))
    # The config lists available stocks/policies; the API still needs an
    # explicit selection. Take the first of each — the public dataset ships
    # exactly one expansion policy and one stock.
    finder.stock.select(finder.stock.items[0])
    finder.expansion_policy.select(finder.expansion_policy.items[0])
    if finder.filter_policy.items:
        finder.filter_policy.select(finder.filter_policy.items[0])

    out: dict[str, dict] = {}
    for cand in candidates:
        smi = cand.get("smiles")
        if not smi or smi in out:
            continue
        try:
            finder.target_smiles = smi
            finder.tree_search()
            finder.build_routes()
            stats = finder.extract_statistics()
            solved = bool(stats.get("is_solved"))
            n_steps = int(stats.get("number_of_steps") or 0)
            top_score = float(stats.get("top_score") or 0.0)
            out[smi] = {
                "solved": solved,
                "n_steps": n_steps,
                "top_score": round(top_score, 4),
                "synthesis_score": synthesis_score(solved, n_steps),
            }
        except Exception as e:  # one bad SMILES must not sink the batch
            out[smi] = {
                "solved": False,
                "n_steps": 0,
                "top_score": 0.0,
                "synthesis_score": 0.15,
                "error": str(e),
            }

    Path(args.output).write_text(json.dumps(out, indent=2))
    solved_n = sum(1 for v in out.values() if v.get("solved"))
    print(f"[aizynth_score] scored {len(out)} unique SMILES "
          f"({solved_n} with a route) → {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
