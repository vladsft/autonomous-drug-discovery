"""
Build an ingestion manifest for a specific P2Rank pocket rank.

The orchestrator's ingestion always picks P2Rank pocket #1 (highest score).
This standalone helper lets a campaign target a different pocket — e.g. a
functionally relevant site that out-scores on residue chemistry rather than on
P2Rank's geometric ligandability — without modifying the auto-selection logic.

It reuses an existing P2Rank run: it does not re-run P2Rank, only re-reads its
predictions CSV and writes a manifest + pocket PDB for the chosen rank.

Usage:
    python scripts/select_pocket.py <processed.pdb> --pocket-rank 2
"""

import csv
import json
import argparse
from pathlib import Path


def parse_predictions(predictions_csv):
    """Return P2Rank pockets as a list of dicts, sorted by rank."""
    pockets = []
    with open(predictions_csv, newline="") as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        reader.fieldnames = [h.strip() for h in (reader.fieldnames or [])]
        for raw in reader:
            row = {(k or "").strip(): (v.strip() if isinstance(v, str) else v)
                   for k, v in raw.items()}
            pockets.append({
                "name": row["name"],
                "rank": int(row["rank"]),
                "score": float(row["score"]),
                "probability": float(row["probability"]),
                "center": [float(row["center_x"]), float(row["center_y"]),
                           float(row["center_z"])],
                "residue_ids": row["residue_ids"],
            })
    return sorted(pockets, key=lambda p: p["rank"])


def write_pocket_pdb(source_pdb, residue_ids, out_pdb):
    """Extract the atoms of the pocket's residues from `source_pdb`."""
    wanted = set()
    for rid in residue_ids.split():
        parts = rid.split("_")
        if len(parts) == 2:
            try:
                wanted.add((parts[0], int(parts[1])))
            except ValueError:
                continue
    n = 0
    with open(source_pdb) as fin, open(out_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    continue
                if (chain, resnum) in wanted:
                    fout.write(line)
                    n += 1
    if n == 0:
        raise RuntimeError(f"No atoms matched pocket residues in {source_pdb}.")
    return n


def main():
    ap = argparse.ArgumentParser(description="Build a manifest for a P2Rank pocket rank")
    ap.add_argument("pdb", help="Processed PDB that was run through P2Rank")
    ap.add_argument("--pocket-rank", type=int, required=True, help="P2Rank pocket rank")
    ap.add_argument("--p2rank-dir", default=None,
                    help="P2Rank output dir (default: <pdb-dir>/<stem>_p2rank)")
    ap.add_argument("--output", default=None,
                    help="Manifest path (default: <pdb-dir>/<stem>_pocket<rank>_manifest.json)")
    args = ap.parse_args()

    pdb = Path(args.pdb).resolve()
    p2rank_dir = Path(args.p2rank_dir) if args.p2rank_dir else pdb.parent / f"{pdb.stem}_p2rank"
    predictions = p2rank_dir / f"{pdb.name}_predictions.csv"
    if not predictions.exists():
        raise FileNotFoundError(f"P2Rank predictions not found: {predictions}")

    pockets = parse_predictions(predictions)
    chosen = next((p for p in pockets if p["rank"] == args.pocket_rank), None)
    if chosen is None:
        raise ValueError(f"No pocket with rank {args.pocket_rank} "
                         f"(found ranks 1-{len(pockets)}).")

    pocket_pdb = p2rank_dir / f"{pdb.stem}_pocket{args.pocket_rank}_atm.pdb"
    n_atoms = write_pocket_pdb(pdb, chosen["residue_ids"], pocket_pdb)

    manifest = {
        "input_pdb": str(pdb),
        "run_pdb": str(pdb),
        "receptor_pdb": str(pdb),
        "pocket_backend": "p2rank",
        "p2rank_out_dir": str(p2rank_dir),
        "pockets_found": len(pockets),
        "selected_pocket_rank": args.pocket_rank,
        "best_pocket": str(pocket_pdb),
        "best_pocket_score": chosen["score"],
        "best_pocket_probability": chosen["probability"],
        "best_pocket_center": chosen["center"],
    }
    out = Path(args.output) if args.output else \
        pdb.parent / f"{pdb.stem}_pocket{args.pocket_rank}_manifest.json"
    out.write_text(json.dumps(manifest, indent=4))
    print(f"Pocket {args.pocket_rank}: score {chosen['score']:.2f}, "
          f"prob {chosen['probability']:.3f}, {n_atoms} atoms → {pocket_pdb.name}")
    print(f"Manifest: {out}")


if __name__ == "__main__":
    main()
