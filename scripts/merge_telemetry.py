"""Merge two telemetry databases into one.

Both the M1/M2 dev-box history and the Hugging Face mirror hold campaign
records, but on disjoint campaign sets. This unions them so no history is
lost. `runs.run_id` is a text UUID primary key, so runs merge directly;
`molecule_scores.score_id` is an autoincrement integer, so those rows are
re-inserted and the merged DB assigns fresh IDs (`run_id` is text, FK intact).

The merge is additive and idempotent-safe: it aborts if a run_id already
exists in the base, rather than silently overwriting.

Usage:
    python scripts/merge_telemetry.py BASE.db INCOMING.db --output MERGED.db
"""

from __future__ import annotations

import argparse
import shutil
import sqlite3
import sys
from pathlib import Path

RUN_COLS = (
    "run_id", "campaign_id", "module_name", "started_at", "completed_at",
    "status", "input_hash", "input_path", "output_path", "parameters",
    "error_trace", "git_commit", "notes",
)
SCORE_COLS = (
    "run_id", "molecule_id", "smiles", "qed", "sa_score", "logp",
    "mol_weight", "docking_score", "passed_triage", "stage_eliminated",
)


def counts(db: sqlite3.Connection) -> dict:
    return {
        "campaigns": db.execute("SELECT COUNT(DISTINCT campaign_id) FROM runs").fetchone()[0],
        "runs": db.execute("SELECT COUNT(*) FROM runs").fetchone()[0],
        "molecule_scores": db.execute("SELECT COUNT(*) FROM molecule_scores").fetchone()[0],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("base", help="Base DB — kept intact, copied to --output")
    ap.add_argument("incoming", help="DB whose rows are merged into the base")
    ap.add_argument("--output", required=True, help="Path for the merged DB")
    args = ap.parse_args()

    base, incoming, out = Path(args.base), Path(args.incoming), Path(args.output)
    for p in (base, incoming):
        if not p.exists():
            print(f"ERROR: {p} not found", file=sys.stderr)
            return 1

    shutil.copyfile(base, out)
    db = sqlite3.connect(out)
    src = sqlite3.connect(incoming)

    print(f"base     {base.name}: {counts(db)}")
    print(f"incoming {incoming.name}: {counts(src)}")

    base_run_ids = {r[0] for r in db.execute("SELECT run_id FROM runs")}
    incoming_runs = src.execute(f"SELECT {','.join(RUN_COLS)} FROM runs").fetchall()

    collisions = [r[0] for r in incoming_runs if r[0] in base_run_ids]
    if collisions:
        print(f"ERROR: {len(collisions)} run_id(s) exist in both DBs, e.g. "
              f"{collisions[:3]} — refusing to merge.", file=sys.stderr)
        return 1

    db.executemany(
        f"INSERT INTO runs ({','.join(RUN_COLS)}) "
        f"VALUES ({','.join('?' * len(RUN_COLS))})",
        incoming_runs,
    )
    # score_id is omitted so the merged DB autoincrements it afresh.
    incoming_scores = src.execute(f"SELECT {','.join(SCORE_COLS)} FROM molecule_scores").fetchall()
    db.executemany(
        f"INSERT INTO molecule_scores ({','.join(SCORE_COLS)}) "
        f"VALUES ({','.join('?' * len(SCORE_COLS))})",
        incoming_scores,
    )
    db.commit()

    merged = counts(db)
    print(f"merged   {out.name}: {merged}")

    # Integrity: every molecule_scores.run_id must resolve to a run.
    orphans = db.execute(
        "SELECT COUNT(*) FROM molecule_scores "
        "WHERE run_id NOT IN (SELECT run_id FROM runs)"
    ).fetchone()[0]
    db.close()
    src.close()
    if orphans:
        print(f"ERROR: {orphans} orphaned molecule_scores rows after merge", file=sys.stderr)
        return 1
    print("OK — foreign keys intact.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
