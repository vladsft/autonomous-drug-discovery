"""Tests for scripts/merge_telemetry.py — telemetry DB union.

This script can silently lose campaign history if the merge logic is wrong, so
we exercise (a) a successful additive merge, (b) the run_id collision abort,
and (c) FK integrity of the merged database.
"""

from __future__ import annotations

import json
import sqlite3
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
MERGE_SCRIPT = REPO_ROOT / "scripts" / "merge_telemetry.py"
sys.path.insert(0, str(REPO_ROOT / "autonomous_drug_discovery"))

from telemetry import TelemetryDB  # noqa: E402


def _seed_db(path: Path, n_runs: int = 2, mol_per_run: int = 3,
             run_id_prefix: str = "") -> list[str]:
    """Drop a tiny telemetry DB with realistic schema + content. Returns run IDs."""
    db = TelemetryDB(str(path))
    run_ids: list[str] = []
    for i in range(n_runs):
        rid = f"{run_id_prefix}{uuid.uuid4()}"
        run_ids.append(rid)
        db._conn.execute(
            "INSERT INTO runs (run_id, campaign_id, module_name, started_at, "
            "status, parameters) VALUES (?, ?, ?, ?, ?, ?)",
            (rid, f"campaign_{i}", "02_generation",
             datetime.now(timezone.utc).isoformat(), "success",
             json.dumps({"mode": "rdkit"})),
        )
        db.log_molecules_batch(rid, [
            {"molecule_id": f"mol_{i}_{j}", "smiles": "CCO", "qed": 0.5}
            for j in range(mol_per_run)
        ])
    db.close()
    return run_ids


def _run_merge(base: Path, incoming: Path, out: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(MERGE_SCRIPT), str(base), str(incoming),
         "--output", str(out)],
        capture_output=True, text=True,
    )


def test_disjoint_merge_unions_rows(tmp_path):
    base, incoming, out = tmp_path / "base.db", tmp_path / "inc.db", tmp_path / "out.db"
    _seed_db(base, n_runs=2, mol_per_run=3, run_id_prefix="base-")
    _seed_db(incoming, n_runs=3, mol_per_run=2, run_id_prefix="inc-")

    proc = _run_merge(base, incoming, out)
    assert proc.returncode == 0, proc.stderr

    conn = sqlite3.connect(str(out))
    runs = conn.execute("SELECT COUNT(*) FROM runs").fetchone()[0]
    mols = conn.execute("SELECT COUNT(*) FROM molecule_scores").fetchone()[0]
    conn.close()
    assert runs == 5            # 2 base + 3 incoming
    assert mols == 2 * 3 + 3 * 2  # base mols + incoming mols


def test_collision_aborts_without_overwrite(tmp_path):
    base, incoming, out = tmp_path / "base.db", tmp_path / "inc.db", tmp_path / "out.db"
    shared_id = f"shared-{uuid.uuid4()}"
    _seed_db(base, run_id_prefix=shared_id + "-base-")
    # Force a run_id collision: seed incoming with the same uuid as a base row.
    db_base = sqlite3.connect(str(base))
    base_rid = db_base.execute("SELECT run_id FROM runs LIMIT 1").fetchone()[0]
    db_base.close()

    _seed_db(incoming, run_id_prefix="inc-")
    # Insert a colliding run_id into the incoming DB.
    db_inc = sqlite3.connect(str(incoming))
    db_inc.execute(
        "INSERT INTO runs (run_id, campaign_id, module_name, started_at, status) "
        "VALUES (?, 'c', 'm', ?, 's')",
        (base_rid, datetime.now(timezone.utc).isoformat()),
    )
    db_inc.commit()
    db_inc.close()

    proc = _run_merge(base, incoming, out)
    assert proc.returncode == 1, "expected non-zero on collision"
    assert "exist in both DBs" in proc.stderr or "refusing to merge" in proc.stderr


def test_foreign_keys_intact_after_merge(tmp_path):
    """Every molecule_scores.run_id must resolve to a run in the merged DB."""
    base, incoming, out = tmp_path / "base.db", tmp_path / "inc.db", tmp_path / "out.db"
    _seed_db(base, n_runs=1, mol_per_run=5, run_id_prefix="base-")
    _seed_db(incoming, n_runs=1, mol_per_run=4, run_id_prefix="inc-")

    proc = _run_merge(base, incoming, out)
    assert proc.returncode == 0, proc.stderr

    conn = sqlite3.connect(str(out))
    orphans = conn.execute(
        "SELECT COUNT(*) FROM molecule_scores "
        "WHERE run_id NOT IN (SELECT run_id FROM runs)"
    ).fetchone()[0]
    conn.close()
    assert orphans == 0


def test_missing_input_returns_nonzero(tmp_path):
    proc = _run_merge(tmp_path / "missing.db", tmp_path / "also-missing.db",
                      tmp_path / "out.db")
    assert proc.returncode == 1
