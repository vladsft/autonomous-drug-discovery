"""
Telemetry Module for the Adaptive Discovery Orchestrator.

Provides a lightweight SQLite-backed database for capturing run metadata,
generation parameters, and per-molecule scoring distributions. This is the
foundation of the compounding data advantage described in north_star.md.

Usage:
    db = TelemetryDB("path/to/telemetry.db")
    run_id = db.start_run("campaign_001", "01_ingestion", "/data/target.pdb", {"param": "value"})
    db.complete_run(run_id, "success", "/data/output.json")
    db.log_molecule(run_id, "mol_001", "CCO", {"qed": 0.8, "sa_score": 2.1})
"""

import sqlite3
import uuid
import json
import hashlib
from datetime import datetime, timezone
from pathlib import Path


_SCHEMA_RUNS = """
CREATE TABLE IF NOT EXISTS runs (
    run_id          TEXT PRIMARY KEY,
    campaign_id     TEXT NOT NULL,
    module_name     TEXT NOT NULL,
    started_at      TEXT NOT NULL,
    completed_at    TEXT,
    status          TEXT NOT NULL,
    input_hash      TEXT,
    input_path      TEXT,
    output_path     TEXT,
    parameters      TEXT,
    error_trace     TEXT,
    git_commit      TEXT,
    notes           TEXT
);
"""

_SCHEMA_MOLECULE_SCORES = """
CREATE TABLE IF NOT EXISTS molecule_scores (
    score_id        INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          TEXT NOT NULL REFERENCES runs(run_id),
    molecule_id     TEXT NOT NULL,
    smiles          TEXT,
    qed             REAL,
    sa_score        REAL,
    logp            REAL,
    mol_weight      REAL,
    docking_score   REAL,
    passed_triage   INTEGER,
    stage_eliminated TEXT
);
"""

_SCHEMA_INDEXES = """
CREATE INDEX IF NOT EXISTS idx_runs_campaign ON runs(campaign_id);
CREATE INDEX IF NOT EXISTS idx_runs_module ON runs(module_name);
CREATE INDEX IF NOT EXISTS idx_runs_status ON runs(status);
CREATE INDEX IF NOT EXISTS idx_molscores_run ON molecule_scores(run_id);
CREATE INDEX IF NOT EXISTS idx_molscores_triage ON molecule_scores(passed_triage);
"""


def _compute_file_hash(filepath: str) -> str | None:
    """Compute SHA256 hash of a file. Returns None if file doesn't exist."""
    try:
        path = Path(filepath)
        if not path.exists():
            return None
        sha256 = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                sha256.update(chunk)
        return sha256.hexdigest()
    except (OSError, TypeError):
        return None


class TelemetryDB:
    """SQLite-backed telemetry database for the Discovery Orchestrator.

    Thread-safety: Each TelemetryDB instance holds its own connection.
    For multi-threaded use, create one instance per thread.
    """

    def __init__(self, db_path: str):
        """Initialize the telemetry database.

        Creates the database file and schema tables if they do not exist.

        Args:
            db_path: Path to the SQLite database file.
        """
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(self.db_path))
        self._conn.row_factory = sqlite3.Row
        self._conn.execute("PRAGMA journal_mode=WAL;")
        self._conn.execute("PRAGMA foreign_keys=ON;")
        self._initialize_schema()

    def _initialize_schema(self):
        """Create tables and indexes if they don't exist."""
        cursor = self._conn.cursor()
        cursor.executescript(_SCHEMA_RUNS + _SCHEMA_MOLECULE_SCORES + _SCHEMA_INDEXES)
        self._conn.commit()

    def start_run(
        self,
        campaign_id: str,
        module_name: str,
        input_path: str,
        parameters: dict,
        git_commit: str | None = None,
        notes: str | None = None,
    ) -> str:
        """Record the start of a module execution.

        Args:
            campaign_id: Groups runs under one target campaign.
            module_name: e.g., "01_ingestion", "02_generation".
            input_path: Path to the primary input file.
            parameters: Dict of all parameters used for this run.
            git_commit: Optional git commit hash of the model repo.
            notes: Optional free-text annotations.

        Returns:
            run_id: A UUID string identifying this run.
        """
        run_id = str(uuid.uuid4())
        input_hash = _compute_file_hash(input_path)
        now = datetime.now(timezone.utc).isoformat()

        self._conn.execute(
            """INSERT INTO runs
               (run_id, campaign_id, module_name, started_at, status,
                input_hash, input_path, parameters, git_commit, notes)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                run_id,
                campaign_id,
                module_name,
                now,
                "running",
                input_hash,
                str(input_path),
                json.dumps(parameters),
                git_commit,
                notes,
            ),
        )
        self._conn.commit()
        return run_id

    def complete_run(
        self,
        run_id: str,
        status: str,
        output_path: str | None = None,
        error_trace: str | None = None,
    ):
        """Record the completion (or failure) of a module execution.

        Args:
            run_id: The UUID returned by start_run.
            status: "success" or "failed".
            output_path: Path to the primary output file.
            error_trace: Full traceback string on failure.
        """
        now = datetime.now(timezone.utc).isoformat()
        self._conn.execute(
            """UPDATE runs
               SET completed_at = ?, status = ?, output_path = ?, error_trace = ?
               WHERE run_id = ?""",
            (now, status, str(output_path) if output_path else None, error_trace, run_id),
        )
        self._conn.commit()

    def log_molecule(
        self,
        run_id: str,
        molecule_id: str,
        smiles: str | None,
        scores: dict,
    ):
        """Log per-molecule scoring data.

        Args:
            run_id: The UUID of the run that produced this molecule.
            molecule_id: Internal identifier (e.g., "mol_001").
            smiles: Canonical SMILES string, or None.
            scores: Dict with optional keys: qed, sa_score, logp, mol_weight,
                    docking_score, passed_triage, stage_eliminated.
        """
        self._conn.execute(
            """INSERT INTO molecule_scores
               (run_id, molecule_id, smiles, qed, sa_score, logp,
                mol_weight, docking_score, passed_triage, stage_eliminated)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                run_id,
                molecule_id,
                smiles,
                scores.get("qed"),
                scores.get("sa_score"),
                scores.get("logp"),
                scores.get("mol_weight"),
                scores.get("docking_score"),
                scores.get("passed_triage"),
                scores.get("stage_eliminated"),
            ),
        )
        self._conn.commit()

    def log_molecules_batch(self, run_id: str, molecules: list[dict]):
        """Log multiple molecules in a single transaction for performance.

        Args:
            run_id: The UUID of the run.
            molecules: List of dicts, each with keys: molecule_id, smiles, and
                       optional score keys (qed, sa_score, logp, etc.).
        """
        rows = []
        for mol in molecules:
            rows.append((
                run_id,
                mol["molecule_id"],
                mol.get("smiles"),
                mol.get("qed"),
                mol.get("sa_score"),
                mol.get("logp"),
                mol.get("mol_weight"),
                mol.get("docking_score"),
                mol.get("passed_triage"),
                mol.get("stage_eliminated"),
            ))
        self._conn.executemany(
            """INSERT INTO molecule_scores
               (run_id, molecule_id, smiles, qed, sa_score, logp,
                mol_weight, docking_score, passed_triage, stage_eliminated)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            rows,
        )
        self._conn.commit()

    def query_runs(
        self,
        campaign_id: str | None = None,
        module_name: str | None = None,
        status: str | None = None,
    ) -> list[dict]:
        """Query run records with optional filters.

        Args:
            campaign_id: Filter by campaign. None = all campaigns.
            module_name: Filter by module. None = all modules.
            status: Filter by status. None = all statuses.

        Returns:
            List of run records as dicts.
        """
        query = "SELECT * FROM runs WHERE 1=1"
        params = []

        if campaign_id is not None:
            query += " AND campaign_id = ?"
            params.append(campaign_id)
        if module_name is not None:
            query += " AND module_name = ?"
            params.append(module_name)
        if status is not None:
            query += " AND status = ?"
            params.append(status)

        query += " ORDER BY started_at DESC"

        rows = self._conn.execute(query, params).fetchall()
        return [dict(row) for row in rows]

    def query_molecules(self, run_id: str) -> list[dict]:
        """Retrieve all molecule scores for a given run.

        Args:
            run_id: The UUID of the run.

        Returns:
            List of molecule score records as dicts.
        """
        rows = self._conn.execute(
            "SELECT * FROM molecule_scores WHERE run_id = ? ORDER BY score_id",
            (run_id,),
        ).fetchall()
        return [dict(row) for row in rows]

    def get_campaign_summary(self, campaign_id: str) -> dict:
        """Get a high-level summary of a campaign's progress.

        Args:
            campaign_id: The campaign to summarize.

        Returns:
            Dict with run counts by status and total molecules logged.
        """
        run_counts = self._conn.execute(
            "SELECT status, COUNT(*) as count FROM runs WHERE campaign_id = ? GROUP BY status",
            (campaign_id,),
        ).fetchall()

        mol_count = self._conn.execute(
            """SELECT COUNT(*) as count FROM molecule_scores ms
               JOIN runs r ON ms.run_id = r.run_id
               WHERE r.campaign_id = ?""",
            (campaign_id,),
        ).fetchone()

        triage_stats = self._conn.execute(
            """SELECT passed_triage, COUNT(*) as count FROM molecule_scores ms
               JOIN runs r ON ms.run_id = r.run_id
               WHERE r.campaign_id = ?
               GROUP BY passed_triage""",
            (campaign_id,),
        ).fetchall()

        return {
            "campaign_id": campaign_id,
            "runs_by_status": {row["status"]: row["count"] for row in run_counts},
            "total_molecules": mol_count["count"] if mol_count else 0,
            "triage_stats": {
                str(row["passed_triage"]): row["count"] for row in triage_stats
            },
        }

    def close(self):
        """Close the database connection."""
        if self._conn:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
