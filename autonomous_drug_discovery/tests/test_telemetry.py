"""
Unit tests for the TelemetryDB module.

Run: python -m pytest tests/test_telemetry.py -v
"""

import os
import sys
import tempfile
import json
import sqlite3
import pytest
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.resolve()))

from telemetry import TelemetryDB


@pytest.fixture
def db_path(tmp_path):
    """Provide a temporary database path."""
    return str(tmp_path / "test_telemetry.db")


@pytest.fixture
def db(db_path):
    """Provide a fresh TelemetryDB instance."""
    database = TelemetryDB(db_path)
    yield database
    database.close()


class TestSchemaInitialization:
    """Tests for database creation and schema."""

    def test_creates_database_file(self, db_path):
        db = TelemetryDB(db_path)
        assert Path(db_path).exists()
        db.close()

    def test_creates_tables(self, db_path):
        db = TelemetryDB(db_path)
        conn = sqlite3.connect(db_path)
        tables = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()
        table_names = {t[0] for t in tables}
        assert "runs" in table_names
        assert "molecule_scores" in table_names
        conn.close()
        db.close()

    def test_creates_parent_directories(self, tmp_path):
        deep_path = str(tmp_path / "a" / "b" / "c" / "test.db")
        db = TelemetryDB(deep_path)
        assert Path(deep_path).exists()
        db.close()

    def test_idempotent_schema(self, db_path):
        """Creating TelemetryDB twice on the same file should not error."""
        db1 = TelemetryDB(db_path)
        db1.close()
        db2 = TelemetryDB(db_path)
        db2.close()


class TestRunLifecycle:
    """Tests for start_run / complete_run round-trip."""

    def test_start_run_returns_uuid(self, db):
        run_id = db.start_run("campaign_1", "01_ingestion", "/input.pdb", {"key": "val"})
        assert isinstance(run_id, str)
        assert len(run_id) == 36  # UUID format

    def test_start_run_creates_record(self, db):
        run_id = db.start_run("campaign_1", "01_ingestion", "/input.pdb", {"key": "val"})
        runs = db.query_runs()
        assert len(runs) == 1
        assert runs[0]["run_id"] == run_id
        assert runs[0]["status"] == "running"
        assert runs[0]["campaign_id"] == "campaign_1"
        assert runs[0]["module_name"] == "01_ingestion"

    def test_complete_run_success(self, db):
        run_id = db.start_run("campaign_1", "01_ingestion", "/input.pdb", {})
        db.complete_run(run_id, "success", "/output.json")
        runs = db.query_runs()
        assert runs[0]["status"] == "success"
        assert runs[0]["output_path"] == "/output.json"
        assert runs[0]["completed_at"] is not None

    def test_complete_run_failure(self, db):
        run_id = db.start_run("campaign_1", "02_generation", "/input.json", {})
        db.complete_run(run_id, "failed", error_trace="Traceback: boom!")
        runs = db.query_runs()
        assert runs[0]["status"] == "failed"
        assert "boom" in runs[0]["error_trace"]

    def test_parameters_stored_as_json(self, db):
        params = {"temperature": 0.8, "steps": 1000, "schedule": "cosine"}
        run_id = db.start_run("campaign_1", "02_generation", "/in.json", params)
        runs = db.query_runs()
        stored_params = json.loads(runs[0]["parameters"])
        assert stored_params == params

    def test_git_commit_stored(self, db):
        run_id = db.start_run(
            "campaign_1", "02_generation", "/in.json", {},
            git_commit="abc123def456"
        )
        runs = db.query_runs()
        assert runs[0]["git_commit"] == "abc123def456"


class TestMoleculeLogging:
    """Tests for per-molecule scoring data."""

    def test_log_molecule(self, db):
        run_id = db.start_run("c1", "03_screening", "/in.sdf", {})
        db.log_molecule(run_id, "mol_001", "CCO", {
            "qed": 0.75,
            "sa_score": 2.1,
            "logp": 0.5,
            "mol_weight": 46.07,
            "passed_triage": 1,
        })
        molecules = db.query_molecules(run_id)
        assert len(molecules) == 1
        assert molecules[0]["smiles"] == "CCO"
        assert molecules[0]["qed"] == 0.75
        assert molecules[0]["passed_triage"] == 1

    def test_log_molecule_rejected(self, db):
        run_id = db.start_run("c1", "03_screening", "/in.sdf", {})
        db.log_molecule(run_id, "mol_002", "CCC(=O)Oc1ccccc1", {
            "passed_triage": 0,
            "stage_eliminated": "pains_hit",
        })
        molecules = db.query_molecules(run_id)
        assert molecules[0]["passed_triage"] == 0
        assert molecules[0]["stage_eliminated"] == "pains_hit"

    def test_log_molecules_batch(self, db):
        run_id = db.start_run("c1", "03_screening", "/in.sdf", {})
        batch = [
            {"molecule_id": f"mol_{i:04d}", "smiles": "CCO", "qed": 0.5 + i * 0.1}
            for i in range(5)
        ]
        db.log_molecules_batch(run_id, batch)
        molecules = db.query_molecules(run_id)
        assert len(molecules) == 5

    def test_log_molecule_partial_scores(self, db):
        """Molecule with only some score fields populated."""
        run_id = db.start_run("c1", "04_docking", "/in.json", {})
        db.log_molecule(run_id, "mol_001", None, {"docking_score": -9.5})
        molecules = db.query_molecules(run_id)
        assert molecules[0]["docking_score"] == -9.5
        assert molecules[0]["smiles"] is None
        assert molecules[0]["qed"] is None


class TestQueryFiltering:
    """Tests for query_runs filtering."""

    def test_filter_by_campaign(self, db):
        db.start_run("campaign_A", "01_ingestion", "/a.pdb", {})
        db.start_run("campaign_B", "01_ingestion", "/b.pdb", {})
        runs = db.query_runs(campaign_id="campaign_A")
        assert len(runs) == 1
        assert runs[0]["campaign_id"] == "campaign_A"

    def test_filter_by_module(self, db):
        db.start_run("c1", "01_ingestion", "/a.pdb", {})
        db.start_run("c1", "02_generation", "/a.json", {})
        runs = db.query_runs(module_name="02_generation")
        assert len(runs) == 1

    def test_filter_by_status(self, db):
        r1 = db.start_run("c1", "01_ingestion", "/a.pdb", {})
        r2 = db.start_run("c1", "02_generation", "/a.json", {})
        db.complete_run(r1, "success", "/out.json")
        db.complete_run(r2, "failed", error_trace="error")
        runs = db.query_runs(status="failed")
        assert len(runs) == 1
        assert runs[0]["module_name"] == "02_generation"


class TestCampaignSummary:
    """Tests for campaign-level summaries."""

    def test_campaign_summary(self, db):
        r1 = db.start_run("c1", "01_ingestion", "/a.pdb", {})
        db.complete_run(r1, "success", "/out.json")
        r2 = db.start_run("c1", "03_screening", "/a.sdf", {})
        db.complete_run(r2, "success", "/out.sdf")
        db.log_molecule(r2, "mol_001", "CCO", {"passed_triage": 1})
        db.log_molecule(r2, "mol_002", "CCC", {"passed_triage": 0, "stage_eliminated": "qed"})

        summary = db.get_campaign_summary("c1")
        assert summary["campaign_id"] == "c1"
        assert summary["runs_by_status"]["success"] == 2
        assert summary["total_molecules"] == 2


class TestContextManager:
    """Tests for context manager usage."""

    def test_context_manager(self, db_path):
        with TelemetryDB(db_path) as db:
            run_id = db.start_run("c1", "01_ingestion", "/a.pdb", {})
            db.complete_run(run_id, "success", "/out.json")
        # Verify data persisted after close
        db2 = TelemetryDB(db_path)
        runs = db2.query_runs()
        assert len(runs) == 1
        db2.close()
