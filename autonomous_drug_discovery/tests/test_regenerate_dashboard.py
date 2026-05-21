"""Tests for scripts/regenerate_dashboard.py — the no-arg `make dashboard` path.

We build a minimal fake campaign on disk (one screened SDF, one docking CSV)
and a telemetry row that points at it, then verify the script produces a
parseable professor_demo.js. The point is to catch the obvious failure modes
(missing telemetry, wrong campaign discovery, JSON corruption); the dataset
assembly itself is exercised by build_demo_dataset's existing path.
"""

from __future__ import annotations

import importlib.util
import json
import sqlite3
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(REPO_ROOT / "autonomous_drug_discovery"))
sys.path.insert(0, str(SCRIPTS_DIR))


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def regen():
    """Load scripts/regenerate_dashboard.py without registering it in sys.modules globally."""
    return _load_module("regenerate_dashboard", SCRIPTS_DIR / "regenerate_dashboard.py")


def _seed_campaign(data_dir: Path, campaign_id: str, target: str, backend: str) -> None:
    """Drop a minimal screening report + docking CSV that build_demo_dataset will consume."""
    campaign = data_dir / campaign_id
    (campaign / "screened").mkdir(parents=True, exist_ok=True)
    (campaign / "results").mkdir(parents=True, exist_ok=True)

    screening = {
        "backend": "molscore",
        "survivors": [
            {
                "molecule_id": "mol_0001",
                "smiles": "CCO",
                "desc_MolWt": 46.07,
                "desc_MolLogP": -0.14,
                "desc_QED": 0.42,
                "desc_SAScore": 1.5,
                "desc_TPSA": 20.2,
                "desc_NumHBD": 1,
                "desc_NumHBA": 1,
                "desc_NumRotatableBonds": 0,
                "desc_NumHeavyAtoms": 3,
                "filter_PAINS": 0,
            }
        ],
        "rejections": [],
    }
    (campaign / "screened" / "screening_report.json").write_text(json.dumps(screening))

    docking_csv = "ligand_id,affinity,smiles\nmol_0001,-7.2,CCO\n"
    (campaign / "results" / "docking_results.csv").write_text(docking_csv)


def _seed_telemetry(db_path: Path, campaign_id: str, target: str, backend: str) -> None:
    conn = sqlite3.connect(str(db_path))
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS runs (
            run_id TEXT PRIMARY KEY,
            campaign_id TEXT NOT NULL,
            module_name TEXT NOT NULL,
            started_at TEXT NOT NULL,
            completed_at TEXT,
            status TEXT NOT NULL,
            input_hash TEXT,
            input_path TEXT,
            output_path TEXT,
            parameters TEXT,
            error_trace TEXT,
            git_commit TEXT,
            notes TEXT
        );
        """
    )
    conn.execute(
        "INSERT INTO runs (run_id, campaign_id, module_name, started_at, status, "
        "input_path, parameters) VALUES (?, ?, ?, ?, ?, ?, ?)",
        (
            str(uuid.uuid4()),
            campaign_id,
            "02_generation",
            datetime.now(timezone.utc).isoformat(),
            "success",
            f"/fake/data/processed/{target}_manifest.json",
            json.dumps({"mode": backend, "num_samples": 1}),
        ),
    )
    conn.commit()
    conn.close()


def test_discovers_latest_campaign(regen, tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    _seed_campaign(data_dir, "campaign_test_rdkit", "1M17", "rdkit")
    _seed_telemetry(data_dir / "telemetry.db", "campaign_test_rdkit", "1M17", "rdkit")

    cid = regen.latest_successful_campaign(data_dir / "telemetry.db", "1M17", "rdkit")
    assert cid == "campaign_test_rdkit"


def test_skips_other_targets(regen, tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    _seed_campaign(data_dir, "campaign_test_2hyy", "2HYY", "rdkit")
    _seed_telemetry(data_dir / "telemetry.db", "campaign_test_2hyy", "2HYY", "rdkit")

    assert regen.latest_successful_campaign(data_dir / "telemetry.db", "1M17", "rdkit") is None


def test_emits_parseable_dashboard(regen, tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    out_dir = tmp_path / "dashboard"
    _seed_campaign(data_dir, "campaign_test_rdkit", "1M17", "rdkit")
    _seed_telemetry(data_dir / "telemetry.db", "campaign_test_rdkit", "1M17", "rdkit")

    rc = regen.assemble(
        target="1M17",
        data_dir=data_dir,
        out_dir=out_dir,
        max_per_backend=30,
        backends_requested=["rdkit"],
    )
    assert rc == 0

    json_path = out_dir / "professor_demo.json"
    js_path = out_dir / "professor_demo.js"
    assert json_path.exists() and js_path.exists()

    doc = json.loads(json_path.read_text())
    # Multi-target schema: even single-target assemble() writes through the
    # targets map so the dashboard JS has one consistent shape to read.
    assert doc["default_target"] == "1M17"
    assert doc["target_order"] == ["1M17"]
    assert "1M17" in doc["targets"]
    t = doc["targets"]["1M17"]
    assert t["pdb"] == "1M17"
    assert t["default_backend"] == "rdkit"
    assert "rdkit" in t["backends"]
    assert t["backends"]["rdkit"]["summary"]["total_generated"] == 1

    # The JS bridge must be a single `window.PROFESSOR_DEMO_DATA = …;` assignment.
    js_text = js_path.read_text()
    assert js_text.startswith("window.PROFESSOR_DEMO_DATA = ")
    assert js_text.rstrip().endswith(";")


def test_missing_telemetry_is_diagnosed(regen, tmp_path):
    rc = regen.assemble(
        target="1M17",
        data_dir=tmp_path / "missing",
        out_dir=tmp_path / "dashboard",
        max_per_backend=30,
        backends_requested=["rdkit"],
    )
    assert rc == 1


def test_discover_targets_enumerates_distinct(regen, tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    _seed_campaign(data_dir, "campaign_a", "1M17", "rdkit")
    _seed_campaign(data_dir, "campaign_b", "2HYY", "rdkit")
    _seed_campaign(data_dir, "campaign_c", "1M17", "targetdiff")  # dupes 1M17
    _seed_telemetry(data_dir / "telemetry.db", "campaign_a", "1M17", "rdkit")
    _seed_telemetry(data_dir / "telemetry.db", "campaign_b", "2HYY", "rdkit")
    _seed_telemetry(data_dir / "telemetry.db", "campaign_c", "1M17", "targetdiff")

    targets = regen.discover_targets(data_dir / "telemetry.db")
    assert set(targets) == {"1M17", "2HYY"}


def test_assemble_all_writes_multi_target_dashboard(regen, tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    out_dir = tmp_path / "dashboard"
    _seed_campaign(data_dir, "campaign_a", "1M17", "rdkit")
    _seed_campaign(data_dir, "campaign_b", "2HYY", "rdkit")
    _seed_telemetry(data_dir / "telemetry.db", "campaign_a", "1M17", "rdkit")
    _seed_telemetry(data_dir / "telemetry.db", "campaign_b", "2HYY", "rdkit")

    rc = regen.assemble_all(
        data_dir=data_dir,
        out_dir=out_dir,
        max_per_backend=30,
        backends_requested=["rdkit"],
    )
    assert rc == 0
    doc = json.loads((out_dir / "professor_demo.json").read_text())
    assert set(doc["targets"]) == {"1M17", "2HYY"}
    assert doc["default_target"] in doc["target_order"]
    for t in ("1M17", "2HYY"):
        assert doc["targets"][t]["pdb"] == t
        assert "rdkit" in doc["targets"][t]["backends"]
