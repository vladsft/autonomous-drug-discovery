"""Tests for the pure-Python logic in scripts/batch_cloud_run.py.

We test only the parts that don't talk to RunPod or R2 — target parsing, the
skip-list query against a synthetic telemetry DB, and the cost estimate.
The networked paths (provision/poll/terminate) are exercised by the action
itself when it fires.
"""

from __future__ import annotations

import importlib.util
import json
import sqlite3
import sys
import uuid
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(REPO_ROOT / "autonomous_drug_discovery"))
sys.path.insert(0, str(SCRIPTS_DIR))


@pytest.fixture
def batch():
    """Import batch_cloud_run without registering it globally — keeps test
    isolation tight even if other tests load it later."""
    spec = importlib.util.spec_from_file_location(
        "batch_cloud_run", SCRIPTS_DIR / "batch_cloud_run.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _telemetry_with(path: Path, target: str, mode: str,
                    started_at: datetime, status: str = "success") -> None:
    """Seed a minimal telemetry row that matches the dispatcher's skip-list query."""
    from telemetry import TelemetryDB
    TelemetryDB(str(path)).close()
    conn = sqlite3.connect(str(path))
    conn.execute(
        "INSERT INTO runs (run_id, campaign_id, module_name, started_at, "
        "status, input_path, parameters) VALUES (?, ?, ?, ?, ?, ?, ?)",
        (str(uuid.uuid4()), f"campaign_{uuid.uuid4().hex[:8]}",
         "02_generation", started_at.isoformat(), status,
         f"/fake/data/processed/{target}_manifest.json",
         json.dumps({"mode": mode})),
    )
    conn.commit()
    conn.close()


class TestParseTargets:
    def test_whitespace_separated(self, batch):
        assert batch.parse_targets("1M17 2HYY 6P3D") == ["1M17", "2HYY", "6P3D"]

    def test_comma_separated(self, batch):
        assert batch.parse_targets("1M17,2HYY,6P3D") == ["1M17", "2HYY", "6P3D"]

    def test_mixed_separators_dedupes(self, batch):
        assert batch.parse_targets("1M17, 2HYY  1M17") == ["1M17", "2HYY"]

    def test_rejects_special_chars(self, batch):
        with pytest.raises(ValueError):
            batch.parse_targets("bad-target!")

    def test_rejects_too_short(self, batch):
        with pytest.raises(ValueError):
            batch.parse_targets("ABC")

    def test_accepts_8_char_codes(self, batch):
        # Future-proof against RCSB extending PDB codes
        assert batch.parse_targets("ABC12345") == ["ABC12345"]


class TestSkipList:
    def test_fresh_campaign_skipped(self, batch, tmp_path):
        db = tmp_path / "telemetry.db"
        _telemetry_with(db, "1M17", "targetdiff", datetime.now(timezone.utc))
        skip = batch.compute_skip_list(db, ["1M17", "2HYY"], "targetdiff")
        assert skip == {"1M17"}

    def test_old_campaign_not_skipped(self, batch, tmp_path):
        db = tmp_path / "telemetry.db"
        _telemetry_with(db, "1M17", "targetdiff",
                        datetime.now(timezone.utc) - timedelta(hours=48))
        skip = batch.compute_skip_list(db, ["1M17"], "targetdiff")
        assert skip == set()

    def test_different_mode_not_skipped(self, batch, tmp_path):
        db = tmp_path / "telemetry.db"
        _telemetry_with(db, "1M17", "rdkit", datetime.now(timezone.utc))
        skip = batch.compute_skip_list(db, ["1M17"], "targetdiff")
        assert skip == set()

    def test_failed_campaign_not_skipped(self, batch, tmp_path):
        db = tmp_path / "telemetry.db"
        _telemetry_with(db, "1M17", "targetdiff", datetime.now(timezone.utc),
                        status="failed")
        skip = batch.compute_skip_list(db, ["1M17"], "targetdiff")
        assert skip == set()

    def test_missing_db_returns_empty(self, batch, tmp_path):
        skip = batch.compute_skip_list(tmp_path / "nope.db", ["1M17"], "targetdiff")
        assert skip == set()

    def test_no_substring_false_match(self, batch, tmp_path):
        """If 1M17 has a fresh campaign, querying for 1M1 must NOT skip."""
        db = tmp_path / "telemetry.db"
        _telemetry_with(db, "1M17", "targetdiff", datetime.now(timezone.utc))
        # Note: 1M1 is 3 chars, parse_targets would reject — but compute_skip_list
        # is called with whatever the caller passes; verify the SQL is anchored.
        skip = batch.compute_skip_list(db, ["1M1X"], "targetdiff")
        assert skip == set()


class TestCostEstimate:
    def test_small_batch(self, batch):
        # 5 targets, all in 1 wave, 90 min ceiling
        # gpu_hours = 1 * 5 * 1.5 = 7.5
        # cost = 7.5 * 0.5 * 1.5 = 5.625
        est = batch.estimate_cost_usd(5, 5, 90)
        assert 5.5 < est < 5.75

    def test_full_50_target_batch(self, batch):
        # 50 targets, 5 parallel = 10 waves, 90 min ceiling
        # gpu_hours = 10 * 5 * 1.5 = 75
        # cost = 75 * 0.5 * 1.5 = 56.25
        est = batch.estimate_cost_usd(50, 5, 90)
        assert 55 < est < 58

    def test_higher_parallelism_doesnt_explode_cost(self, batch):
        # Same 50 targets at parallelism=10: 5 waves, same total GPU-hours.
        # 5 waves × 10 pods × 1.5h = 75 gpu-hours — *same* as parallelism=5.
        est5 = batch.estimate_cost_usd(50, 5, 90)
        est10 = batch.estimate_cost_usd(50, 10, 90)
        assert abs(est5 - est10) < 0.01


class TestEnvLoading:
    def test_missing_env_exits(self, batch, monkeypatch):
        # Clear every required key
        for k in ("R2_BUCKET", "RUNPOD_API_KEY", "RUNPOD_NETWORK_VOLUME_ID", "IMAGE"):
            monkeypatch.delenv(k, raising=False)
        with pytest.raises(SystemExit):
            batch.load_env()

    def test_full_env_returns_dict(self, batch, monkeypatch):
        required = {
            "R2_BUCKET": "agent-harness",
            "RCLONE_CONFIG_R2_TYPE": "s3",
            "RCLONE_CONFIG_R2_PROVIDER": "Cloudflare",
            "RCLONE_CONFIG_R2_ACCESS_KEY_ID": "key",
            "RCLONE_CONFIG_R2_SECRET_ACCESS_KEY": "sec",
            "RCLONE_CONFIG_R2_ENDPOINT": "https://x.r2.cloudflarestorage.com",
            "RUNPOD_API_KEY": "tok",
            "RUNPOD_NETWORK_VOLUME_ID": "vol",
            "IMAGE": "ghcr.io/x/y:latest",
        }
        for k, v in required.items():
            monkeypatch.setenv(k, v)
        env = batch.load_env()
        for k, v in required.items():
            assert env[k] == v
