#!/usr/bin/env python3
"""Batch driver for Phase 1.5 — fire-and-forget over 30-50 targets.

The GitHub Action (.github/workflows/batch.yml) is the only intended caller.
Local invocation works for debugging but expects the same env vars (R2 +
RunPod credentials) that the action wires in from repo Secrets.

Flow (per plan.md Phase 1.5):
  1. Validate inputs, fetch missing PDBs from RCSB.
  2. Pull telemetry.db from R2 → compute skip-list (24h idempotency).
  3. Cost-guard: refuse to start if RunPod balance < parallelism × est_hours × $0.30 × 1.5.
  4. ThreadPoolExecutor with `parallelism` workers; each worker provisions
     one RunPod pod, polls R2 sentinels, retries up to 3x on failure.
  5. After all pods drain: pull merged campaigns from R2, run merge_telemetry,
     write dashboard JSON, return the work-list status to the caller.

The Action handles dashboard regeneration + commit + Pages deploy after we
return — that needs `git` + GH token which doesn't belong in this script.

This script intentionally has zero hard dependencies beyond the stdlib +
requests (already in the conda base env). No aiohttp, no asyncio — sync code
with a ThreadPoolExecutor is plenty for an I/O-bound pool of ~10 pods.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Iterable

try:
    import requests
except ImportError:
    print("ERROR: `requests` not available — install it in the orchestrator env.",
          file=sys.stderr)
    sys.exit(1)

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "autonomous_drug_discovery" / "data"
TELEMETRY_PATH = DATA_DIR / "telemetry.db"

RUNPOD_GRAPHQL = "https://api.runpod.io/graphql"
RCSB_BASE = "https://files.rcsb.org/download"

# Same regex as fetch_pdb.py — keep them aligned.
TARGET_PATTERN = re.compile(r"^[A-Za-z0-9]{4,8}$")

# A polling tick of 60 s is generous enough to avoid hammering R2's S3 API and
# tight enough that a 30-min campaign isn't bottlenecked by polling latency.
SENTINEL_POLL_INTERVAL_S = 60
# Sentinel files older than this are stale (left over from a prior crashed run)
# and ignored to prevent a poisoned skip-list.
SENTINEL_STALENESS_S = 6 * 3600

# How long an individual campaign is allowed to keep a pod alive before the
# dispatcher cancels it. Equals the per-pod ceiling in pod_campaign.sh + 5 min
# slack for R2 push.
POD_TIMEOUT_MIN_DEFAULT = 90
# A pod that fails this many times in a row gets dropped from the work-list.
MAX_RETRIES_PER_TARGET = 3
# Crash-loop detection. If, after a startup grace period (long enough to pull
# the multi-GB image), the container has never stayed alive longer than
# CRASH_LOOP_UPTIME_S, its command is dying on startup — abort the attempt
# rather than waiting out the whole timeout.
CRASH_LOOP_GRACE_S = 8 * 60
CRASH_LOOP_UPTIME_S = 150

# GPU types tried in order when RUNPOD_GPU_TYPE isn't set, falling through on
# supply shortage. All are Ampere+ (fine for TargetDiff); for CPU-bound rdkit
# any of them just provides a host. Ordered cheap-and-common first to maximise
# the odds of landing a pod in a tight market.
DEFAULT_GPU_PRIORITY = ",".join([
    "NVIDIA GeForce RTX 3090",
    "NVIDIA RTX A5000",
    "NVIDIA RTX A4000",
    "NVIDIA RTX A4500",
    "NVIDIA GeForce RTX 4090",
    "NVIDIA L4",
])
# Idempotency window: targets with a successful campaign newer than this are
# skipped unless --force.
SKIP_WINDOW_HOURS = 24

# Default cost coefficients used by the pre-flight check. RunPod RTX 3090
# pricing is variable; pad with 1.5x.
EST_USD_PER_GPU_HOUR = 0.50
COST_SAFETY_MULTIPLIER = 1.5


# ── env loading ──────────────────────────────────────────────────────────────
def load_env() -> dict[str, str]:
    """Return required env vars, exiting with a clear message if anything's missing."""
    required = [
        "R2_BUCKET", "RCLONE_CONFIG_R2_TYPE", "RCLONE_CONFIG_R2_PROVIDER",
        "RCLONE_CONFIG_R2_ACCESS_KEY_ID", "RCLONE_CONFIG_R2_SECRET_ACCESS_KEY",
        "RCLONE_CONFIG_R2_ENDPOINT",
        "RUNPOD_API_KEY", "IMAGE",
    ]
    missing = [k for k in required if not os.environ.get(k)]
    if missing:
        print(f"ERROR: missing required env: {' '.join(missing)}", file=sys.stderr)
        print("Inside the Action these come from repo Secrets; locally from .env.",
              file=sys.stderr)
        sys.exit(1)
    # RUNPOD_NETWORK_VOLUME_ID is optional: when set, pods attach it (and are
    # thus region-locked to its datacenter). Left empty, pods run in ANY
    # datacenter with free GPUs and use ephemeral container disk — R2 is the
    # source of truth either way, so nothing is lost by going volumeless.
    optional = ["RUNPOD_GPU_TYPE", "RUNPOD_TIMEOUT_MIN",
                "RUNPOD_CONTAINER_REGISTRY_AUTH_ID", "RUNPOD_NETWORK_VOLUME_ID"]
    return {k: os.environ[k] for k in required + optional if os.environ.get(k)}


# ── RunPod GraphQL ───────────────────────────────────────────────────────────
def runpod_gql(api_key: str, query: str, variables: dict | None = None,
               retries: int = 3, retry_backoff_s: float = 2.0) -> dict:
    """POST a GraphQL operation. Retries on transient network errors."""
    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            resp = requests.post(
                f"{RUNPOD_GRAPHQL}?api_key={api_key}",
                json={"query": query, "variables": variables or {}},
                timeout=30,
            )
            resp.raise_for_status()
            body = resp.json()
            if body.get("errors"):
                # GraphQL semantic errors are not retryable — surface immediately.
                raise RuntimeError(f"RunPod GraphQL error: {body['errors']}")
            return body["data"]
        except (requests.RequestException, ValueError) as e:
            last_exc = e
            if attempt < retries - 1:
                wait = retry_backoff_s * (2 ** attempt)
                print(f"  RunPod transient error ({e}); retry in {wait}s")
                time.sleep(wait)
    raise RuntimeError(f"RunPod GraphQL exhausted retries: {last_exc}")


def query_balance(api_key: str) -> float | None:
    """Best-effort balance check. Returns None if the field schema has shifted —
    we never want a balance-API regression to block a real batch."""
    try:
        data = runpod_gql(api_key,
                          "query { myself { currentSpendPerHr clientBalance } }")
        return float(data["myself"]["clientBalance"])
    except Exception as e:
        print(f"[batch] balance query failed: {e}; skipping cost guard")
        return None


def provision_pod(api_key: str, image: str, gpu_type: str,
                  network_volume_id: str | None,
                  pod_name: str, env_vars: dict[str, str],
                  registry_auth_id: str | None = None) -> str:
    """Create an on-demand pod. Returns pod_id.

    `network_volume_id`, when set, attaches a persistent volume at the data
    dir — but it region-locks the pod to the volume's datacenter. Left None,
    the pod runs in any datacenter with a free GPU and uses ephemeral
    container disk; the pipeline still syncs through R2.

    `registry_auth_id`, when set, references a RunPod Container Registry Auth
    credential so a *private* GHCR image can be pulled. Leave it None for a
    public image.
    """
    env_list = [{"key": k, "value": v} for k, v in env_vars.items()]
    pod_input = {
        "cloudType": "ALL", "gpuCount": 1, "gpuTypeId": gpu_type,
        "name": pod_name, "imageName": image,
        "dockerArgs": "bash /app/scripts/pod_campaign.sh",
        # Ephemeral working disk. The image lives outside this; data/ holds a
        # handful of SDFs/CSVs pulled from R2, so 30 GB is ample headroom.
        "containerDiskInGb": 30,
        "env": env_list,
    }
    if network_volume_id:
        pod_input["networkVolumeId"] = network_volume_id
        pod_input["volumeMountPath"] = "/app/autonomous_drug_discovery/data"
    if registry_auth_id:
        pod_input["containerRegistryAuthId"] = registry_auth_id
    variables = {"input": pod_input}
    data = runpod_gql(api_key,
                      """mutation Deploy($input: PodFindAndDeployOnDemandInput!) {
                           podFindAndDeployOnDemand(input: $input) { id }
                         }""",
                      variables)
    pod_id = data.get("podFindAndDeployOnDemand", {}).get("id")
    if not pod_id:
        raise RuntimeError(f"RunPod returned no pod id: {data}")
    return pod_id


def query_pod_uptime(api_key: str, pod_id: str) -> int | None:
    """Return the pod's current container uptime in seconds, or None if the
    pod is gone / not yet running. A pod that keeps reporting near-zero uptime
    while wall-clock advances is crash-looping (its command exits immediately
    and RunPod restarts the container)."""
    try:
        data = runpod_gql(
            api_key,
            "query Status($input: PodFilter!) { pod(input: $input) "
            "{ runtime { uptimeInSeconds } } }",
            {"input": {"podId": pod_id}}, retries=1,
        )
        rt = (data.get("pod") or {}).get("runtime") or {}
        return rt.get("uptimeInSeconds")
    except Exception:
        return None


def terminate_pod(api_key: str, pod_id: str) -> None:
    """Best-effort terminate. Logged on failure but never raises — orphaned
    pods are visible on the RunPod console and the cost guard catches them
    on the next batch."""
    try:
        runpod_gql(api_key,
                   "mutation Kill($input: PodTerminateInput!) { podTerminate(input: $input) }",
                   {"input": {"podId": pod_id}}, retries=1)
    except Exception as e:
        print(f"  WARNING: terminate failed for pod {pod_id}: {e}")


# ── R2 via rclone ────────────────────────────────────────────────────────────
def rclone(args: list[str], check: bool = True, capture: bool = True) -> subprocess.CompletedProcess:
    cmd = ["rclone", *args]
    return subprocess.run(cmd, capture_output=capture, text=True, check=check)


def r2_remote(env: dict[str, str]) -> str:
    return f"r2:{env['R2_BUCKET']}"


def pull_telemetry(env: dict[str, str], dest: Path) -> bool:
    """Try to pull telemetry.db from R2. Returns True on success — a missing
    file is *not* an error (fresh bucket → no skip-list)."""
    remote = f"{r2_remote(env)}/telemetry.db"
    res = rclone(["copyto", remote, str(dest)], check=False)
    if res.returncode == 0:
        return True
    if "directory not found" in (res.stderr or "").lower() or "not found" in (res.stderr or "").lower():
        return False
    print(f"  WARNING: telemetry pull non-zero ({res.returncode}): {res.stderr}")
    return False


def list_sentinels(env: dict[str, str]) -> dict[str, dict]:
    """Return {sentinel_key: {outcome, modtime}} from r2:bucket/sentinels/."""
    remote = f"{r2_remote(env)}/sentinels"
    res = rclone(["lsjson", remote], check=False)
    if res.returncode != 0:
        # No sentinel directory yet → empty.
        return {}
    try:
        entries = json.loads(res.stdout or "[]")
    except json.JSONDecodeError:
        return {}
    out: dict[str, dict] = {}
    for e in entries:
        name = e["Name"]
        if not (name.endswith(".done") or name.endswith(".failed")):
            continue
        key, outcome = name.rsplit(".", 1)
        out[key] = {"outcome": outcome, "modtime": e.get("ModTime", "")}
    return out


def read_sentinel(env: dict[str, str], key: str, outcome: str) -> str:
    """Fetch a sentinel's body so we can include error tails in the report."""
    remote = f"{r2_remote(env)}/sentinels/{key}.{outcome}"
    res = rclone(["cat", remote], check=False)
    return res.stdout if res.returncode == 0 else ""


def remove_sentinel(env: dict[str, str], key: str, outcome: str) -> None:
    remote = f"{r2_remote(env)}/sentinels/{key}.{outcome}"
    rclone(["delete", remote], check=False)


# ── Skip-list computation ────────────────────────────────────────────────────
def compute_skip_list(telemetry_path: Path, targets: Iterable[str], mode: str,
                      window_hours: int = SKIP_WINDOW_HOURS) -> set[str]:
    """Return the subset of targets that have a successful generation campaign
    in the given mode newer than `window_hours` ago.

    Public so the test suite can pin behaviour without spinning up the rest of
    the dispatcher.
    """
    if not telemetry_path.exists():
        return set()
    cutoff = (datetime.now(timezone.utc) - timedelta(hours=window_hours)).isoformat()
    conn = sqlite3.connect(str(telemetry_path))
    try:
        skip: set[str] = set()
        for target in targets:
            patterns = (f"%/{target}_manifest.json", f"%/{target}_pocket%_manifest.json")
            rows = conn.execute(
                """SELECT parameters FROM runs
                   WHERE module_name = '02_generation'
                     AND status = 'success'
                     AND started_at > ?
                     AND (input_path LIKE ? OR input_path LIKE ?)""",
                (cutoff, *patterns),
            ).fetchall()
            for (params_json,) in rows:
                try:
                    if json.loads(params_json or "{}").get("mode") == mode:
                        skip.add(target)
                        break
                except json.JSONDecodeError:
                    continue
        return skip
    finally:
        conn.close()


# ── Cost guard ───────────────────────────────────────────────────────────────
def estimate_cost_usd(num_targets: int, parallelism: int,
                      timeout_min: int = POD_TIMEOUT_MIN_DEFAULT) -> float:
    """Worst-case headroom: assume every pod runs to ceiling, every pod is one
    GPU-hour at the most-expensive rate. Padded by COST_SAFETY_MULTIPLIER."""
    waves = (num_targets + parallelism - 1) // parallelism
    gpu_hours = waves * parallelism * (timeout_min / 60.0)
    return gpu_hours * EST_USD_PER_GPU_HOUR * COST_SAFETY_MULTIPLIER


# ── Worker ───────────────────────────────────────────────────────────────────
def run_one_target(target: str, mode: str, num_samples: int | None,
                   env: dict[str, str], dry_run: bool,
                   timeout_min: int) -> dict:
    """Provision one pod and wait for its sentinel. Returns an outcome dict."""
    pod_id: str | None = None
    sentinel_key = f"{target}-{mode}-{uuid.uuid4().hex[:8]}"
    pod_name = f"agent-harness-{target}-{int(time.time())}"

    pod_env = {
        "R2_BUCKET": env["R2_BUCKET"],
        "TARGET": target,
        "MODE": mode,
        "NUM": str(num_samples) if num_samples else "",
        "SENTINEL_KEY": sentinel_key,
        "RUNPOD_TIMEOUT_MIN": env.get("RUNPOD_TIMEOUT_MIN", str(timeout_min)),
        "RCLONE_CONFIG_R2_TYPE": env["RCLONE_CONFIG_R2_TYPE"],
        "RCLONE_CONFIG_R2_PROVIDER": env["RCLONE_CONFIG_R2_PROVIDER"],
        "RCLONE_CONFIG_R2_ACCESS_KEY_ID": env["RCLONE_CONFIG_R2_ACCESS_KEY_ID"],
        "RCLONE_CONFIG_R2_SECRET_ACCESS_KEY": env["RCLONE_CONFIG_R2_SECRET_ACCESS_KEY"],
        "RCLONE_CONFIG_R2_ENDPOINT": env["RCLONE_CONFIG_R2_ENDPOINT"],
    }

    if dry_run:
        print(f"[dry-run] would provision pod for {target} (sentinel={sentinel_key})")
        return {"target": target, "sentinel_key": sentinel_key,
                "outcome": "dry-run", "pod_id": None}

    try:
        # Try each acceptable GPU type in order, falling through on supply
        # shortages. For a CPU-bound rdkit run any of these works; the list is
        # also fine for TargetDiff (all Ampere+). Override via RUNPOD_GPU_TYPE
        # (comma-separated to set your own priority order).
        gpu_types = [g.strip() for g in
                     env.get("RUNPOD_GPU_TYPE", DEFAULT_GPU_PRIORITY).split(",")
                     if g.strip()]
        pod_id = None
        last_supply_err: str | None = None
        for gpu in gpu_types:
            try:
                pod_id = provision_pod(
                    env["RUNPOD_API_KEY"], env["IMAGE"], gpu,
                    env.get("RUNPOD_NETWORK_VOLUME_ID"),
                    pod_name, pod_env,
                    registry_auth_id=env.get("RUNPOD_CONTAINER_REGISTRY_AUTH_ID"),
                )
                print(f"[pool] {target}: got {gpu}")
                break
            except Exception as e:
                if "SUPPLY_CONSTRAINT" in str(e) or "no longer any instances" in str(e):
                    print(f"[pool] {target}: {gpu} out of stock, trying next")
                    last_supply_err = str(e)
                    continue
                # A non-supply provision failure (RunPod 5xx, quota, bad input)
                # must not crash the whole batch — record it and let the pool
                # carry on with the other targets.
                print(f"[pool] {target}: provision failed — {e}")
                return {"target": target, "sentinel_key": sentinel_key,
                        "outcome": "provision_error", "pod_id": None,
                        "tail": str(e)}
        if pod_id is None:
            print(f"[pool] {target}: no GPU type had stock ({gpu_types})")
            return {"target": target, "sentinel_key": sentinel_key,
                    "outcome": "no_capacity", "pod_id": None,
                    "tail": last_supply_err or "all GPU types out of stock"}
        print(f"[pool] {target}: pod {pod_id} provisioned")

        # Poll for the sentinel (the source of truth for "did it finish?").
        # Alongside it, watch pod uptime so a crash-looping container — one
        # whose command dies on startup and gets restarted forever, never
        # producing a sentinel — is caught in minutes instead of hanging for
        # the full timeout. (This is exactly the failure a bad --mode caused.)
        start = time.monotonic()
        deadline = start + timeout_min * 60
        max_uptime = 0
        while time.monotonic() < deadline:
            time.sleep(SENTINEL_POLL_INTERVAL_S)
            sentinels = list_sentinels(env)
            if sentinel_key in sentinels:
                outcome = sentinels[sentinel_key]["outcome"]
                tail = read_sentinel(env, sentinel_key, outcome)
                print(f"[pool] {target}: sentinel={outcome}\n{tail}")
                return {"target": target, "sentinel_key": sentinel_key,
                        "outcome": outcome, "pod_id": pod_id, "tail": tail}

            uptime = query_pod_uptime(env["RUNPOD_API_KEY"], pod_id) or 0
            max_uptime = max(max_uptime, uptime)
            elapsed = time.monotonic() - start
            # Crash-loop heuristic: after a generous startup grace (image pull
            # can take minutes), if the container has never stayed alive longer
            # than CRASH_LOOP_UPTIME_S, its command is dying on startup.
            if elapsed > CRASH_LOOP_GRACE_S and max_uptime < CRASH_LOOP_UPTIME_S:
                print(f"[pool] {target}: crash-loop detected "
                      f"(elapsed {int(elapsed)}s, max container uptime {max_uptime}s, "
                      f"no sentinel) — aborting this attempt")
                return {"target": target, "sentinel_key": sentinel_key,
                        "outcome": "crash_loop", "pod_id": pod_id}

        print(f"[pool] {target}: outer timeout ({timeout_min} min) hit")
        return {"target": target, "sentinel_key": sentinel_key,
                "outcome": "timeout", "pod_id": pod_id}
    finally:
        if pod_id is not None and not dry_run:
            terminate_pod(env["RUNPOD_API_KEY"], pod_id)


# ── Telemetry merge after the batch drains ──────────────────────────────────
def merge_back_telemetry(env: dict[str, str], local_db: Path,
                         results: list[dict]) -> None:
    """For every successful campaign, pull its per-campaign DB from R2 and
    merge into the canonical local db. The script then re-uploads the merged
    DB at the very end (caller handles that)."""
    if not local_db.exists():
        # Bootstrap: a fresh DB created on first batch.
        from telemetry import TelemetryDB  # noqa: WPS433
        sys.path.insert(0, str(REPO_ROOT / "autonomous_drug_discovery"))
        TelemetryDB(str(local_db)).close()

    merge_script = REPO_ROOT / "scripts" / "merge_telemetry.py"
    for r in results:
        if r["outcome"] != "done":
            continue
        # Each campaign's telemetry is colocated with its outputs in R2 at
        # campaign_<id>/telemetry.db (the pod writes a private DB then sync).
        # We don't actually have that today — the pod writes to data/telemetry.db
        # directly, and rclone copies it back as the canonical one. For now,
        # successful campaigns share the canonical DB which was already
        # synchronised at script start, so this is a no-op. Left here as the
        # extension point when per-campaign DBs land.
        _ = merge_script  # placeholder
    print("[batch] telemetry merge: no per-campaign DBs to fold in (canonical DB is authoritative)")


# ── Entry point ──────────────────────────────────────────────────────────────
def parse_targets(raw: str) -> list[str]:
    """Accept whitespace- or comma-separated codes; de-dupe preserving order."""
    tokens = re.split(r"[\s,]+", raw.strip())
    seen: list[str] = []
    for t in tokens:
        if not t:
            continue
        if not TARGET_PATTERN.match(t):
            raise ValueError(f"invalid target code: {t!r}")
        if t not in seen:
            seen.append(t)
    return seen


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--targets", required=True, type=str,
                    help='Whitespace- or comma-separated PDB codes')
    ap.add_argument("--mode", default="targetdiff",
                    choices=["rdkit", "targetdiff", "pocket2mol", "simulation"])
    ap.add_argument("--num-samples", type=int, default=30, dest="num_samples")
    ap.add_argument("--parallelism", type=int, default=5)
    ap.add_argument("--force", action="store_true",
                    help="Disable the 24h idempotency skip-list")
    ap.add_argument("--timeout-min", type=int, default=POD_TIMEOUT_MIN_DEFAULT)
    ap.add_argument("--dry-run", action="store_true",
                    help="Print plan + GraphQL bodies; provision no pods")
    args = ap.parse_args()

    targets = parse_targets(args.targets)
    if not targets:
        print("ERROR: no valid targets", file=sys.stderr)
        return 1

    env = load_env()

    # 1. Skip-list from telemetry (pulled from R2).
    if not args.force:
        if pull_telemetry(env, TELEMETRY_PATH):
            skip = compute_skip_list(TELEMETRY_PATH, targets, args.mode)
            if skip:
                print(f"[batch] skipping (fresh successful campaign within "
                      f"{SKIP_WINDOW_HOURS}h): {' '.join(sorted(skip))}")
                targets = [t for t in targets if t not in skip]

    if not targets:
        print("[batch] nothing to do (everything in the skip-list)")
        return 0

    # 2. Fetch missing PDBs.
    print(f"[batch] ensuring PDBs for {len(targets)} target(s) ...")
    fetch_cmd = [sys.executable, str(REPO_ROOT / "scripts" / "fetch_pdb.py"),
                 "--targets", *targets]
    res = subprocess.run(fetch_cmd, check=False)
    if res.returncode != 0:
        print("[batch] fetch_pdb failed; aborting", file=sys.stderr)
        return 1

    # 3. Cost guard.
    est = estimate_cost_usd(len(targets), args.parallelism, args.timeout_min)
    print(f"[batch] estimated worst-case spend: ${est:.2f}")
    balance = query_balance(env["RUNPOD_API_KEY"])
    if balance is not None:
        if balance < est:
            print(f"ERROR: RunPod balance ${balance:.2f} below worst-case "
                  f"estimate ${est:.2f}. Top up or shrink the target list.",
                  file=sys.stderr)
            return 1
        print(f"[batch] balance OK: ${balance:.2f} ≥ ${est:.2f}")

    if args.dry_run:
        print("[batch] dry-run — would dispatch:")
        for t in targets:
            print(f"  {t}")
        return 0

    # 4. Seed R2 with PDBs.
    print("[batch] seeding R2 with input PDBs ...")
    rclone(["copy", str(DATA_DIR / "processed"), f"{r2_remote(env)}/processed",
            "--include", "*.pdb"], check=False)

    # 5. Worker pool. Each target has its own retry budget; on .failed the
    #    worker re-submits up to MAX_RETRIES_PER_TARGET times.
    results: dict[str, dict] = {}
    attempts: dict[str, int] = {t: 0 for t in targets}
    work_queue = list(targets)

    def run_with_retries(target: str) -> dict:
        last: dict | None = None
        while attempts[target] < MAX_RETRIES_PER_TARGET:
            attempts[target] += 1
            last = run_one_target(target, args.mode, args.num_samples,
                                  env, args.dry_run, args.timeout_min)
            if last["outcome"] == "done":
                return last
            # A crash-loop is deterministic (bad image/config/args) — retrying
            # just burns money on the same failure. no_capacity means every GPU
            # type was out of stock; an instant retry won't conjure supply.
            # Both are non-retryable here.
            if last["outcome"] in ("crash_loop", "no_capacity"):
                print(f"[pool] {target}: {last['outcome']} is non-retryable — giving up")
                return last
            print(f"[pool] {target}: outcome={last['outcome']}, "
                  f"attempt {attempts[target]}/{MAX_RETRIES_PER_TARGET}")
        return last or {"target": target, "outcome": "exhausted"}

    print(f"[batch] dispatching {len(work_queue)} target(s) with "
          f"parallelism={args.parallelism} ...")
    with ThreadPoolExecutor(max_workers=args.parallelism) as ex:
        futures = {ex.submit(run_with_retries, t): t for t in work_queue}
        for fut in as_completed(futures):
            r = fut.result()
            results[r["target"]] = r

    # 6. Pull final state from R2.
    print("[batch] syncing campaign outputs back from R2 ...")
    rclone(["copy", r2_remote(env), str(DATA_DIR), "--progress"], check=False)

    # 7. (placeholder) merge per-campaign telemetry. Today this is a no-op
    #    because the pod writes to the canonical DB and we just rclone'd it back.
    merge_back_telemetry(env, TELEMETRY_PATH, list(results.values()))

    # 8. Re-upload canonical telemetry (last-writer-wins; safe because we
    #    pulled at the start and nothing else writes during the batch).
    rclone(["copyto", str(TELEMETRY_PATH), f"{r2_remote(env)}/telemetry.db"],
           check=False)

    # 9. Report.
    successes = [t for t, r in results.items() if r["outcome"] == "done"]
    failures = [(t, r["outcome"]) for t, r in results.items() if r["outcome"] != "done"]
    print()
    print("=" * 70)
    print(f"BATCH SUMMARY: {len(successes)}/{len(results)} successful")
    for t in successes:
        print(f"  ✓ {t}")
    for t, outcome in failures:
        print(f"  ✗ {t} ({outcome})")
    print("=" * 70)

    # 10. Write a machine-readable summary the Action can post to the run page.
    out_dir = REPO_ROOT / "data"
    out_dir.mkdir(exist_ok=True)
    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": args.mode, "num_samples": args.num_samples,
        "parallelism": args.parallelism,
        "results": [
            {"target": t, **{k: v for k, v in r.items() if k != "target"}}
            for t, r in results.items()
        ],
    }
    (out_dir / "batch_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[batch] wrote {out_dir / 'batch_summary.json'}")

    return 0 if not failures else 2


if __name__ == "__main__":
    sys.exit(main())
