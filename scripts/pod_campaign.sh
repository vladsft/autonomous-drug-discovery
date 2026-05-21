#!/usr/bin/env bash
# Runs inside the GPU pod (baked into the image; invoked by cloud_run.sh or
# scripts/batch_cloud_run.py).
#
# The pod has no state of its own: it pulls everything it needs from R2, runs
# one campaign, pushes the results back, and writes a sentinel file at
# r2:bucket/sentinels/<SENTINEL_KEY>.{done,failed}. The dispatcher polls those
# sentinels to know which pods finished — pod-status APIs are eventually
# consistent and lie on ungraceful exit; a file the pod itself wrote just
# before exiting is authoritative.
#
# Expected environment:
#   R2_BUCKET, TARGET, MODE              required
#   NUM                                  optional, --num_samples override
#   SENTINEL_KEY                         optional; defaults to "${TARGET}-${MODE}-$(date +%s)"
#   RUNPOD_TIMEOUT_MIN                   inner-fuse minutes (default 360)
#   RCLONE_CONFIG_R2_*                   rclone config via env

set -euo pipefail

cd /app/autonomous_drug_discovery

: "${R2_BUCKET:?}" "${TARGET:?}" "${MODE:?}"
remote="r2:${R2_BUCKET}"
SENTINEL_KEY="${SENTINEL_KEY:-${TARGET}-${MODE}-$(date +%s)}"

# Write a sentinel even on abnormal exit so the dispatcher never has to guess.
# Captures the exit code + the last 50 log lines. Idempotent: a re-run of this
# script for the same SENTINEL_KEY overwrites the previous outcome.
write_sentinel() {
    local outcome="$1" rc="$2"
    local tmp; tmp="$(mktemp)"
    {
        echo "sentinel_key: ${SENTINEL_KEY}"
        echo "target: ${TARGET}"
        echo "mode: ${MODE}"
        echo "exit_code: ${rc}"
        echo "finished_at: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "--- last 50 log lines ---"
        [ -f data/run.log ] && tail -n 50 data/run.log || echo "(no run.log)"
    } > "${tmp}"
    rclone copyto "${tmp}" "${remote}/sentinels/${SENTINEL_KEY}.${outcome}" || \
        echo "[pod] WARNING: sentinel upload failed; dispatcher may time us out"
    rm -f "${tmp}"
}

# trap fires on any non-zero exit (set -e), or on signals. We always want a
# sentinel — even if rclone or the orchestrator dies partway through.
on_exit() {
    local rc=$?
    if [ "${rc}" -ne 0 ]; then
        write_sentinel "failed" "${rc}"
    fi
}
trap on_exit EXIT

# Inner fuse: enforce a hard wallclock on the pipeline itself so a hung
# stage (P2Rank deadlock, vina lockup, mamba network stall) can never keep
# the pod billing past the outer cloud_run.sh timeout. Leave 5 min headroom
# so R2 push still lands before cloud_run.sh tears the pod down.
TIMEOUT_MIN="${RUNPOD_TIMEOUT_MIN:-360}"
PIPELINE_TIMEOUT_MIN=$(( TIMEOUT_MIN - 5 ))
[ "${PIPELINE_TIMEOUT_MIN}" -lt 5 ] && PIPELINE_TIMEOUT_MIN=5

echo "[pod] === Agent Harness campaign: ${TARGET} / ${MODE} / num=${NUM:-default} ==="
echo "[pod] === pipeline fuse: ${PIPELINE_TIMEOUT_MIN} min ==="

echo "[pod] Pulling input state (PDBs, telemetry) from ${remote} ..."
rclone copy "${remote}" data --progress

num_flag=()
[ -n "${NUM:-}" ] && num_flag=(--num_samples "${NUM}")

echo "[pod] Running the pipeline ..."
set +e
# --foreground so signals propagate to the orchestrator and its children;
# without it `timeout` sends SIGTERM only to itself and the pipeline lingers.
# --kill-after=60s upgrades to SIGKILL if SIGTERM is ignored.
timeout --foreground --kill-after=60s "${PIPELINE_TIMEOUT_MIN}m" \
    python orchestrator.py run "data/processed/${TARGET}.pdb" \
    --mode "${MODE}" "${num_flag[@]}" 2>&1 | tee "data/run.log"
rc=${PIPESTATUS[0]}
set -e

if [ "${rc}" -eq 124 ]; then
    echo "[pod] FATAL: pipeline exceeded ${PIPELINE_TIMEOUT_MIN}-minute fuse." | tee -a "data/run.log"
fi

echo "[pod] Pipeline exit code: ${rc}. Pushing results to ${remote} ..."
rclone copy data "${remote}" --progress

# Outcome is decided by the pipeline rc; sentinel is written explicitly here so
# the .done case doesn't depend on the EXIT trap (which fires only on rc != 0).
if [ "${rc}" -eq 0 ]; then
    write_sentinel "done" "${rc}"
fi

echo "[pod] === done (rc=${rc}, sentinel=${SENTINEL_KEY}) ==="
exit "${rc}"
