#!/usr/bin/env bash
# Runs inside the GPU pod (baked into the image; invoked by cloud_run.sh).
#
# The pod has no state of its own: it pulls everything it needs from R2, runs
# one campaign, and pushes the results back. R2 is the source of truth — the
# pod is disposable. rclone is configured entirely through RCLONE_CONFIG_R2_*
# environment variables that cloud_run.sh passes in.
#
# Expected environment: R2_BUCKET, TARGET, MODE, NUM (optional), and the
# RCLONE_CONFIG_R2_* set.
set -euo pipefail

cd /app/autonomous_drug_discovery

: "${R2_BUCKET:?}" "${TARGET:?}" "${MODE:?}"
remote="r2:${R2_BUCKET}"

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

echo "[pod] === done (rc=${rc}) ==="
exit "${rc}"
