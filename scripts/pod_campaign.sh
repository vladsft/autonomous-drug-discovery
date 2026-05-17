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

echo "[pod] === Agent Harness campaign: ${TARGET} / ${MODE} / num=${NUM:-default} ==="

echo "[pod] Pulling input state (PDBs, telemetry) from ${remote} ..."
rclone copy "${remote}" data --progress

num_flag=()
[ -n "${NUM:-}" ] && num_flag=(--num_samples "${NUM}")

echo "[pod] Running the pipeline ..."
set +e
python orchestrator.py run "data/processed/${TARGET}.pdb" \
    --mode "${MODE}" "${num_flag[@]}" 2>&1 | tee "data/run.log"
rc=${PIPESTATUS[0]}
set -e

echo "[pod] Pipeline exit code: ${rc}. Pushing results to ${remote} ..."
rclone copy data "${remote}" --progress

echo "[pod] === done (rc=${rc}) ==="
exit "${rc}"
