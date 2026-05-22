#!/usr/bin/env bash
# cloud_run.sh — run one campaign on a rented RunPod GPU, hands-free.
#
#   make cloud-run TARGET=2HYY MODE=targetdiff NUM=30
#
# What it does:
#   1. Seeds Cloudflare R2 with the input PDB + telemetry.db.
#   2. Provisions an on-demand GPU pod from the published GHCR image, with a
#      persistent network volume mounted at the data directory.
#   3. The pod pulls state from R2, runs the full pipeline, pushes results back
#      (see scripts/pod_campaign.sh).
#   4. Polls until the pod's command exits, then TERMINATES the pod — a trap
#      guarantees teardown even if this script is interrupted, so a GPU is
#      never left billing.
#   5. Pulls the results back into the local data/ directory.
#
# Credentials and configuration come from .env (see .env.example).
#
# NOTE: this script talks to the live RunPod GraphQL API. It cannot be tested
# without a funded RunPod account; expect to iterate on the first real run.
set -euo pipefail

# --- Locate the repo and load configuration ----------------------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${REPO_ROOT}/autonomous_drug_discovery/data"
ENV_FILE="${REPO_ROOT}/.env"

[ -f "${ENV_FILE}" ] || { echo "ERROR: ${ENV_FILE} not found — copy .env.example to .env."; exit 1; }
set -a; # shellcheck disable=SC1090
source "${ENV_FILE}"; set +a

# --- Run parameters (passed through by the Makefile) --------------------------
TARGET="${TARGET:-1M17}"
MODE="${MODE:-targetdiff}"
NUM="${NUM:-}"
IMAGE="${IMAGE:?IMAGE not set}"
GPU_TYPE="${RUNPOD_GPU_TYPE:-NVIDIA GeForce RTX 3090}"
# Polling ceiling. Configurable via .env so a long batch run can bump it
# without editing this script. Defaults to 6h; the polling loop ticks every
# 20s, so total iterations = TIMEOUT_MIN × 3.
TIMEOUT_MIN="${RUNPOD_TIMEOUT_MIN:-360}"

# --- Preflight ----------------------------------------------------------------
for tool in curl jq rclone; do
    command -v "${tool}" >/dev/null || { echo "ERROR: '${tool}' is required."; exit 1; }
done
: "${RUNPOD_API_KEY:?set RUNPOD_API_KEY in .env}"
: "${RUNPOD_NETWORK_VOLUME_ID:?set RUNPOD_NETWORK_VOLUME_ID in .env}"
: "${R2_BUCKET:?set R2_BUCKET in .env}"
: "${RCLONE_CONFIG_R2_ACCESS_KEY_ID:?set the R2 credentials in .env}"

PDB="${DATA_DIR}/processed/${TARGET}.pdb"
[ -f "${PDB}" ] || { echo "ERROR: input PDB not found: ${PDB}"; exit 1; }

GQL="https://api.runpod.io/graphql?api_key=${RUNPOD_API_KEY}"

# gql <query> [variables-json] — run a GraphQL operation, print the response.
# Using GraphQL *variables* (not string interpolation) keeps every value
# correctly typed and escaped.
gql() {
    local query="$1" vars="${2:-{}}"
    curl -sf -X POST "${GQL}" -H 'Content-Type: application/json' \
        -d "$(jq -n --arg q "${query}" --argjson v "${vars}" '{query: $q, variables: $v}')"
}

# --- 1. Seed R2 with the input state -----------------------------------------
echo "==> Seeding r2:${R2_BUCKET} with ${TARGET}.pdb + telemetry.db ..."
rclone copy "${DATA_DIR}/processed" "r2:${R2_BUCKET}/processed" --include "${TARGET}.pdb*"
[ -f "${DATA_DIR}/telemetry.db" ] && rclone copy "${DATA_DIR}/telemetry.db" "r2:${R2_BUCKET}/"

# --- 2. Provision the pod -----------------------------------------------------
# The pod env carries the R2 credentials so the pod's rclone needs no config
# file. dockerArgs runs pod_campaign.sh through the image's python entrypoint.
POD_NAME="agent-harness-${TARGET}-$(date +%s)"
ENV_JSON="$(jq -n \
    --arg bucket "${R2_BUCKET}" --arg target "${TARGET}" \
    --arg mode "${MODE}" --arg num "${NUM}" \
    --arg timeout "${TIMEOUT_MIN}" \
    --arg t "${RCLONE_CONFIG_R2_TYPE:-s3}" \
    --arg prov "${RCLONE_CONFIG_R2_PROVIDER:-Cloudflare}" \
    --arg ak "${RCLONE_CONFIG_R2_ACCESS_KEY_ID}" \
    --arg sk "${RCLONE_CONFIG_R2_SECRET_ACCESS_KEY:-}" \
    --arg ep "${RCLONE_CONFIG_R2_ENDPOINT:-}" \
    '[ {key:"R2_BUCKET",value:$bucket}, {key:"TARGET",value:$target},
       {key:"MODE",value:$mode}, {key:"NUM",value:$num},
       {key:"RUNPOD_TIMEOUT_MIN",value:$timeout},
       {key:"RCLONE_CONFIG_R2_TYPE",value:$t},
       {key:"RCLONE_CONFIG_R2_PROVIDER",value:$prov},
       {key:"RCLONE_CONFIG_R2_ACCESS_KEY_ID",value:$ak},
       {key:"RCLONE_CONFIG_R2_SECRET_ACCESS_KEY",value:$sk},
       {key:"RCLONE_CONFIG_R2_ENDPOINT",value:$ep} ]')"

# When the GHCR image is private, RunPod needs a registry credential to pull
# it. Register a GitHub PAT (read:packages) under RunPod → Settings → Container
# Registry Auth, then put the resulting auth ID in RUNPOD_CONTAINER_REGISTRY_AUTH_ID.
# Left empty ⇒ the image is assumed public and no auth is attached.
AUTH_ID="${RUNPOD_CONTAINER_REGISTRY_AUTH_ID:-}"
INPUT_JSON="$(jq -n \
    --arg name "${POD_NAME}" --arg image "${IMAGE}" --arg gpu "${GPU_TYPE}" \
    --arg vol "${RUNPOD_NETWORK_VOLUME_ID}" --argjson env "${ENV_JSON}" \
    --arg authid "${AUTH_ID}" \
    '{ cloudType: "ALL", gpuCount: 1, gpuTypeId: $gpu,
       name: $name, imageName: $image,
       dockerArgs: "bash /app/scripts/pod_campaign.sh",
       containerDiskInGb: 20,
       networkVolumeId: $vol,
       volumeMountPath: "/app/autonomous_drug_discovery/data",
       env: $env }
     + (if $authid == "" then {} else {containerRegistryAuthId: $authid} end)')"

DEPLOY='mutation Deploy($input: PodFindAndDeployOnDemandInput!) {
  podFindAndDeployOnDemand(input: $input) { id machineId }
}'

echo "==> Provisioning ${GPU_TYPE} pod (${POD_NAME}) ..."
RESP="$(gql "${DEPLOY}" "$(jq -n --argjson i "${INPUT_JSON}" '{input: $i}')")"
POD_ID="$(echo "${RESP}" | jq -r '.data.podFindAndDeployOnDemand.id // empty')"
if [ -z "${POD_ID}" ]; then
    echo "ERROR: pod creation failed. RunPod said:"; echo "${RESP}" | jq .; exit 1
fi
echo "    pod id: ${POD_ID}"

# Guarantee teardown — even on Ctrl-C or an error below.
terminate() {
    echo "==> Terminating pod ${POD_ID} ..."
    gql 'mutation Kill($input: PodTerminateInput!) { podTerminate(input: $input) }' \
        "$(jq -n --arg p "${POD_ID}" '{input: {podId: $p}}')" >/dev/null \
        && echo "    terminated." || echo "    WARNING: terminate call failed — check runpod.io!"
}
trap terminate EXIT

# --- 3. Wait for the pod's command to finish ---------------------------------
# The pod's runtime is null before the container starts and again after its
# command exits. Wait for it to start, then wait for it to stop. The ceiling
# is configured via RUNPOD_TIMEOUT_MIN in .env (default 6h).
echo "==> Waiting for the campaign to finish (polling every 20s, ${TIMEOUT_MIN}m ceiling) ..."
STATUS='query Status($input: PodFilter!) { pod(input: $input) { runtime { uptimeInSeconds } } }'
STATUS_VARS="$(jq -n --arg p "${POD_ID}" '{input: {podId: $p}}')"
started=0
ticks=$((TIMEOUT_MIN * 3))   # 3 ticks/min × TIMEOUT_MIN
timed_out=1
for _ in $(seq 1 "${ticks}"); do
    sleep 20
    UPTIME="$(gql "${STATUS}" "${STATUS_VARS}" | jq -r '.data.pod.runtime.uptimeInSeconds // empty')"
    if [ -n "${UPTIME}" ]; then
        started=1
        printf '\r    running — uptime %ss      ' "${UPTIME}"
    elif [ "${started}" -eq 1 ]; then
        echo; echo "==> Pod command exited."
        timed_out=0
        break
    fi
done
if [ "${started}" -eq 1 ] && [ "${timed_out}" -eq 1 ]; then
    echo; echo "WARNING: pod still running after ${TIMEOUT_MIN} minutes — tearing it down."
    echo "         Inspect runpod.io for partial output; bump RUNPOD_TIMEOUT_MIN if this was premature."
fi

# --- 4. Pull results back -----------------------------------------------------
echo "==> Syncing results from r2:${R2_BUCKET} into ${DATA_DIR} ..."
rclone copy "r2:${R2_BUCKET}" "${DATA_DIR}" --progress

echo
echo "==> --- run.log tail -------------------------------------------------"
[ -f "${DATA_DIR}/run.log" ] && tail -n 25 "${DATA_DIR}/run.log" || echo "(no run.log retrieved)"
echo "==> --------------------------------------------------------------------"
echo "Done. Campaign output is in ${DATA_DIR}/. (Pod teardown follows.)"
