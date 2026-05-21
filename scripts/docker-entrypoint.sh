#!/usr/bin/env bash
# Container entrypoint for the Agent Harness image.
#
# Runs whatever it is handed as a Python program inside the `base` conda
# environment, from the autonomous_drug_discovery/ directory. So:
#
#   docker run --rm agent-harness orchestrator.py run data/processed/1M17.pdb
#   docker run --rm agent-harness benchmark.py
#
# both Just Work. The orchestrator itself shells out to `conda run -n
# targetdiff_env ...` for the diffusion stage, so only `base` is needed here.
#
# Escape hatches:
#   docker run --entrypoint bash -it agent-harness      # interactive shell
#   docker run agent-harness bash -lc '...'             # ad-hoc command
set -euo pipefail

# Allow `bash`/`sh` as a first arg to drop into a shell without --entrypoint.
# We shift the shell name off so it isn't doubled as an arg to itself.
case "${1:-}" in
    bash|sh)
        shell="$1"; shift
        exec conda run --no-capture-output -n base "${shell}" "$@"
        ;;
esac

exec conda run --no-capture-output -n base python "$@"
