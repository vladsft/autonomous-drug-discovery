#!/usr/bin/env bash
# Apply in-repo patches to the TargetDiff submodule.
#
# Why this exists: TargetDiff (guanjq/targetdiff @ 142f1eb) was written
# against NumPy <1.20 and uses removed aliases (np.long, np.bool). The
# submodule is pinned at the upstream SHA — we don't fork — so each
# fresh checkout has to re-apply the compatibility patches stored in
# modules/02_generation/targetdiff_patches/.
#
# Idempotent: skips patches that are already applied.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TARGETDIFF_DIR="$REPO_ROOT/autonomous_drug_discovery/modules/02_generation/targetdiff"
PATCH_DIR="$REPO_ROOT/autonomous_drug_discovery/modules/02_generation/targetdiff_patches"

if [ ! -d "$TARGETDIFF_DIR/.git" ] && [ ! -f "$TARGETDIFF_DIR/.git" ]; then
    echo "TargetDiff submodule not initialized. Run:" >&2
    echo "  git submodule update --init --recursive" >&2
    exit 1
fi

shopt -s nullglob
patches=("$PATCH_DIR"/*.patch)
if [ ${#patches[@]} -eq 0 ]; then
    echo "No patches found in $PATCH_DIR — nothing to do."
    exit 0
fi

for patch in "${patches[@]}"; do
    name=$(basename "$patch")
    if git -C "$TARGETDIFF_DIR" apply --check --reverse "$patch" >/dev/null 2>&1; then
        echo "skip   $name (already applied)"
    elif git -C "$TARGETDIFF_DIR" apply --check "$patch" >/dev/null 2>&1; then
        git -C "$TARGETDIFF_DIR" apply "$patch"
        echo "apply  $name"
    else
        echo "ERROR  $name does not apply cleanly. Inspect $patch and the submodule state." >&2
        exit 1
    fi
done
