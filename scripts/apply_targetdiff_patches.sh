#!/usr/bin/env bash
# Apply in-repo patches to the TargetDiff submodule.
#
# Why this exists: TargetDiff (guanjq/targetdiff @ 142f1eb) was written
# against NumPy <1.20 and uses removed aliases (np.long, np.bool). The
# submodule is pinned at the upstream SHA — we don't fork — so each fresh
# checkout has to re-apply the compatibility patches stored in
# modules/02_generation/targetdiff_patches/.
#
# Uses `patch`, not `git apply`, on purpose: inside the Docker image the
# submodule's `.git` file dangles (the parent .git is not in the build
# context), which makes every `git` invocation in that directory fail with
# "not a git repository". `patch` has no repository dependency, so this script
# behaves identically on a dev checkout and in the container.
#
# Idempotent: a patch that is already applied is detected (reverse dry-run)
# and skipped.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TARGETDIFF_DIR="$REPO_ROOT/autonomous_drug_discovery/modules/02_generation/targetdiff"
PATCH_DIR="$REPO_ROOT/autonomous_drug_discovery/modules/02_generation/targetdiff_patches"

# The submodule must actually be checked out. Test for a source file every
# patch expects — not for `.git`, which is unreliable inside the image.
if [ ! -f "$TARGETDIFF_DIR/datasets/protein_ligand.py" ]; then
    echo "TargetDiff submodule not checked out at $TARGETDIFF_DIR. Run:" >&2
    echo "  git submodule update --init --recursive" >&2
    exit 1
fi

shopt -s nullglob
patches=("$PATCH_DIR"/*.patch)
if [ ${#patches[@]} -eq 0 ]; then
    echo "No patches found in $PATCH_DIR — nothing to do."
    exit 0
fi

for patch_file in "${patches[@]}"; do
    name=$(basename "$patch_file")
    # -p1 strips the a/ b/ prefix; -d runs relative to the submodule root.
    if patch -p1 -R --dry-run --force -d "$TARGETDIFF_DIR" < "$patch_file" >/dev/null 2>&1; then
        echo "skip   $name (already applied)"
    elif patch -p1 --dry-run --force -d "$TARGETDIFF_DIR" < "$patch_file" >/dev/null 2>&1; then
        patch -p1 -d "$TARGETDIFF_DIR" < "$patch_file"
        echo "apply  $name"
    else
        echo "ERROR  $name does not apply cleanly. Inspect $patch_file and the submodule state." >&2
        exit 1
    fi
done
