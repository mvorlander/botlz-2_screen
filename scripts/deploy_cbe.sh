#!/bin/bash
set -euo pipefail

TARGET_ROOT="${1:-/resources/AF2_PPI_tools/boltz}"
BASE_IMAGE="${BASE_IMAGE:-$TARGET_ROOT/containers/boltz_screen_base}"
BUILD_ROOT="${BOLTZ_BUILD_ROOT:-$TARGET_ROOT/tmp/$USER}"
APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$BUILD_ROOT/apptainer-tmp}"
APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$BUILD_ROOT/apptainer-cache}"
BUILD_TAG="${BOLTZ_BUILD_TAG:-$(date +%Y%m%d_%H%M%S)}"
NEW_CONTAINER="$TARGET_ROOT/containers/boltz_screen_${BUILD_TAG}.sif"

mkdir -p "$TARGET_ROOT"
rsync -av --delete \
  --exclude '.git' \
  --exclude '__pycache__' \
  ./ "$TARGET_ROOT/"

mkdir -p "$TARGET_ROOT/containers"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
chmod 700 "$BUILD_ROOT" "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR" 2>/dev/null || true

if [ ! -e "$BASE_IMAGE" ]; then
  if [ -d "$TARGET_ROOT/containers/current" ] || [ -f "$TARGET_ROOT/containers/current" ]; then
    cp -a "$TARGET_ROOT/containers/current" "$BASE_IMAGE"
  else
    echo "Base image not found: $BASE_IMAGE" >&2
    echo "Set BASE_IMAGE explicitly for the first deployment." >&2
    exit 1
  fi
fi

sed "s|^From: .*|From: $BASE_IMAGE|" \
  "$TARGET_ROOT/containers/boltz_screen.def" > "$TARGET_ROOT/containers/boltz_screen.cbe.def"

rm -f "$NEW_CONTAINER"
TMPDIR="$APPTAINER_TMPDIR" \
APPTAINER_TMPDIR="$APPTAINER_TMPDIR" \
APPTAINER_CACHEDIR="$APPTAINER_CACHEDIR" \
apptainer build --fakeroot \
  "$NEW_CONTAINER" \
  "$TARGET_ROOT/containers/boltz_screen.cbe.def"

ln -sfn "$(basename "$NEW_CONTAINER")" "$TARGET_ROOT/containers/current"

chmod +x \
  "$TARGET_ROOT/boltz-screen.sh" \
  "$TARGET_ROOT/boltz-analysis.sh" \
  "$TARGET_ROOT/boltz-fetch-ptms.sh" \
  "$TARGET_ROOT/scripts/deploy_cbe.sh"

echo "Installed to $TARGET_ROOT"
echo "Container target: $TARGET_ROOT/containers/current"
echo "Apptainer temp dir: $APPTAINER_TMPDIR"
echo "Apptainer cache dir: $APPTAINER_CACHEDIR"
