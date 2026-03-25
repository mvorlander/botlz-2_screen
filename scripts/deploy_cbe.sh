#!/bin/bash
set -euo pipefail

TARGET_ROOT="${1:-/resources/AF2_PPI_tools/boltz}"
BASE_IMAGE="${BASE_IMAGE:-/groups/plaschka/shared/software/boltz-2/boltz2:develop}"

mkdir -p "$TARGET_ROOT"
rsync -av --delete \
  --exclude '.git' \
  --exclude '__pycache__' \
  ./ "$TARGET_ROOT/"

mkdir -p "$TARGET_ROOT/containers"

sed "s|^From: .*|From: $BASE_IMAGE|" \
  "$TARGET_ROOT/containers/boltz_screen.def" > "$TARGET_ROOT/containers/boltz_screen.cbe.def"

apptainer build --fakeroot --sandbox \
  "$TARGET_ROOT/containers/boltz_screen" \
  "$TARGET_ROOT/containers/boltz_screen.cbe.def"

ln -sfn boltz_screen "$TARGET_ROOT/containers/current"

chmod +x \
  "$TARGET_ROOT/boltz-screen.sh" \
  "$TARGET_ROOT/boltz-analysis.sh" \
  "$TARGET_ROOT/boltz-fetch-ptms.sh" \
  "$TARGET_ROOT/scripts/deploy_cbe.sh"

echo "Installed to $TARGET_ROOT"
echo "Container target: $TARGET_ROOT/containers/current"
