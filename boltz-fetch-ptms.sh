#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_IMAGE="$ROOT_DIR/containers/current"
if [ -e "$DEFAULT_IMAGE" ]; then
  DEFAULT_IMAGE="$(readlink -f "$DEFAULT_IMAGE")"
fi
export BOLTZ_APPTAINER_IMAGE="${BOLTZ_APPTAINER_IMAGE:-$DEFAULT_IMAGE}"
export BOLTZ_CONTAINER_SITEPKGS="${BOLTZ_CONTAINER_SITEPKGS:-$ROOT_DIR/sitepkgs_bundle}"

if [ ! -e "$BOLTZ_APPTAINER_IMAGE" ]; then
  echo "Container image not found: $BOLTZ_APPTAINER_IMAGE" >&2
  exit 1
fi

exec apptainer exec --cleanenv "$BOLTZ_APPTAINER_IMAGE" \
  python "$ROOT_DIR/boltz_fetch_ptms.py" "$@"
