#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export BOLTZ_APPTAINER_IMAGE="${BOLTZ_APPTAINER_IMAGE:-$ROOT_DIR/containers/current}"
export BOLTZ_CONTAINER_SITEPKGS="${BOLTZ_CONTAINER_SITEPKGS:-$ROOT_DIR/sitepkgs_bundle}"
export BOLTZ_HOST_PYTHON="${BOLTZ_HOST_PYTHON:-$ROOT_DIR/containers/current/usr/local/apps/pyenv/versions/miniforge3-24.11.3-2/envs/boltz-conda/bin/python}"

if [ ! -x "$BOLTZ_HOST_PYTHON" ]; then
  echo "Host Python not found: $BOLTZ_HOST_PYTHON" >&2
  echo "Expected bundled runtime inside $ROOT_DIR/containers/current" >&2
  exit 1
fi

exec "$BOLTZ_HOST_PYTHON" "$ROOT_DIR/boltz-2_wrapper.py" "$@"
