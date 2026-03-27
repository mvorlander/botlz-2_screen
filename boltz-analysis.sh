#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(pwd -P)"
DEFAULT_IMAGE="$ROOT_DIR/containers/current"
if [ -e "$DEFAULT_IMAGE" ]; then
  DEFAULT_IMAGE="$(readlink -f "$DEFAULT_IMAGE")"
fi
ANALYSIS_IMAGE="${BOLTZ_APPTAINER_IMAGE:-$ROOT_DIR/containers/boltz_screen}"
export BOLTZ_CONTAINER_SITEPKGS="${BOLTZ_CONTAINER_SITEPKGS:-$ROOT_DIR/sitepkgs_bundle}"
export APPTAINERENV_PYTHONNOUSERSITE=1

if [ -d "$BOLTZ_CONTAINER_SITEPKGS" ]; then
  export APPTAINERENV_PYTHONPATH="$BOLTZ_CONTAINER_SITEPKGS${PYTHONPATH:+:$PYTHONPATH}"
fi

exec apptainer exec --cleanenv --no-mount hostfs \
  --bind "$ROOT_DIR:$ROOT_DIR" \
  --bind "$WORK_DIR:$WORK_DIR" \
  "$ANALYSIS_IMAGE" \
  /usr/local/apps/pyenv/versions/miniforge3-24.11.3-2/envs/boltz-conda/bin/python \
  "$ROOT_DIR/boltz_analysis.py" "$@"
