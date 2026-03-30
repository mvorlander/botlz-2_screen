#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(pwd -P)"
DEFAULT_IMAGE="$ROOT_DIR/containers/current"
if [ -e "$DEFAULT_IMAGE" ]; then
  DEFAULT_IMAGE="$(readlink -f "$DEFAULT_IMAGE")"
fi
PREP_IMAGE="${BOLTZ_PREPARE_IMAGE:-$DEFAULT_IMAGE}"
export BOLTZ_APPTAINER_IMAGE="${BOLTZ_APPTAINER_IMAGE:-$DEFAULT_IMAGE}"
export BOLTZ_CONTAINER_SITEPKGS="${BOLTZ_CONTAINER_SITEPKGS:-$ROOT_DIR/sitepkgs_bundle}"
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE
unset APPTAINERENV_PYTHONPATH APPTAINERENV_PYTHONHOME APPTAINERENV_PYTHONUSERBASE
export APPTAINERENV_PYTHONNOUSERSITE=1
if [ -d "$BOLTZ_CONTAINER_SITEPKGS" ]; then
  export APPTAINERENV_PYTHONPATH="$BOLTZ_CONTAINER_SITEPKGS"
fi

if [ ! -e "$PREP_IMAGE" ]; then
  echo "Preparation container not found: $PREP_IMAGE" >&2
  exit 1
fi

exec apptainer exec --cleanenv --no-mount hostfs \
  --bind "$ROOT_DIR:$ROOT_DIR" \
  --bind "$WORK_DIR:$WORK_DIR" \
  "$PREP_IMAGE" \
  /bin/bash --noprofile --norc -lc 'export PATH="/usr/local/apps/pyenv/versions/miniforge3-24.11.3-2/envs/boltz-conda/bin:$PATH"; exec python -I "$@"' \
  _ "$ROOT_DIR/boltz_fetch_ptms.py" "$@"
