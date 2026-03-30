#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(pwd -P)"
runtime_base_default() {
  if [ -n "${BOLTZ_RUNTIME_BASE:-}" ]; then
    printf '%s\n' "$BOLTZ_RUNTIME_BASE"
  elif [ -n "${USER:-}" ] && [ -d "/scratch-cbe/users/$USER" ] && [ -w "/scratch-cbe/users/$USER" ]; then
    printf '/scratch-cbe/users/%s/.boltz_runtime\n' "$USER"
  elif [ -n "${HOME:-}" ] && [ -d "$HOME" ] && [ -w "$HOME" ]; then
    printf '%s/.cache/boltz_runtime\n' "$HOME"
  else
    printf '%s/boltz_runtime_%s\n' "${TMPDIR:-/tmp}" "${USER:-$(id -u)}"
  fi
}
RUNTIME_BASE="$(runtime_base_default)"
RUNTIME_HOME="$RUNTIME_BASE/home"
RUNTIME_TMP="$RUNTIME_BASE/tmp"
mkdir -p "$RUNTIME_HOME" "$RUNTIME_TMP"
TARGET_ROOT="${1:-$WORK_DIR}"
if [ -d "$TARGET_ROOT" ]; then
  TARGET_ROOT="$(cd "$TARGET_ROOT" && pwd -P)"
fi
DEFAULT_IMAGE="$ROOT_DIR/containers/current"
if [ -e "$DEFAULT_IMAGE" ]; then
  DEFAULT_IMAGE="$(readlink -f "$DEFAULT_IMAGE")"
fi
ANALYSIS_IMAGE="${BOLTZ_ANALYSIS_APPTAINER_IMAGE:-${BOLTZ_APPTAINER_IMAGE:-$DEFAULT_IMAGE}}"
export BOLTZ_CONTAINER_SITEPKGS="${BOLTZ_CONTAINER_SITEPKGS:-$ROOT_DIR/sitepkgs_bundle}"
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE
unset APPTAINERENV_PYTHONPATH APPTAINERENV_PYTHONHOME APPTAINERENV_PYTHONUSERBASE
export APPTAINERENV_PYTHONNOUSERSITE=1
export APPTAINERENV_TMPDIR="$RUNTIME_TMP"
export APPTAINERENV_TMP="$RUNTIME_TMP"
export APPTAINERENV_TEMP="$RUNTIME_TMP"

if [ -d "$BOLTZ_CONTAINER_SITEPKGS" ]; then
  export APPTAINERENV_PYTHONPATH="$BOLTZ_CONTAINER_SITEPKGS"
fi

exec apptainer exec --cleanenv --no-mount hostfs \
  --home "$RUNTIME_HOME" \
  --pwd "$WORK_DIR" \
  --bind "$ROOT_DIR:$ROOT_DIR" \
  --bind "$WORK_DIR:$WORK_DIR" \
  --bind "$TARGET_ROOT:$TARGET_ROOT" \
  --bind "$RUNTIME_TMP:$RUNTIME_TMP" \
  "$ANALYSIS_IMAGE" \
  /bin/bash --noprofile --norc -lc 'export PATH="/usr/local/apps/pyenv/versions/miniforge3-24.11.3-2/envs/boltz-conda/bin:$PATH"; exec python -I "$@"' \
  _ "$ROOT_DIR/boltz_analysis.py" "$@"
