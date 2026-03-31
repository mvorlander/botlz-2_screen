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
export APPTAINERENV_TMPDIR="$RUNTIME_TMP"
export APPTAINERENV_TMP="$RUNTIME_TMP"
export APPTAINERENV_TEMP="$RUNTIME_TMP"
export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1
if [ -d "$BOLTZ_CONTAINER_SITEPKGS" ]; then
  export APPTAINERENV_PYTHONPATH="$BOLTZ_CONTAINER_SITEPKGS"
fi

declare -a EXTRA_BINDS=()
declare -A EXTRA_BIND_SEEN=()
add_bind_target() {
  local raw="${1:-}"
  local bind_path=""
  [ -n "$raw" ] || return 0
  case "$raw" in
    /*)
      if [ -d "$raw" ]; then
        bind_path="$(cd "$raw" && pwd -P)"
      elif [ -e "$raw" ]; then
        bind_path="$(cd "$(dirname "$raw")" && pwd -P)"
      else
        bind_path="$raw"
        while [ ! -e "$bind_path" ] && [ "$bind_path" != "/" ]; do
          bind_path="$(dirname "$bind_path")"
        done
        [ -d "$bind_path" ] || return 0
        bind_path="$(cd "$bind_path" && pwd -P)"
      fi
      if [ -z "${EXTRA_BIND_SEEN[$bind_path]+x}" ]; then
        EXTRA_BIND_SEEN["$bind_path"]=1
        EXTRA_BINDS+=(--bind "$bind_path:$bind_path")
      fi
      ;;
  esac
}
add_bind_value() {
  local raw="${1:-}"
  local item
  local IFS=','
  for item in $raw; do
    add_bind_target "$item"
  done
}
path_arg=""
for arg in "$@"; do
  if [ -n "$path_arg" ]; then
    add_bind_value "$arg"
    path_arg=""
    continue
  fi
  case "$arg" in
    --screen|--bait|-c|--chain-names|--out|-o|--msa|-m)
      path_arg="$arg"
      ;;
    --screen=*|--bait=*|--chain-names=*|--out=*|--msa=*|-c=*|-o=*|-m=*)
      add_bind_value "${arg#*=}"
      ;;
    /*)
      add_bind_target "$arg"
      ;;
  esac
done

TMP_LOG="$(mktemp "$RUNTIME_TMP/boltz-screen.XXXXXX.log")"
trap 'rm -f "$TMP_LOG"' EXIT

if [ ! -e "$PREP_IMAGE" ]; then
  echo "Preparation container not found: $PREP_IMAGE" >&2
  exit 1
fi

if [ ! -e "$BOLTZ_APPTAINER_IMAGE" ]; then
  echo "Container image not found: $BOLTZ_APPTAINER_IMAGE" >&2
  exit 1
fi

cmd=(
  env
  APPTAINERENV_BOLTZ_PREPARE_ONLY=1
  APPTAINERENV_BOLTZ_APPTAINER_IMAGE="$BOLTZ_APPTAINER_IMAGE"
  APPTAINERENV_BOLTZ_CONTAINER_SITEPKGS="$BOLTZ_CONTAINER_SITEPKGS"
  apptainer exec --cleanenv --no-mount hostfs
  --home "$RUNTIME_HOME"
  --pwd "$WORK_DIR"
  --bind "$ROOT_DIR:$ROOT_DIR"
  --bind "$WORK_DIR:$WORK_DIR"
)
if [ "${#EXTRA_BINDS[@]}" -gt 0 ]; then
  cmd+=("${EXTRA_BINDS[@]}")
fi
cmd+=(
  --bind "$RUNTIME_TMP:$RUNTIME_TMP"
  "$PREP_IMAGE"
  /bin/bash --noprofile --norc -lc 'export PATH="/usr/local/apps/pyenv/versions/miniforge3-24.11.3-2/envs/boltz-conda/bin:$PATH"; exec python -I "$@"'
  _
  "$ROOT_DIR/boltz-2_wrapper.py"
  "$@"
)
"${cmd[@]}" | tee "$TMP_LOG"

DISPATCH_ENV="$(awk -F= '/^BOLTZ_DISPATCH_ENV=/{print $2}' "$TMP_LOG" | tail -n 1)"
if [ -z "$DISPATCH_ENV" ]; then
  exit 0
fi

# shellcheck disable=SC1090
source "$DISPATCH_ENV"

ARRAY_RES="$(sbatch "$ARRAY_SLURM")"
ARRAY_ID="${ARRAY_RES##* }"
echo "🗄  array job $ARRAY_ID (${JOB_COUNT} tasks)"

ANALYSIS_SUBMIT="${ANALYSIS_SLURM%.slurm}.submit.slurm"
sed "s/__ARRAY_ID__/$ARRAY_ID/g" "$ANALYSIS_SLURM" > "$ANALYSIS_SUBMIT"
sbatch "$ANALYSIS_SUBMIT" >/dev/null

echo "✅  Screening dispatch complete."
