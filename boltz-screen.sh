#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_IMAGE="$ROOT_DIR/containers/current"
if [ -e "$DEFAULT_IMAGE" ]; then
  DEFAULT_IMAGE="$(readlink -f "$DEFAULT_IMAGE")"
fi
export BOLTZ_APPTAINER_IMAGE="${BOLTZ_APPTAINER_IMAGE:-$DEFAULT_IMAGE}"
export BOLTZ_CONTAINER_SITEPKGS="${BOLTZ_CONTAINER_SITEPKGS:-$ROOT_DIR/sitepkgs_bundle}"
TMP_LOG="$(mktemp "${TMPDIR:-/tmp}/boltz-screen.XXXXXX.log")"
trap 'rm -f "$TMP_LOG"' EXIT

if [ ! -e "$BOLTZ_APPTAINER_IMAGE" ]; then
  echo "Container image not found: $BOLTZ_APPTAINER_IMAGE" >&2
  exit 1
fi

env BOLTZ_PREPARE_ONLY=1 \
  apptainer exec --cleanenv "$BOLTZ_APPTAINER_IMAGE" \
  python "$ROOT_DIR/boltz-2_wrapper.py" "$@" | tee "$TMP_LOG"

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
