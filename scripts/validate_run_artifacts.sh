#!/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: validate_run_artifacts.sh RUN_ROOT

Checks that a Boltz screening run produced the expected analysis artifacts:
  - summary_metrics.csv
  - analytics/slurm_metrics.csv
  - plots/scatter_dashboard.html
  - at least one plot in plots/dot_plots
  - at least one plot in plots/matrix_plots
  - at least one successful row in summary_metrics.csv
EOF
}

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
  usage
  exit 0
fi

if [ "$#" -ne 1 ]; then
  usage >&2
  exit 2
fi

ROOT="$1"
SUMMARY="$ROOT/summary_metrics.csv"
SLURM_METRICS="$ROOT/analytics/slurm_metrics.csv"
SCATTER="$ROOT/plots/scatter_dashboard.html"
DOT_DIR="$ROOT/plots/dot_plots"
MATRIX_DIR="$ROOT/plots/matrix_plots"

require_file() {
  local path="$1"
  if [ ! -f "$path" ]; then
    echo "missing file: $path" >&2
    exit 1
  fi
}

require_nonempty_dir() {
  local path="$1"
  if [ ! -d "$path" ]; then
    echo "missing directory: $path" >&2
    exit 1
  fi
  if ! find "$path" -maxdepth 1 -type f | grep -q .; then
    echo "no files found in: $path" >&2
    exit 1
  fi
}

require_file "$SUMMARY"
require_file "$SLURM_METRICS"
require_file "$SCATTER"
require_nonempty_dir "$DOT_DIR"
require_nonempty_dir "$MATRIX_DIR"

if ! awk -F, 'NR > 1 && $2 == "ok" { found = 1 } END { exit(found ? 0 : 1) }' "$SUMMARY"; then
  echo "no successful predictions recorded in: $SUMMARY" >&2
  exit 1
fi

echo "validation ok: $ROOT"
echo "  summary: $SUMMARY"
echo "  slurm metrics: $SLURM_METRICS"
echo "  scatter: $SCATTER"
echo "  dot plots: $(find "$DOT_DIR" -maxdepth 1 -type f | wc -l | tr -d ' ')"
echo "  matrix plots: $(find "$MATRIX_DIR" -maxdepth 1 -type f | wc -l | tr -d ' ')"
