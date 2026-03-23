#!/bin/bash
# Nomogram one-click entry
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="$SCRIPT_DIR/scripts/Nomogram.R"

if [[ ! -f "$R_SCRIPT" ]]; then
  echo "Error: script not found: $R_SCRIPT" >&2
  exit 1
fi

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
Usage: ./run_nomogram.sh [Nomogram.R args]

Examples:
  ./run_nomogram.sh -d test_data/GSE126124.dat.csv -g test_data/GSE126124.group.csv -G test_data/04.gene.csv
  ./run_nomogram.sh -d data.csv -g group.csv -G genes.csv -c "Control,Tumor" -o results/Nomogram
EOF
  exit 0
fi

cd "$SCRIPT_DIR"
Rscript "$R_SCRIPT" "$@"
