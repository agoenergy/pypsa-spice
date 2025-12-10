# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

# coding: utf-8

#!/usr/bin/env bash
# render_dag.sh — render a Snakemake job DAG to SVG/PNG
#
# What it does:
#   - Runs `snakemake -n --dag` to emit a DOT graph of job dependencies.
#   - Trims any preamble before 'digraph' (some Snakemake versions print info lines).
#   - Pipes DOT into Graphviz `dot` to produce an image file.
#
# Output:
#   - The image is ALWAYS written to a directory named $OUTPUT_DIR_NAME
#     located next to this script. The directory is created if missing.
#
# Prerequisites:
#   - bash >= 4.0 (for lowercase ${var,,} parameter expansion)
#   - snakemake (to generate the DAG in DOT format)
#   - graphviz (provides `dot` for rendering the DAG to SVG/PNG)
#
# Examples:
#   ./render_dag.sh
#   ./render_dag.sh -f png -l portrait -o my_run
#   ./render_dag.sh --check

# exit on errors, exit on unset vars, and fail pipelines if any step fails.
set -euo pipefail

FORMAT="png"                # svg|png (output format)
LAYOUT="portrait"           # landscape|portrait (graph orientation)
OUTPUT_DIR_NAME="output"    # folder placed next to this script (created if missing)
OUTPUT_NAME=""              # base name (no extension). Defaults to 'pypsa-spice_dag'

usage() {
  cat <<EOF

Render Snakemake job DAG 🐍
===========================

Options:
  -c, --check    Check prerequisites
  -f, --format   svg|png             (default: png)
  -l, --layout   landscape|portrait  (default: portrait)
  -o, --output   Output base name    (default: pypsa-spice_dag)
  -h, --help     Show this help
EOF
}

check_prereqs() {
  echo 
  echo "🔍 Checking prerequisites..."
  echo "============================="
  echo "OS:        $(uname -o)"
  echo "Bash:      $BASH_VERSION"
  if command -v snakemake >/dev/null 2>&1; then
    echo "Snakemake: $(snakemake --version  | grep -v 'Restricted license')"
  else
    echo "Snakemake: MISSING ❌"
  fi
  if command -v dot >/dev/null 2>&1; then
    echo "Graphviz:  $(dot -V 2>&1)"
  else
    echo "Graphviz:  MISSING ❌"
  fi
}

# Parse args
while [[ $# -gt 0 ]]; do
  case "${1:-}" in
    -f|--format)   FORMAT="${2,,}"; shift 2 ;;
    -l|--layout)   LAYOUT="${2,,}"; shift 2 ;;
    -o|--output)   OUTPUT_NAME="$2"; shift 2 ;;
    -c|--check)    check_prereqs; exit 0 ;;
    -h|--help)     usage; exit 0 ;;
    *) echo; echo "Unknown Option: $1 ❌"; usage; exit 1 ;;
  esac
done

# Directory of this script
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"

# Dependency checks
command -v snakemake >/dev/null 2>&1 || { echo "Error: snakemake not found" >&2; exit 1; }
command -v dot       >/dev/null 2>&1 || { echo "Error: graphviz 'dot' not found" >&2; exit 1; }

# Normalize options
case "$FORMAT" in svg|png) : ;; *) echo "Error: --format must be svg or png"; exit 1 ;; esac
case "$LAYOUT" in
  landscape|land|lr) RANKDIR="LR" ;;   # Left-to-Right
  portrait|port|tb)  RANKDIR="TB" ;;   # Top-to-Bottom
  *) echo "Error: --layout must be landscape (lr) or portrait (tb)"; exit 1 ;;
esac

# Defaults
[[ -z "$OUTPUT_NAME" ]] && OUTPUT_NAME="pypsa-spice_dag"

# Remove everything from the first '.' onwards and add a timestamp
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
OUTPUT_BASE="${TIMESTAMP}_${OUTPUT_NAME%%.*}"

# Enforced output directory next to the script
OUT_DIR="${SCRIPT_DIR}/${OUTPUT_DIR_NAME}"
mkdir -p -- "$OUT_DIR"

# Final output path
OUT_PATH="${OUT_DIR}/${OUTPUT_BASE}.${FORMAT}"

# Render:
PYTHONWARNINGS=ignore snakemake -n --dag \
  | awk 'BEGIN{p=0} /^digraph/{p=1} p' \
  | dot -Grankdir="$RANKDIR" -T"$FORMAT" > "$OUT_PATH"

echo "🎉 DAG created in: '$OUT_PATH'"
