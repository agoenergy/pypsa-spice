#!/usr/bin/env bash
# Before running the script, make sure that you edit scenarios list in the `compile_csvs_into_excel.py` to ensure it includes the desired scenarios.
set -euo pipefail

# Repo root = directory containing this bash file
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$REPO_ROOT/data/mod-data-th-pdp/TH_PDP_release"

# Run from script directory so relative paths in the Python script work there
(
  cd "$SCRIPT_DIR"
  python3 compile_csvs_into_excel.py
)
