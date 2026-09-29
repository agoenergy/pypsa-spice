#!/usr/bin/env bash
# Before running the script, make sure that you edit SELECT_SCENARIO in the `generate_hourly_gen_tables_for_flourish.py` to ensure it includes the desired scenarios.
set -euo pipefail

# Repo root = directory containing this bash file
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$REPO_ROOT/data/mod-data-th-pdp/TH_PDP_release"

# Run from script directory so relative paths in the Python script work there
(
  cd "$SCRIPT_DIR"
  python3 generate_hourly_gen_tables_for_flourish.py
)
