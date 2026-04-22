#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Repository setup
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$REPO_ROOT"

# ---------------------------------------------------------------------------
# Python entrypoint resolution
# ---------------------------------------------------------------------------

PYTHON_BIN="${PYTHON_BIN:-python3}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  if command -v python >/dev/null 2>&1; then
    PYTHON_BIN="python"
  else
    echo "ERROR: python3 or python is required on PATH." >&2
    exit 1
  fi
fi

# Hand off all workflow arguments to the Python worker
exec "$PYTHON_BIN" "$SCRIPT_DIR/02_build_dna_quality_tables.py" "$@"
