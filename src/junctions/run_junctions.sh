#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$REPO_ROOT"

CONFIG_FILE="${CONFIG_FILE:-${REPO_ROOT}/config/junctions.env}"
if [ -f "$CONFIG_FILE" ]; then
  set -a
  # shellcheck source=/dev/null
  source "$CONFIG_FILE"
  set +a
fi

echo "[1/4] Build PRNP junction FASTA"
Rscript "$SCRIPT_DIR/01_build_prnp_junction_fasta.R"

echo "[2/4] Build padded PRNP BED"
Rscript "$SCRIPT_DIR/02_make_prnp_bed.R"

echo "[3/4] Extract + realign PRNP reads to junction reference"
bash "$SCRIPT_DIR/03_process_bam.sh"

echo "[4/4] Count junction-spanning fragments"
Rscript "$SCRIPT_DIR/04_count_prnp_junctions.R"

echo "Done."
