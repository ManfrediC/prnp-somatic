#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$REPO_ROOT"

# Run the analytical stages in the order documented for the original pipeline.
bash "$SCRIPT_DIR/1_controls_mutect2_no_pon.sh"
bash "$SCRIPT_DIR/2_controls_postprocess_no_pon.sh"
bash "$SCRIPT_DIR/3_controls_readcount_qc_no_pon.sh"
bash "$SCRIPT_DIR/5_controls_variant_qc_no_pon.sh"
bash "$SCRIPT_DIR/7_controls_create_pon.sh"
bash "$SCRIPT_DIR/8_cjd_dilutions_mutect2_with_pon.sh"
bash "$SCRIPT_DIR/9_cjd_dilutions_postprocess_with_pon.sh"
bash "$SCRIPT_DIR/10_cjd_dilutions_readcount_qc_with_pon.sh"
bash "$SCRIPT_DIR/11_cjd_dilutions_readcount_to_tsv_with_pon.sh"
bash "$SCRIPT_DIR/12_cjd_dilutions_variant_qc_with_pon.sh"
