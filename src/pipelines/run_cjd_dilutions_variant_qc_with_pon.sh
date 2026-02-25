#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$REPO_ROOT"

# Publication-path entrypoint for CJD+dilutions variant table/QC regeneration.
bash "$SCRIPT_DIR/12_cjd_dilutions_variant_qc_with_pon.sh" "$@"
