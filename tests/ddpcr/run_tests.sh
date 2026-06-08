#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

Rscript tests/ddpcr/test_fractional_abundance.R
Rscript tests/ddpcr/test_dpcp_input_helpers.R
Rscript tests/ddpcr/test_umbrella_input_helpers.R
