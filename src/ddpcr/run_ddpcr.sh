#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$REPO_ROOT"

log_dir="results/ddPCR/logs"
mkdir -p "$log_dir"

run_r() {
  local script_path="$1"
  local log_name
  log_name="$(basename "$script_path" .R).log"
  echo "==> Rscript $script_path"
  Rscript "$script_path" 2>&1 | tee "$log_dir/$log_name"
}

run_py() {
  local script_path="$1"
  local log_name
  log_name="$(basename "$script_path" .py).log"
  echo "==> python $script_path"
  python "$script_path" 2>&1 | tee "$log_dir/$log_name"
}

run_r "$SCRIPT_DIR/create_snv_dataframe.R"
run_r "$SCRIPT_DIR/figures/ddpcr_fractional_abundance.R"
run_py "$SCRIPT_DIR/figures/fix_snv_all_mutations_legend_bottom_final_svg.py"
run_r "$SCRIPT_DIR/figures/ddpcr_fractional_abundance_pooled.R"
run_r "$SCRIPT_DIR/ddpcr_samples_results_tbl.R"
run_r "$SCRIPT_DIR/ddpcr_sample_number.R"
run_r "$SCRIPT_DIR/estimate_haploid_genomes_surveyed.R"
run_r "$SCRIPT_DIR/figures/create_ddpcr_scatterplots.R"
run_py "$SCRIPT_DIR/figures/build_ddpcr_gating_svg_panels.py"
