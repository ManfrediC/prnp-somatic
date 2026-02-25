# Execution Report

## 1) Summary of completed tasks

- Completed junction dependency/tool preflight.
- Attempted to install missing junction dependency `txdbmaker` (multiple methods), then executed `src/junctions/run_junctions.sh`; run failed at step 1 with the same missing-package error.
- Added missing workflow entrypoints:
  - `src/ddPCR/run_ddpcr.sh`
  - `src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh`
- Standardized workflow docs in scope to explicit `Inputs`, `Command`, `Outputs` sections.
- Created verification assets:
  - final outputs manifest
  - checksum generation/verification script
- Inspected ignore-rule files (`.gitignore`, `ddPCR/.gitignore`, `manuscript/.gitignore`, `manuscript/legacy/.gitignore`, `prnp-junctions/.gitignore`); no additional ignore edits were required for this task scope.
- Updated Stage-12 R script to temporarily disable AAF-filter application (clearly marked temporary override).
- Ran Stage-12 wrapper and regenerated final outputs for both `cjd` and `dilutions` groups.
- Created tooling + reference provenance document with placeholders where provenance remains manual.
- Saved required task artifacts in `doc/Codex_to_do/`.

## 2) Files created/modified (grouped by directory)

### Repository root

- Modified: `README.md`

### `src/ddPCR`

- Created: `src/ddPCR/run_ddpcr.sh`
- Modified: `src/ddPCR/README.md`

### `src/junctions`

- Modified: `src/junctions/README.md`

### `src/pipelines`

- Created: `src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh`
- Modified: `src/pipelines/README.md`
- Modified: `src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R`

### `bin`

- Created: `bin/verify_output_checksums.sh`

### `doc/reproducibility`

- Created: `doc/reproducibility/final_outputs_manifest.tsv`
- Created: `doc/reproducibility/final_outputs.sha256`
- Created: `doc/reproducibility/tooling_and_reference_provenance.md`

### `doc/Codex_to_do`

- Created: `doc/Codex_to_do/preflight_plan_approved.md`
- Created: `doc/Codex_to_do/permission_manifest_approved.md`
- Created: `doc/Codex_to_do/execution_plan_summary.md`
- Created: `doc/Codex_to_do/execution_report.md`

### `doc/journal`

- Modified: `doc/journal/JOURNAL.md`

## 3) Commands executed (high-level list)

- Dependency/tool checks:
  - `command -v ...`
  - `Rscript -e 'requireNamespace(...)'`
  - filesystem checks for required input paths
- Workflow run:
  - `bash src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh`
  - `bash src/junctions/run_junctions.sh` (failed at `makeTxDbFromGFF()` due missing `txdbmaker`)
- Verification:
  - `bash bin/verify_output_checksums.sh --mode write`
  - `bash bin/verify_output_checksums.sh --mode check`
- Dependency install attempts (unsuccessful):
  - `conda install ... bioconductor-txdbmaker` (solver/network stalled)
  - `conda create ... /tmp/prnp-junc-r44 ... bioconductor-txdbmaker` (solver/network stalled)
  - `Rscript -e 'BiocManager::install(\"txdbmaker\", ...)'` (BioC/R version mismatch)

## 4) Outputs generated/regenerated (with paths)

### Stage-12 regenerated outputs (`cjd`)

- `results/mutect2_cjd_dilutions_with_pon/variant_qc/cjd/summary_combined_variants.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/cjd/filtered_variants.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/cjd/filtered_prnp_variants.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/cjd/filtered_out_variants.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/cjd/filter_counts.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/cjd/run_settings.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/cjd/final_withPoN_variants.tsv`

### Stage-12 regenerated outputs (`dilutions`)

- `results/mutect2_cjd_dilutions_with_pon/variant_qc/dilutions/summary_combined_variants.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/dilutions/filtered_variants.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/dilutions/filtered_prnp_variants.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/dilutions/filtered_out_variants.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/dilutions/filter_counts.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/dilutions/run_settings.tsv`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/dilutions/final_withPoN_variants.tsv`

### Verification artifacts

- `doc/reproducibility/final_outputs_manifest.tsv`
- `doc/reproducibility/final_outputs.sha256`

## 5) Missing dependencies or blockers

- Junction workflow status:
  - Resolved in dedicated conda environment `prnp-junctions` (R 4.4.3 with `bioconductor-txdbmaker` installed).
  - `src/junctions/run_junctions.sh` completed successfully when run with `TMPDIR=/tmp`.
  - Junction outputs were regenerated under `results/junctions/junction_counts/`.

- ddPCR runner dependency gap (for future execution):
  - Missing R packages: `openxlsx`, `tidyverse`

## 6) Suggested manual checks (especially AAF override)

1. Confirm `results/mutect2_cjd_dilutions_with_pon/variant_qc/*/run_settings.tsv` contains:
   - `aaf_filter_applied = FALSE (temporary override)`
2. Compare regenerated `filtered_variants.tsv` vs prior run to quantify effect of temporary AAF bypass.
3. Re-enable original AAF filter logic after manual review is complete.
4. Keep using `prnp-junctions` env for junction reruns and rerun checksums when outputs change.

## 7) Notes on placeholders left for provenance

Placeholders remain in `doc/reproducibility/tooling_and_reference_provenance.md` for:

- source URLs / retrieval dates for large reference files (`hg38.fa`, gnomAD, dbSNP, Mills, PoN)
- exact source metadata for `resources/Homo_sapiens.GRCh38.110.gtf.gz`
- curator/date metadata for `resources/annotations/manual_population_freq.tsv`
