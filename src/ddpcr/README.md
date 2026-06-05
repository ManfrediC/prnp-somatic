# ddPCR Workflow

This folder contains the official *PRNP* somatic-mutation ddPCR workflow.
It rebuilds the ddPCR analysis from immutable source exports under
`raw/ddpcr`, including `.ddpcr` archives, exported droplet-amplitude JSONs,
canonical CSV exports, layout workbooks, and raw provenance manifests.

## Scripts

- `create_snv_dataframe.R`
- `ddpcr_raw_import_helpers.R`
- `ddpcr_fractional_abundance.R`
- `ddpcr_fractional_abundance_pooled.R`
- `ddpcr_samples_results_tbl.R`
- `ddpcr_sample_number.R`
- `estimate_haploid_genomes_surveyed.R`
- `create_ddpcr_scatterplots.R`

## Inputs

Expected under the repository root:

- Raw ddPCR database: `raw/ddpcr/{ddpcr_archive,csv_export,archive_contents,layout_xlsx,manifests}/`
- Sample metadata: `raw/ddpcr/sample_details.xlsx`

## Command

Create a dedicated conda environment (reviewer-side). Then the workflow command from repository root is:

```bash
bash src/ddpcr/run_ddpcr.sh
```

Environment setup example:

```bash
conda env create -f env/ddpcr.environment.yml
conda activate prnp-somatic-ddpcr
```

Windows/RStudio usage:

1. Open this repository in RStudio, or set the working directory anywhere inside the repository.
2. From the repository root, run `source("src/ddpcr/create_snv_dataframe.R")`, or open the file and click Source in RStudio.
3. Then run the downstream scripts in the order listed in `run_ddpcr.sh`.

Optional Docker path (if Docker daemon is available):

```bash
docker run --rm -v "$(pwd)":/work -w /work rocker/tidyverse:4.4.3 bash -lc \
  "Rscript -e 'options(repos=c(CRAN=\"https://cloud.r-project.org\")); install.packages(c(\"openxlsx\",\"binom\",\"readxl\",\"jsonlite\"))' && bash src/ddpcr/run_ddpcr.sh"
```

## Outputs

Written to `results/ddPCR/`:

- `SNV_data_final.xlsx` (long format)
- `SNV_pooled_participant.xlsx`
- `p0_fallback.csv`
- `ddpcr_partition_counts_by_sample_assay.csv`
- `validation/`
- `scatterplots/`
- `haploid_genomes_surveyed/`

The workflow also updates ddPCR-related manuscript tables and figures under
`manuscript/tables` and `manuscript/figures`.

### Output File Meanings

- `SNV_data_final.xlsx`:
  - main long-format analysis table produced by the ddPCR workflow (sample/assay-level SNV results used for downstream interpretation).

- `SNV_pooled_participant.xlsx`:
  - participant-level pooled summary derived from the long-format table (aggregates replicate/sample-level rows into one consolidated participant view where applicable).

- `p0_fallback.csv`:
  - fallback blank-rate (`p0`) values used when a plate lacks usable blank controls, supporting LoB/LoBFA-based background classification.

- `ddpcr_partition_counts_by_sample_assay.csv`:
  - intermediate sample/assay partition-count table used by the haploid-genome supplementary workflow.

- `validation/`:
  - raw import checks, old-versus-new comparison outputs where prior results were archived, LoB sensitivity outputs, and scatterplot manifests.

- `scatterplots/`:
  - rendered droplet-amplitude scatterplots from exported raw JSON data.

- `haploid_genomes_surveyed/`:
  - dedicated supplementary output directory containing sample-region, participant-pooled, participant-review TeX, and run-summary outputs for droplet counts and Poisson-corrected haploid genome-equivalent estimates.
  - cohort-level pooled droplet summaries are not emitted.

These expected output paths are also listed in:

- `doc/reproducibility/final_outputs_manifest.tsv`

## R packages

The script uses:

- `readr`
- `readxl`
- `jsonlite`
- `ggplot2`
- `tidyverse`
- `openxlsx`
- `magrittr`
- `binom`
