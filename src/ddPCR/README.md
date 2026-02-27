# ddPCR reproducible script

This folder contains the *PRNP* somatic mutations ddPCR dataframe workflow.

## Script

- `create_snv_dataframe.R`

## Inputs

Expected in the repository root:

- Raw ddPCR CSV files: `ddPCR/*.csv`
- Sample metadata: `ddPCR/sample_details.xlsx`

## Command

Create a dedicated conda environment (reviewer-side). Then the workflow command from repository root is:

```bash
bash src/ddPCR/run_ddpcr.sh
```

Environment setup example:

```bash
conda env create -f env/ddpcr.environment.yml
conda activate prnp-somatic-ddpcr
```

Optional Docker path (if Docker daemon is available):

```bash
docker run --rm -v "$(pwd)":/work -w /work rocker/tidyverse:4.4.3 bash -lc \
  "Rscript -e 'options(repos=c(CRAN=\"https://cloud.r-project.org\")); install.packages(c(\"openxlsx\",\"binom\"))' && bash src/ddPCR/run_ddpcr.sh"
```

## Outputs

Written to `results/ddPCR/`:

- `SNV_data_final.xlsx` (long format)
- `SNV_pooled_participant.xlsx`
- `p0_fallback.csv`

### Output File Meanings

- `SNV_data_final.xlsx`:
  - main long-format analysis table produced by the ddPCR workflow (sample/assay-level SNV results used for downstream interpretation).

- `SNV_pooled_participant.xlsx`:
  - participant-level pooled summary derived from the long-format table (aggregates replicate/sample-level rows into one consolidated participant view where applicable).

- `p0_fallback.csv`:
  - fallback blank-rate (`p0`) values used when a plate lacks usable blank controls, supporting LoB/LoBFA-based background classification.

These expected output paths are also listed in:

- `doc/reproducibility/final_outputs_manifest.tsv`

## R packages

The script uses:

- `readr`
- `tidyverse`
- `openxlsx`
- `magrittr`
- `binom`
