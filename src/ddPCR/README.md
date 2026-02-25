# ddPCR reproducible script

This folder contains a portable version of the legacy ddPCR dataframe workflow.

## Script

- `create_snv_dataframe.R`

This script is intentionally kept close to the legacy structure and annotations, with minimal changes for reproducibility:
- no hardcoded Windows `setwd(...)` paths,
- repo-relative input/output paths,
- metadata input renamed to `sample_details.xlsx`.

## Inputs

Expected in the repository root:

- Raw ddPCR CSV files: `ddPCR/*.csv`
- Sample metadata: `ddPCR/sample_details.xlsx`

## Command

Create a dedicated conda environment (reviewer-side), then run from the repository root:

```bash
conda --no-plugins create -n prnp-somatic-ddpcr -c conda-forge -y \
  r-base=4.3 r-readr r-tidyverse r-openxlsx r-magrittr r-binom
conda activate prnp-somatic-ddpcr
bash src/ddPCR/run_ddpcr.sh
```

This repository does not provide or version-control the R tools themselves.

Optional Docker path (if Docker daemon is available):

```bash
docker run --rm -v "$(pwd)":/work -w /work rocker/tidyverse:4.3.3 bash -lc \
  "Rscript -e 'options(repos=c(CRAN=\"https://cloud.r-project.org\")); install.packages(c(\"openxlsx\",\"binom\"))' && bash src/ddPCR/run_ddpcr.sh"
```

## Outputs

Written to `results/ddPCR/`:

- `SNV_data_final.xlsx` (long format)
- `SNV_pooled_participant.xlsx`
- `p0_fallback.csv`

These expected output paths are also listed in:

- `doc/reproducibility/final_outputs_manifest.tsv`

## R packages

The script uses:

- `readr`
- `tidyverse`
- `openxlsx`
- `magrittr`
- `binom`
