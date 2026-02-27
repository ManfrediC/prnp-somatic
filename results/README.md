# Results

This directory contains generated outputs created by running the reproducible workflows.

Most contents under `results/` are intentionally untracked in git and can be regenerated.
Only placeholders and selected manifests are kept in version control.

## What you will find after running workflows

- ddPCR workflow (`bash src/ddPCR/run_ddpcr.sh`):
  - `results/ddPCR/SNV_data_final.xlsx`
  - `results/ddPCR/SNV_pooled_participant.xlsx`
  - `results/ddPCR/p0_fallback.csv`

- SNV pipeline (`src/pipelines`):
  - `results/final_bam/*.bam` (+ indexes)
  - `results/mutect2_controls_no_pon/variant_qc/*`
  - `results/mutect2_controls_pon/panel_of_normals/*`
  - `results/mutect2_cjd_dilutions_with_pon/variant_qc/*`

- Junction workflow (`bash src/junctions/run_junctions.sh`):
  - `results/junctions/junction_align/*`
  - `results/junctions/junction_counts/prnp_junction_counts.tsv`
  - `results/junctions/junction_counts/prnp_junction_summary.tsv`

## Notes

- Outputs in this directory are derived artefacts and may be deleted/rebuilt.
- Runtime intermediates are primarily written under `runs/`; final and summary outputs are written here.
