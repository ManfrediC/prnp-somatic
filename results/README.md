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

- Sequencing QC metrics workflow:
  - `results/qc/<QC_RUN>/sequencing_metrics_per_sample.tsv` (current Makefile path via `make qc_metrics`)
  - producing script: `authoritative_files/compute_sequencing_metrics.py`

- DNA quality proxy workflow:
  - `results/dna_quality/<RUN_ID>/file_inventory.tsv`
  - `results/dna_quality/<RUN_ID>/library_qc.tsv`
  - `results/dna_quality/<RUN_ID>/sample_quality_master.tsv`
  - `results/dna_quality/<RUN_ID>/dna_quality_scorecard.tsv`
  - producing script: `src/dna_quality/run_dna_quality_analysis.ps1`

## Notes

- Outputs in this directory are derived artefacts and may be deleted/rebuilt.
- Runtime intermediates are primarily written under `runs/`; final and summary outputs are written here.
