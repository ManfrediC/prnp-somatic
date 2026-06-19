# Results

Generated workflow outputs live here. Most contents are gitignored and can be
rebuilt; selected manifests and lightweight review artefacts may be tracked.

## Main outputs

- ddPCR (`bash src/ddpcr/run_ddpcr.sh`):
  - `results/ddPCR/SNV_data_final.xlsx`
  - `results/ddPCR/SNV_pooled_participant.xlsx`
  - `results/ddPCR/p0_fallback.csv`

- SNV pipeline (`src/pipelines`):
  - `results/final_bam/*.bam` plus indexes
  - `results/mutect2_controls_no_pon/variant_qc/*`
  - `results/mutect2_controls_pon/panel_of_normals/*`
  - `results/mutect2_cjd_dilutions_with_pon/variant_qc/*`

- Junction workflow (`bash src/junctions/run_junctions.sh`):
  - `results/junctions/junction_align/*`
  - `results/junctions/junction_counts/prnp_junction_counts.tsv`
  - `results/junctions/junction_counts/prnp_junction_summary.tsv`
  - `results/junctions/junction_counts/prnp_junction_qc_table.tsv`

- Repeat workflow (`bash src/repeats/01_run_prnp_orr.sh`):
  - `results/repeats/sample_manifest.tsv`
  - `results/repeats/sample_calls.tsv`
  - `results/repeats/candidate_calls.tsv`
  - `results/repeats/cohort_summary.tsv`
  - `results/repeats/subclonal_read_support.tsv`
  - `results/repeats/gangstr_calls.tsv`
  - `results/repeats/somatic_screen.tsv`
  - `results/repeats/run_settings.tsv`

- QC helpers:
  - `results/sequencing_qc/sequencing_metrics_per_sample.tsv` from `make qc_metrics`
  - `results/dna_quality/sample_quality_evidence_table.tsv` from `make dna_quality`

Runtime intermediates belong under `runs/`.
