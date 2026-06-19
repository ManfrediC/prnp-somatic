# Runs

Pipeline intermediates and run-level artefacts are written here. Final summary
outputs belong under `results/`.

Expected top-level workflow directories:

- `preprocessing/`
- `mutect2_controls_no_pon/`
- `mutect2_cjd_dilutions_with_pon/`

## `preprocessing/`

Per-sample outputs from `src/pipelines/preprocessing.sh`, organised as:

`runs/preprocessing/<batch_id>/<sample_id>/`

Each sample directory contains trimmed FASTQs, alignment/preprocessing BAMs,
Picard/BQSR artefacts, `RUN_META.txt`, and logs. Batch IDs come from
`config/preprocessing_samples.tsv`.

## Mutect2 run directories

`mutect2_controls_no_pon/` stores control-only tumour-only Mutect2 runs without
a panel of normals:

- `vcf/`
- `f1r2/`
- downstream stage folders such as `orientation/`, `filtered/`, `annot/`,
  `annot_with_gnomad/`, `readcount_qc/`, and `variant_qc/`

`mutect2_cjd_dilutions_with_pon/` stores PoN-enabled CJD and dilution runs:

- `cjd/<stage>/`
- `dilutions/<stage>/`
- `logs/`
