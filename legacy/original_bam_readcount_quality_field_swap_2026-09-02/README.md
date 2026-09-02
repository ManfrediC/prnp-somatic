# Original bam-readcount quality-field swap

The original converter labelled bam-readcount quality fields in the wrong order.
It wrote `metrics[1]` as `MEAN_BQ` and `metrics[2]` as `MEAN_MQ`. The corrected
mapping is `metrics[1]` to `MEAN_MQ` and `metrics[2]` to `MEAN_BQ`, matching
`base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:...`.

The affected files are the 39 per-sample `*_metrics.tsv` tables, the 30 controls,
CJD and dilution variant-QC tables written by the downstream R stage, Table 6,
the combined tables PDF, and the provenance inventories that record the
converter and declared outputs. `MANIFEST.tsv` records all 73 preserved files'
original repository-relative paths, byte sizes, SHA-256 digests and modification
times.

The active repository revision before correction was
`621d2cf16da3fe1502e816766e5ca60645334135`. The working tree already contained
unrelated changes under `src2/sequencing2/`: six modified files and two untracked
files. Those changes were preserved and were not used for this correction.

```text
 M src2/sequencing2/12_cjd_dilutions_variant_qc_with_pon.sh
 M src2/sequencing2/1_controls_mutect2_no_pon.sh
 M src2/sequencing2/4_readcount_to_tsv.py
 M src2/sequencing2/5_controls_variant_qc_no_pon.sh
 M src2/sequencing2/README.md
 M src2/sequencing2/sequencing2.env.example
?? src2/sequencing2/archive_sequencing2_intermediates.ps1
?? src2/sequencing2/run_sequencing2.sh
```

The active `runs/` intermediates had already been removed. The erroneous metrics
were recovered from the retained original-pipeline VM backup at
`C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\VM backups\retired_add_disk_20260709_211229_files\prnp-somatic`.
They were copied into this legacy directory without altering the backup. The 30
active tracked QC outputs, Table 6 TeX source, combined tables PDF, script
inventory and final-output checksum inventory were moved here. The manifest was
verified in full before regeneration.

The current R wrappers were not used for the error-only replay because they now
apply an additional CJD AAF filter and add a display-label column. Instead, the
preserved QC tables supplied the exact original schemas, row order, membership
and filtering state; only `MEAN_BQ` and `MEAN_MQ` were replaced from the
corrected metrics by the original five-column join key.

Upstream caller and readcount outputs were not regenerated because the error was
introduced only when the existing raw bam-readcount records were converted to
labelled TSV columns. BAMs, raw readcounts, Mutect2 outputs, post-processing
outputs, annotations, PASS VCFs and VariantsToTable inputs do not depend on that
labelling step.
