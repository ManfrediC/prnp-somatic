# Sequencing2 numerical provenance

This audit is isolated from the historical numerical-consistency audit. It does not replace or modify the original record.

## Empirical LoD

`results2/spikein/read_recovery/empirical_lod.tsv`

-> `26/3891 = 0.006682086867129272`, with `FilterMutectCalls` status `PASS`

-> `src2/sequencing2/variant_table_qc_common.R`

-> CJD and control AAF filtering recorded in each sequencing2 `run_settings.tsv`

## Final CJD calls

`results2/sequencing2/results/mutect2_cjd_dilutions_with_pon/variant_qc/cjd/filtered_prnp_variants.tsv`

-> `manuscript2/tables2/make_prnp_somatic_snv_summary_tex.py`

-> `manuscript2/tables2/table_prnp_somatic_snv_summary.tex` and `.pdf`

The same CJD TSV also feeds:

-> `manuscript2/figures2/snv_lollipop/make_snv_lollipop.R`

-> `manuscript2/figures2/snv_lollipop/build_snv_lollipop_svg.py`

-> the isolated copies of the original R/trackViewer baseline and reviewed SVG
template

-> `SNV_lollipop_data.csv`, `.svg`, `.pdf`, edit-diff JSON and the isolated TeX
wrapper

## Spike-in LoD table

The marker definition, direct read recovery and final exact-allele calls under `results2/spikein/` feed:

-> `manuscript2/tables2/make_spikein_snp_recovery_table.R`

-> `manuscript2/tables2/table_spikein_snp_recovery.tex` and `.pdf`

## Isolation boundary

No file under `manuscript/`, `audit/`, `reproducibility/`, `tests/regression/`, `src/` or `results/` is generated or updated by this workflow.
