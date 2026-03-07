# Controls No-PoN Results

This directory stores final result tables for the controls-only Mutect2 no-PoN workflow.

## `variant_qc/`

- `summary_combined_variants.tsv`: combined per-variant table before final QC filtering (includes annotation fields, readcount metrics, and population-frequency fields).
- `filtered_variants.tsv`: final QC-passing variants after all configured filters.
- `filtered_prnp_variants.tsv`: subset of `filtered_variants.tsv` restricted to `PRNP`.
- `filtered_out_variants.tsv`: variants removed by one or more QC filters.
- `filter_counts.tsv`: stepwise row counts (`before`, `after`, `removed`) for each QC criterion.
- `run_settings.tsv`: threshold/settings used for this QC run (for reproducibility).

Legacy-style convenience exports (kept for continuity with older analyses):

- `final_noPoN_variants.tsv`: same final set as `filtered_variants.tsv`.
- `noPoN_PRNP_PASS.tsv`: pre-final-filter `PRNP` subset from the combined table.
- `noPoN_TET2_PASS.tsv`: pre-final-filter `TET2` subset from the combined table.
- `noPoN_TTN_PASS.tsv`: pre-final-filter `TTN` subset from the combined table.
