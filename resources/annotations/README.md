This directory stores small, versioned annotation helper tables used by downstream QC steps.

## Files

- `manual_population_freq.tsv`
  - Manually curated population frequencies used by the controls no-PoN QC step.
  - Values were transcribed from the historical Windows analysis code and journal notes.
  - The table includes both `dbsnp_id`-based and `variant_id`-based lookups.

## Notes

- `population_frequency` is treated as missing (`NA`) when blank.
- For variants with no reported population frequency in the original curation, the field is intentionally left blank.
