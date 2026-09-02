# Legacy and corrected output comparison

The corrected converter was checked against every allele block in all 39 raw
bam-readcount files retained in the original-pipeline VM backup. Regenerating
the metrics in a temporary directory produced byte-identical copies of all 39
active metrics tables. All 520 rows assign `avg_mapping_quality` to `MEAN_MQ`
and `avg_basequality` to `MEAN_BQ`.

## Metrics tables

| Group | Files | Data rows | Files with value changes | Rows with value changes |
| --- | ---: | ---: | ---: | ---: |
| Controls | 6 | 62 | 5 | 25 |
| CJD | 26 | 291 | 16 | 120 |
| Dilutions | 7 | 167 | 5 | 64 |
| Total | 39 | 520 | 26 | 209 |

Exactly 209 `MEAN_BQ` cells and 209 `MEAN_MQ` cells changed. No other metrics
column changed. Thirteen metrics files contain only equal BQ/MQ values and are
therefore byte-identical after regeneration.

## Error-only QC replay

The 30 preserved QC tables supplied the original schemas, row order, filtering
state and membership. Corrected quality values were joined by
`(sample_name, CHROM, POS, REF, ALT)` after uppercasing the alleles, matching the
original R join. Only exact non-zero metric matches were accepted. Blank quality
values were preserved only where no exact match existed.

| Result | Count |
| --- | ---: |
| QC files reconstructed | 30 |
| Quality-bearing QC files | 24 |
| Files with corrected values | 14 |
| Byte-identical files | 16 |
| Repeated row occurrences changed | 101 |
| `MEAN_BQ` cells changed | 101 |
| `MEAN_MQ` cells changed | 101 |

Every header, row count, row order and non-quality cell matches the preserved
file. The six filter-count and run-settings files are byte-identical. The final
variant membership is unchanged: controls 0 rows, CJD 4 rows and dilutions 0
rows. The PRNP/TET2/TTN gene tables retain 5/0/0, 22/1/2 and 12/0/1 rows for
controls, CJD and dilutions respectively.

## Quality gate

The predicate is `MEAN_BQ >= 20 AND MEAN_MQ >= 20`. Both operands use the same
threshold, so exchanging their labels cannot change the conjunction. Direct
comparison confirmed zero gate changes across all 43 summary rows: 5 controls,
25 CJD and 13 dilution rows. Filter counts and variant membership therefore
remain exactly as in the original outputs.

## Table 6 and combined tables PDF

Table 6 retains its four CJD rows. Only the displayed base- and mapping-quality
values changed:

| Sample | Position | Legacy BQ/MQ | Corrected BQ/MQ |
| --- | ---: | --- | --- |
| CJD23 | 4694249 | 30.01 / 20.62 | 20.62 / 30.01 |
| CJD23 | 4691920 | 38.49 / 21.00 | 21.00 / 38.49 |
| CJD2 | 4691920 | 37.90 / 22.01 | 22.01 / 37.90 |
| CJD6 | 4691920 | 33.44 / 20.10 | 20.10 / 33.44 |

The combined tables PDF was rebuilt twice with portable MiKTeX. It remains an
11-page document. Poppler text extraction and 120 dpi rendering both show that
only page 6, which contains Table 6, differs from the preserved PDF.

## Provenance

The script inventory retains its original 66-row scope and changes only the row
for `src/pipelines/4_readcount_to_tsv.py`. Sixteen declared final-output checksum
entries were recomputed; nine hashes changed and seven remained identical. No
path or description in `final_outputs_manifest.tsv` changed.

All 16 correction-set checksums verify. The complete repository checksum check
also ran and reported four pre-existing stale hashes for unrelated figure-source
TeX files. Those entries were not refreshed because figures are outside this
correction's dependency boundary.

The fields after `avg_se_mapping_quality` retain their legacy header names. This
correction deliberately changes only the two quality fields in scope; downstream
consumers select only `MEAN_BQ` and `MEAN_MQ` from the metrics tables.
