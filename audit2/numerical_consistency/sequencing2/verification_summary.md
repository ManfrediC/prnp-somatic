# Sequencing2 downstream verification

The isolated downstream artefacts use the final sequencing2 CJD PRNP table as their single call source.

- Table membership: three sample-level calls at two loci.
- Figure membership: identical to the table input.
- Figure design: original reviewed SVG layout, R head paths, axes, vertical
  labels, branching stems, sample legend and exon track retained; CJD6 removed.
- LoD: exact internal value `26/3891 = 0.006682086867129272`; displayed as `0.668%`.
- Controls: zero final calls.
- Dilutions: AAF filtering is recorded as disabled.
- Quality mapping: `MEAN_BQ` is base quality and `MEAN_MQ` is mapping quality.
- Canonical manuscript: unchanged.

Machine checks and visual PDF review are recorded in the task handoff and enforced by `tests2/regression/test_sequencing2_variant_table_provenance.R`.
