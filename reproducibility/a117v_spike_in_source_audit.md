# A117V Spike-In Source Audit

Checked on: 2026-05-21
Scope: Checkpoint 1 of `plans/a117v_spike_in_remediation_plan.md`

## Purpose

Freeze the source-file evidence for relabelling the undiluted A117V spike-in
samples before any registry, script, output-table, or manuscript edits are made.
This note separates observed evidence from interpretation and preserves the
original sequencing labels as provenance aliases.

Original external source files were read only. They were not modified.

## Conclusion

The sample submitted and sequenced as `NA995A05_undil_G01` should not be used as
a nominal 0.5% A117V LoD anchor. The preparation workbook row carrying the 0.5%
label used `0.05` as the A117V heterozygous-source fraction. Given the recorded
formula-derived inputs, this corresponds to approximately 199 ng NA source DNA
plus 10 ng heterozygous A117V-source DNA, or an expected mutant allele fraction
of about 2.39% if the A117V source is heterozygous and all other inputs are wild
type at A117V. In round terms, this is the approximately 2.4-2.5% VAF expected
from a 5% heterozygous-source spike-in.

The final display labels to use downstream are:

| Final display label | Raw sample alias | Sequencing stem | Interpretation |
|---|---|---|---|
| `A117V_low` | `NA99A1_undil` | `20210203.0-o23928_1_2-NA99A1_undil_D02` | Provenance-valid intended 1% A117V spike-in; observed 35/4315 reads = 0.81% VAF in the current with-PoN summary. |
| `A117V_high` | `NA995A05_undil` | `20210203.0-o23928_1_7-NA995A05_undil_G01` | Higher-input A117V spike-in by preparation calculation; 5% heterozygous-source input, expected about 2.4-2.5% VAF; observed 92/6961 reads = 1.32% VAF in the no-PoN manuscript-table source. |

## Original Source Evidence

### Preparation Calculation Workbook

Path:
`C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect\2021-01-22 prepare final mosaics\2021-01-22 mosaic calculation.xlsx`

Relevant sheet: `Sheet1`

Relevant rows:

| Row | Label or field | Evidence |
|---|---|---|
| 2 | Headers | `fraction of NA`; `fraction of A117V (hetero)`; `NA ul`; `A117V ul`; `low TE` |
| 5 | `NA (new) 99.5% A 0.5% undil probes` | `fraction of NA = 0.995`; `fraction of A117V (hetero) = 0.05`; `NA ul = 34.310344827586206`; `A117V ul = 0.28901734104046239`; `low TE = 15.400637831373331` |
| 12 | A117V-source concentration | `A117V 1:20 = 34.6 ng/ul` |
| 13 | NA-source concentration | `NA new batch = 5.8 ng/ul` |

Formula evidence from the same row:

| Cell | Formula | Cached value | Interpretation |
|---|---|---:|---|
| `D5` | `=(200/$B$13)*B5` | 34.310344827586206 | Approximately 199 ng NA-source DNA. |
| `E5` | `=(200/$B$12)*C5` | 0.28901734104046239 | Approximately 10 ng heterozygous A117V-source DNA. |

Calculation check:

- NA source: `34.310344827586206 ul * 5.8 ng/ul = 199 ng`
- A117V heterozygous source: `0.28901734104046239 ul * 34.6 ng/ul = 10 ng`
- Expected mutant allele fraction from these inputs:
  `10 ng * 0.5 / (199 ng + 10 ng) = 2.39%`

Interpretation: the workbook row label says 0.5%, but the A117V
heterozygous-source fraction is `0.05`, corresponding to an expected mutant VAF
of approximately 2.4-2.5% before empirical spike-in and sequencing effects.

### Submission Workbook

Path:
`C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect\2021-01-25 test mosaic submission\2021-01-25 test mosaic sample dilutions.xlsx`

Relevant sheet: `Sheet1`

| Row | Submitted label | Submitted sample name | Well | Molarity-like field |
|---:|---|---|---|---:|
| 3 | `NA 99% A 1% undil probes` | `NA99A1_undil_D02` | `D02` | 22.8 |
| 4 | `NA 99% A 1% 1:5 probes` | `NA99A1_1to5_E01` | `E01` | 6.44 |
| 8 | `NA 99.5% A 0.5% undil probes` | `NA995A05_undil_G01` | `G01` | 13.4 |

Interpretation: the submission labels preserve the intended labels passed
forward for sequencing, but they do not correct the upstream preparation
calculation inconsistency for `NA995A05_undil_G01`.

### Test Mosaic List

Paths checked:

- `C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect\List_of_test_mosaics.xlsx`
- `C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect\List_of_test_mosaics-WS123363.xlsx`

Relevant rows in both copies:

| Label | Well | Mutation | Listed mutant percentage | c129-M | c129-V | Probe dilution | Molarity-like field |
|---|---|---|---:|---:|---:|---|---:|
| `NA 99% A 1% undil probes` | `D02` | `A117V` | 1 | 99 | 1 | none | 22.8 |
| `NA 99% A 1% 1:5 probes` | `E01` | `A117V` | 1 | 99 | 1 | 1:5 | 6.44 |
| `NA 99.5% A 0.5% undil probes` | `G01` | `A117V` | 0.5 | 99.5 | 0.5 | none | 13.4 |

Interpretation: the list workbook documents the submitted intended labels. It
does not supersede the separate preparation calculation row showing `0.05` as
the A117V heterozygous-source fraction for the sample submitted as
`NA995A05_undil_G01`.

### Targeted Workbook Scan

A read-only scan of the original SureSelect `.xlsx` files found the low-sample
labels in the submission and test-mosaic-list workbooks. The explicit
preparation-calculation row found for the high sample was present in these
calculation workbook copies:

- `2021-01-22 mosaic calculation.xlsx`
- `2021-01-22 mosaic calculation-WS123363.xlsx`
- `2021-01-22 mosaic calculation-WS123363-2.xlsx`

The targeted scan did not find an equivalent `NA 99% A 1% undil probes`
preparation-calculation row in the checked 22 Jan calculation workbooks.
Therefore, `A117V_low` is treated as provenance-valid based on the submission
and test-mosaic-list labels plus the absence of contradictory source evidence in
this scan, not because a matching low-sample preparation-calculation row was
identified.

## Sequencing Metadata Evidence

Original dataset metadata path:
`C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Analysis\Sureselect_pipeline\pipeline\dataset-1.tsv`

Relevant rows:

| Dataset name | Read 1 | Read 2 | B-Fabric sample ID | Order ID | Read count |
|---|---|---|---:|---:|---:|
| `o23928_1_2-NA99A1_undil_D02` | `p3111/MiSeq_210203_MS567_o23928/20210203.0-o23928_1_2-NA99A1_undil_D02_R1.fastq.gz` | `p3111/MiSeq_210203_MS567_o23928/20210203.0-o23928_1_2-NA99A1_undil_D02_R2.fastq.gz` | 282860 | 23928 | 2452630 |
| `o23928_1_3-NA99A1_1to5_E01` | `p3111/MiSeq_210203_MS567_o23928/20210203.0-o23928_1_3-NA99A1_1to5_E01_R1.fastq.gz` | `p3111/MiSeq_210203_MS567_o23928/20210203.0-o23928_1_3-NA99A1_1to5_E01_R2.fastq.gz` | 282862 | 23928 | 3318867 |
| `o23928_1_7-NA995A05_undil_G01` | `p3111/MiSeq_210203_MS567_o23928/20210203.0-o23928_1_7-NA995A05_undil_G01_R1.fastq.gz` | `p3111/MiSeq_210203_MS567_o23928/20210203.0-o23928_1_7-NA995A05_undil_G01_R2.fastq.gz` | 282865 | 23928 | 4247919 |

Interpretation: the sequencing metadata consistently carries the submitted
sample names through to the FASTQ stems. This supports preserving the raw names
as aliases, rather than rewriting raw provenance labels.

Current repo registry rows checked:

- `config/preprocessing_samples.tsv`
- `config/sample_manifest.tsv`

These files currently use the raw aliases `NA99A1_undil`, `NA99A1_1to5`, and
`NA995A05_undil`. They have not yet been remediated.

## Observed A117V Evidence

The observed A117V signal is higher in the sample submitted as
`NA995A05_undil_G01` than in the provenance-valid intended 1% sample. This
pattern is reproducible across legacy outputs and is consistent with the source
preparation calculation error.

| Source output | Sample label in output | REF | ALT | Read depth | REF count | ALT count | VAF / AAF |
|---|---|---|---|---:|---:|---:|---:|
| `legacy/src/windows/2025-03-15_dilution_force_call_summaries/original/mutation_csv/A117V_summary_table.csv` | `20210203.0-o23928_1_2-NA99A1_undil_D02` | `CA` | `TG,TA` | 3129 | 3105 | 24 | 0.767% |
| same | `20210203.0-o23928_1_3-NA99A1_1to5_E01` | `CA` | `TG,TA` | 1546 | 1540 | 6 | 0.388% |
| same | `20210203.0-o23928_1_7-NA995A05_undil_G01` | `CA` | `TG,TA` | 4622 | 4562 | 60 | 1.298% |
| `legacy/src/windows/2025-03-15_dilution_force_call_summaries/bwa/mutation_csv/A117V_summary_table.csv` | `NA99A1_undil` | `CA` | `TG,TA` | 4237 | 4201 | 36 | 0.850% |
| same | `NA99A1_1to5` | `CA` | `TG,TA` | 2332 | 2317 | 15 | 0.643% |
| same | `NA995A05_undil` | `CA` | `TG,TA` | 6899 | 6800 | 99 | 1.435% |
| `legacy/src/windows/2025-05-03_dilution_QC_analysis/A117V_table.csv` | `A117V 1% spike-in` | `C` | `T` | 4315 | 4280 | 35 | 0.811% |
| same | `A117V 0.5% spike-in` | `C` | `T` | 6961 | 6869 | 92 | 1.322% |
| `results/mutect2_cjd_dilutions_with_pon/variant_qc/dilutions/summary_combined_variants.tsv` | `NA99A1_undil` | `CA` | `TG` | 4315 | 4280 | 35 | 0.811% |

Interpretation:

- `NA99A1_undil` should become `A117V_low`.
- `NA995A05_undil` should become `A117V_high`.
- The read counts and VAF/AAF values are evidence to preserve, not values to
  change during label remediation.
- The current with-PoN summary contains the provenance-valid low A117V row. The
  no-PoN manuscript-table source contains both historical display labels and is
  the source of the 92/6961 high-sample value currently used in manuscript
  tables.

## Evidence Versus Interpretation

Evidence:

- The preparation calculation workbook row labelled
  `NA (new) 99.5% A 0.5% undil probes` has `fraction of A117V (hetero) = 0.05`.
- Submission and test-mosaic-list workbooks label `NA99A1_undil_D02` as a 1%
  A117V undiluted spike-in and `NA995A05_undil_G01` as a 0.5% A117V undiluted
  spike-in.
- Sequencing metadata preserves the submitted names and wells through FASTQ
  stems and B-Fabric sample IDs.
- Legacy and current outputs show `NA995A05_undil`/historical `A117V 0.5%
  spike-in` with higher observed A117V signal than
  `NA99A1_undil`/historical `A117V 1% spike-in`.

Interpretation:

- The most defensible remediation is not to claim a sequencing/sample-name swap.
  The stronger evidence points to a preparation-calculation inconsistency for
  the sample submitted as `NA995A05_undil_G01`.
- `NA995A05_undil` should be preserved as a raw provenance alias but removed as
  a final scientific display label.
- The final display label should be `A117V_high`, with a note that it had 5%
  heterozygous-source input, expected about 2.4-2.5% VAF from the calculation,
  and observed 1.32% VAF.
- `NA99A1_undil` should be preserved as a raw provenance alias and displayed as
  `A117V_low`, with a note that it is the provenance-valid intended 1% sample
  and had observed 0.81% VAF.

## Privacy And Scope Check

This note uses only project sample IDs, workbook row labels, sequencing stems,
and aggregate read-count/VAF evidence. It does not add patient-identifying
details beyond the sample IDs and provenance labels already present in the
project repository and original project filesystem.
