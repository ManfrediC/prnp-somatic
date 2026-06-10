# Config

Committed:
- `preprocessing_samples.tsv` - sample manifest for preprocessing (`batch_id`, `sample_id`, `r1`, `r2`, optional `display_label`)
- `preprocessing.env.example` - template for local configuration
- `junctions.env.example` - template for exon-exon junction workflow configuration
- `repeats.env.example` - template for PRNP ORR repeat workflow configuration

Not committed (local machine only):
- `preprocessing.env` - local paths and settings (ignored by `.gitignore`)
- `junctions.env` - local paths and settings for `src/junctions/*` scripts
- `repeats.env` - local paths and runtime settings for `src/repeats/*`

Rationale: `preprocessing.env` may include machine-specific paths and runtime settings; large data directories (`raw/fastq/`, `runs/`, `results/final_bam/`) are intentionally not committed.

## Controls post-processing: Settings and paths

For controls post-processing and annotation, set these in `preprocessing.env`:

- `FUNCOTATOR_DS`
- `REF_FASTA`
- `GNOMAD_AF_VCF`
- `VARIANT_QC_ROOT`
- `VARIANT_QC_VCF_DIR`
- `VARIANT_QC_READCOUNT_METRICS_DIR`
- `VARIANT_QC_OUTPUT_DIR`
- `VARIANT_QC_R_SCRIPT`
- `MANUAL_POP_FREQ_TSV`
- `ENABLE_AAF_FILTER`
- `AAF_THRESHOLD`

Recommended `FUNCOTATOR_DS` value:

- `resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38`

Recommended `GNOMAD_AF_VCF` value:

- `resources/somatic-hg38_af-only-gnomad.hg38.vcf.gz`

Recommended `MANUAL_POP_FREQ_TSV` value:

- `resources/annotations/manual_population_freq.tsv`

Recommended `VARIANT_QC_OUTPUT_DIR` value:

- `results/mutect2_controls_no_pon/variant_qc`

Recommended `VARIANT_QC_R_SCRIPT` value:

- `src/pipelines/6_controls_variant_table_qc_no_pon.R`

## Repeat workflow settings

For `src/repeats/01_run_prnp_orr.sh`, these are the main override variables:

- `REPEAT_BAM_DIR`
- `REPEAT_RESULTS_ROOT`
- `REPEAT_REF_FASTA`
- `PRNP_ORR_BED`
- `PRNP_EH_CATALOG`
- `RUN_GANGSTR`
- `GANGSTR_BIN`
- `GANGSTR_REGIONS_BED`
- `GANGSTR_READLENGTH`
- `GANGSTR_INSERTMEAN`
- `GANGSTR_INSERTSDEV`
- `GANGSTR_MAX_PROC_READ`
- `PRNP_TOTAL_REFERENCE_REPEATS`
- `PRNP_VARIABLE_REPEAT_OFFSET`
- `REPEAT_THREADS`
- `ORR_MIN_READS`
- `FORCE`
- `ARCHIVE_EXISTING_RUN`
- `ARCHIVE_RUN_LABEL`
- `EXPECTED_SAMPLE_COUNT`

Recommended values:

- `REPEAT_BAM_DIR="results/final_bam"`
- `REPEAT_RESULTS_ROOT="results/repeats"`
- `REPEAT_REF_FASTA="resources/chr2_chr4_chr20.fasta"`
- `PRNP_ORR_BED="resources/prnp_orr.hg38.bed"`
- `PRNP_EH_CATALOG="resources/repeats/prnp_orr.expansionhunter.json"`
- `RUN_GANGSTR=0`
- `GANGSTR_BIN=""`
- `GANGSTR_REGIONS_BED="resources/repeats/prnp_orr.gangstr.bed"`
- `GANGSTR_READLENGTH=""`
- `GANGSTR_INSERTMEAN=""`
- `GANGSTR_INSERTSDEV=""`
- `GANGSTR_MAX_PROC_READ=1000000`
- `PRNP_TOTAL_REFERENCE_REPEATS=5`
- `PRNP_VARIABLE_REPEAT_OFFSET=3`
- `REPEAT_THREADS=4`
- `ORR_MIN_READS=10`
- `FORCE=0`
- `ARCHIVE_EXISTING_RUN=0`
- `ARCHIVE_RUN_LABEL=""`
- `EXPECTED_SAMPLE_COUNT=32`

GangSTR notes:

- when `RUN_GANGSTR=1`, the workflow runs GangSTR in `--targeted --nonuniform` mode as an orthogonal local STR genotyper for the PRNP ORR mutable repeat block
- leave `GANGSTR_BIN=""` to use `GangSTR` from the active `prnp-repeats` environment, or set an explicit binary path if it is installed elsewhere
- leave `GANGSTR_READLENGTH`, `GANGSTR_INSERTMEAN`, and `GANGSTR_INSERTSDEV` empty to auto-estimate them from each BAM at runtime, or set explicit values if you need to override that behavior
- `GANGSTR_MAX_PROC_READ` is passed through to GangSTR to cap how many reads it evaluates per sample

Clean rerun behavior:

- if prior outputs are present in `results/repeats`, the workflow now stops by default instead of mixing runs
- set `ARCHIVE_EXISTING_RUN=1` to move the current live run outputs into `results/repeats/old_runs/<timestamp>/` before starting a fresh rerun
- optionally set `ARCHIVE_RUN_LABEL` to control the archive directory name
- `EXPECTED_SAMPLE_COUNT=32` enforces the current cohort size and helps catch partial manifests or missing BAMs before the long run proceeds

Manual mosaic review helpers:

- `src/repeats/06_manual_mosaic_prnp_orr.py`
- `src/repeats/07_run_manual_mosaic_prnp_orr_cohort.py`
- `src/repeats/08_filter_manual_mosaic_prnp_orr_cohort.py`

These do not read `config/repeats.env`; they take their thresholds and paths from
CLI arguments. By default they reuse:

- `results/final_bam`
- `resources/chr2_chr4_chr20.fasta`
- `resources/repeats/prnp_orr_manual_panel.tsv`
- `results/repeats/sample_calls.tsv`
- `results/repeats/gangstr_calls.tsv`

The corresponding index (`.tbi` or `.csi`) must also be present.

## Variant QC toggle

`src/pipelines/5_controls_variant_qc_no_pon.sh` supports temporary disabling of the final AAF
threshold filter:

- `ENABLE_AAF_FILTER=1` applies `AAF > AAF_THRESHOLD` (default, recommended)
- `ENABLE_AAF_FILTER=0` keeps all rows after the upstream QC filters for review

## Controls PoN creation settings

For `src/pipelines/7_controls_create_pon.sh`, set:

- `PON_INPUT_DIR`
- `PON_OUTPUT_ROOT`
- `PON_MERGED_VCF`
- `PON_VCF`
- `PON_CONTROLS`

Recommended values:

- `PON_INPUT_DIR="runs/mutect2_controls_no_pon/filtered"`
- `PON_OUTPUT_ROOT="results/mutect2_controls_pon/panel_of_normals"`
- `PON_MERGED_VCF="results/mutect2_controls_pon/panel_of_normals/controls_multisample.filtered.vcf.gz"`
- `PON_VCF="results/mutect2_controls_pon/panel_of_normals/CJD_controls_PoN.vcf.gz"`
- `PON_CONTROLS="Ctrl1 Ctrl2 Ctrl3 Ctrl4 Ctrl5 Ctrl7"`

## CJD + dilution Mutect2 with PoN settings

For `src/pipelines/8_cjd_dilutions_mutect2_with_pon.sh`, set:

- `MUTECT2_WITH_PON_OUT_ROOT`
- `PON_VCF`
- `DILUTION_SAMPLES`

Recommended values:

- `MUTECT2_WITH_PON_OUT_ROOT="runs/mutect2_cjd_dilutions_with_pon"`
- `PON_VCF="results/mutect2_controls_pon/panel_of_normals/CJD_controls_PoN.vcf.gz"`
- `DILUTION_SAMPLES="NA100_undil NA100_1to10 NA99A1_undil A100_1to2 NA99A1_1to5 NA995A05_undil NA100_1to2"`

Note: `DILUTION_SAMPLES` uses raw BAM/sample IDs for file lookup and stable
pipeline ordering. Manuscript-facing labels for the A117V spike-ins are
assigned through sample-registry `display_label` values
(`NA99A1_undil` -> `A117V low`; `NA995A05_undil` -> `A117V high`).

## CJD + dilution post-processing with PoN settings

For `src/pipelines/9_cjd_dilutions_postprocess_with_pon.sh`, set:

- `MUTECT2_WITH_PON_POSTPROCESS_ROOT`
- `WITH_PON_GROUPS`
- `FUNCOTATOR_DS`
- `REF_FASTA`
- `INTERVALS`
- `GNOMAD_AF_VCF`

Recommended values:

- `MUTECT2_WITH_PON_POSTPROCESS_ROOT="runs/mutect2_cjd_dilutions_with_pon"`
- `WITH_PON_GROUPS="cjd dilutions"`
- `FUNCOTATOR_DS="resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38"`
- `REF_FASTA="resources/chr2_chr4_chr20.fasta"`
- `INTERVALS="resources/capture_targets.interval_list"`
- `GNOMAD_AF_VCF="resources/somatic-hg38_af-only-gnomad.hg38.vcf.gz"`

## CJD + dilution readcount collection settings

For `src/pipelines/10_cjd_dilutions_readcount_qc_with_pon.sh`, set:

- `MUTECT2_WITH_PON_READCOUNT_ROOT`
- `WITH_PON_READCOUNT_GROUPS`
- `WITH_PON_READCOUNT_BAM_DIR`
- `WITH_PON_READCOUNT_REF_FASTA`

Recommended values:

- `MUTECT2_WITH_PON_READCOUNT_ROOT="runs/mutect2_cjd_dilutions_with_pon"`
- `WITH_PON_READCOUNT_GROUPS="cjd dilutions"`
- `WITH_PON_READCOUNT_BAM_DIR="results/final_bam"`
- `WITH_PON_READCOUNT_REF_FASTA="resources/chr2_chr4_chr20.fasta"`

## CJD + dilution readcount parsing settings

For `src/pipelines/11_cjd_dilutions_readcount_to_tsv_with_pon.sh`, set:

- `READCOUNT_TO_TSV_PY`

Recommended value:

- `READCOUNT_TO_TSV_PY="src/pipelines/4_readcount_to_tsv.py"`

## CJD + dilution variant-table + QC settings

For `src/pipelines/12_cjd_dilutions_variant_qc_with_pon.sh`, set:

- `MUTECT2_WITH_PON_VARIANT_QC_ROOT`
- `WITH_PON_VARIANT_QC_GROUPS`
- `WITH_PON_VARIANT_QC_RESULTS_ROOT`
- `WITH_PON_VARIANT_QC_R_SCRIPT`
- `MANUAL_POP_FREQ_TSV`

Recommended values:

- `MUTECT2_WITH_PON_VARIANT_QC_ROOT="runs/mutect2_cjd_dilutions_with_pon"`
- `WITH_PON_VARIANT_QC_GROUPS="cjd dilutions"`
- `WITH_PON_VARIANT_QC_RESULTS_ROOT="results/mutect2_cjd_dilutions_with_pon/variant_qc"`
- `WITH_PON_VARIANT_QC_R_SCRIPT="src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R"`
- `MANUAL_POP_FREQ_TSV="resources/annotations/manual_population_freq.tsv"`
