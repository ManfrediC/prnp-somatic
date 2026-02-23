# Config

Committed:
- `preprocessing_samples.tsv` - sample manifest for preprocessing (`batch_id`, `sample_id`, `r1`, `r2`)
- `preprocessing.env.example` - template for local configuration

Current note:
- `preprocessing_samples.tsv` includes all `first_CJD_seq` CJD samples used in the authoritative manifest (`CJD1`, `CJD2`, `CJD6`, `CJD13`, `CJD22`, `CJD23`, `CJD25`, `CJD27`).

Not committed (local machine only):
- `preprocessing.env` - local paths and settings (ignored by `.gitignore`)

Rationale: `preprocessing.env` may include machine-specific paths and runtime settings; large data directories (`fastq/`, `runs/`, `results/final_bam/`) are intentionally not committed.

## Controls post-processing settings

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

The corresponding index (`.tbi` or `.csi`) must also be present.

For controls post-processing, `gnomAD_exome` and `gnomAD_genome` should remain excluded from
the active datasource tree to avoid requester-pays access requirements.
Archived copies are kept under:

- `resources/backup/funcotator_excluded_datasources/`

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

Note:

- Stage 12 always disables the final AAF filter for the `dilutions` group.
- `ENABLE_AAF_FILTER` and `AAF_THRESHOLD` remain configurable for the `cjd` group.
