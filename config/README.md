# Config

Committed:
- `preprocessing_samples.tsv` - sample manifest for preprocessing (`batch_id`, `sample_id`, `R1`, `R2`)
- `preprocessing.env.example` - template for local configuration

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
