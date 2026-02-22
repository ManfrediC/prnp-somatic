# Pipelines

This folder contains executable pipeline entry points.

## Available scripts

- `preflight_preprocessing.sh`: checks inputs/tools for preprocessing
- `preprocessing.sh`: FASTQ to final BAM preprocessing
- `1_controls_mutect2_no_pon.sh`: controls-only Mutect2 tumour-only calling (Stage 1)
- `2_controls_postprocess_no_pon.sh`: controls-only post-processing (Stages 2-7)
- `3_controls_readcount_qc_no_pon.sh`: controls-only read/base quality metrics from annotated variants

## Preprocessing

### Inputs

- `config/preprocessing_samples.tsv` (committed)
- `config/preprocessing.env` (gitignored)
- `fastq/<batch_id>/...` (not committed)
- `resources/` references (committed)

### Outputs

- `runs/preprocessing/<batch_id>/<sample_id>/` (intermediates/logs)
- `results/final_bam/<sample_id>.bam` (+ index)

### Run

1. `conda activate prnp-somatic`
2. `src/pipelines/preflight_preprocessing.sh`
3. `DRY_RUN=1 src/pipelines/preprocessing.sh`
4. `DRY_RUN=0 src/pipelines/preprocessing.sh`

## Controls Variant Calling (No PoN)

Use this when Stage 1 Mutect2 has been run for controls and you want orientation/filtering/annotation.

Post-processing stage order is linear:

- Stage 2 `f1r2 -> orientation`
- Stage 3 `raw vcf + orientation + stats -> filtered (PASS)`
- Stage 4 `filtered -> scores`
- Stage 5 `scores -> norm`
- Stage 6 `norm -> annot`
- Stage 7 `annot -> annot_with_gnomad`

The script fails fast if required per-sample inputs are missing at any stage.

### Stage 1

- `src/pipelines/1_controls_mutect2_no_pon.sh`

Writes:

- `runs/mutect2_controls_no_pon/vcf/`
- `runs/mutect2_controls_no_pon/f1r2/`

### Stages 2-7

- `src/pipelines/2_controls_postprocess_no_pon.sh`

Writes:

- `runs/mutect2_controls_no_pon/orientation/`
- `runs/mutect2_controls_no_pon/filtered/`
- `runs/mutect2_controls_no_pon/scores/`
- `runs/mutect2_controls_no_pon/norm/`
- `runs/mutect2_controls_no_pon/annot/`
- `runs/mutect2_controls_no_pon/annot_with_gnomad/`

### Readcount QC

- `src/pipelines/3_controls_readcount_qc_no_pon.sh`

Writes:

- `runs/mutect2_controls_no_pon/readcount_qc/beds/`
- `runs/mutect2_controls_no_pon/readcount_qc/bam_work/`
- `runs/mutect2_controls_no_pon/readcount_qc/bam_nodup/`
- `runs/mutect2_controls_no_pon/readcount_qc/readcounts/`
- `runs/mutect2_controls_no_pon/readcount_qc/metrics/`

### Config keys used

Set in `config/preprocessing.env` (or via environment variables):

- `MUTECT2_CONTROLS_OUT_ROOT`
- `REF_FASTA`
- `FUNCOTATOR_DS`
- `GNOMAD_AF_VCF`
- `JAVA_MEM_GB`
- `DRY_RUN`

### Funcotator files required

`FUNCOTATOR_DS` must point to the datasource directory used by `gatk Funcotator --data-sources-path`.

Recommended value:

- `resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38`

Expected contents under that directory include source folders such as:

- `clinvar/`
- `dbsnp/`
- `gencode/`
- `hgnc/`

These are separate from the reference FASTA (`REF_FASTA`).

Note: `gnomAD_exome` and `gnomAD_genome` are intentionally excluded from the active
datasource tree to avoid requester-pays bucket access during Funcotator. Archived copies are in
`resources/backup/funcotator_excluded_datasources/`.

### gnomAD AF resource required

`GNOMAD_AF_VCF` is used in Stage 7 (`bcftools annotate`) to copy population AF into
`INFO/GNOMAD_AF`.

Recommended value:

- `resources/somatic-hg38_af-only-gnomad.hg38.vcf.gz`

The corresponding index (`.tbi` or `.csi`) must be present.
