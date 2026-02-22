# Pipelines

This folder contains executable pipeline entry points.

## Available scripts

- `preflight_preprocessing.sh`: checks inputs/tools for preprocessing
- `preprocessing.sh`: FASTQ to final BAM preprocessing
- `mutect2_controls_no_pon.sh`: controls-only Mutect2 tumour-only calling (Stage 1)
- `mutect2_controls_postprocess_no_pon.sh`: controls-only post-processing (Stages 2-6)

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

### Stage 1

- `src/pipelines/mutect2_controls_no_pon.sh`

Writes:

- `runs/mutect2_controls_no_pon/vcf/`
- `runs/mutect2_controls_no_pon/f1r2/`

### Stages 2-6

- `src/pipelines/mutect2_controls_postprocess_no_pon.sh`

Writes:

- `runs/mutect2_controls_no_pon/orientation/`
- `runs/mutect2_controls_no_pon/filtered/`
- `runs/mutect2_controls_no_pon/scores/`
- `runs/mutect2_controls_no_pon/norm/`
- `runs/mutect2_controls_no_pon/annot/`

### Config keys used

Set in `config/preprocessing.env` (or via environment variables):

- `MUTECT2_CONTROLS_OUT_ROOT`
- `REF_FASTA`
- `FUNCOTATOR_DS`
- `JAVA_MEM_GB`
- `DRY_RUN`
