# Runs Directory

This directory stores pipeline run outputs (intermediate and run-level artefacts).

Current top-level contents:

- `preprocessing/`
- `mutect2_controls_no_pon/`
- `mutect2_cjd_dilutions_with_pon/` (created by `src/pipelines/8_cjd_dilutions_mutect2_with_pon.sh`)

## `preprocessing/`

Per-sample outputs from `src/pipelines/preprocessing.sh`, organised as:

`runs/preprocessing/<batch_id>/<sample_id>/`

Current batches and sample counts:

- `CJD_16_samples` (`16` samples)
- `CJD_8_samples` (`8` samples)
- `first_CJD_seq` (`8` samples)
- `sequencing_of_dilutions` (`7` samples)

Current total sample directories: `39`.

Each sample directory currently contains:

- Trimmed and unpaired FASTQ files:
- `<sample>.R1.trimmed.fastq.gz`
- `<sample>.R2.trimmed.fastq.gz`
- `<sample>.R1.unpaired.fastq.gz`
- `<sample>.R2.unpaired.fastq.gz`
- Alignment and preprocessing BAMs:
- `<sample>.bwa.sorted.bam` and `<sample>.bwa.sorted.bam.bai`
- `<sample>.bwa.picard.bam`
- `<sample>.bwa.picard.markedDup.bam`
- `<sample>.bwa.picard.markedDup.recal.bam` and `<sample>.bwa.picard.markedDup.recal.bam.bai`
- Additional run artefacts:
- `<sample>.bwa.picard.markedDup.metrics`
- `<sample>.bwa.recal_data.table`
- `RUN_META.txt`
- `logs/`

## `mutect2_controls_no_pon/`

Outputs from control-only Mutect2 tumour-only runs without a panel of normals:

- `vcf/`
- `f1r2/`

Current controls present:

- `Ctrl1`
- `Ctrl2`
- `Ctrl3`
- `Ctrl4`
- `Ctrl5`
- `Ctrl7`

Current file patterns:

- `f1r2/<sample>.f1r2.tar.gz`
- `vcf/<sample>.raw.vcf.gz`
- `vcf/<sample>.raw.vcf.gz.tbi`
- `vcf/<sample>.raw.vcf.gz.stats`

## `mutect2_cjd_dilutions_with_pon/`

Stage-8 run directory for PoN-enabled Mutect2 over CJD and dilution BAMs.

Generated structure:

- `cjd/vcf/<sample>.raw.vcf.gz`
- `cjd/f1r2/<sample>.f1r2.tar.gz`
- `dilutions/vcf/<sample>.raw.vcf.gz`
- `dilutions/f1r2/<sample>.f1r2.tar.gz`
- `logs/stage8_mutect2_with_pon.log`

Current outputs:

- CJD: `26` raw VCFs and `26` F1R2 tar files
- Dilutions: `7` raw VCFs and `7` F1R2 tar files

This directory stores run/intermediate artefacts. Downstream result tables should be written under `results/`.
