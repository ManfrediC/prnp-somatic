# prnp-somatic

PRNP somatic variant analysis:
FASTQ -> BAM (GATK best practices) -> Mutect2 -> QC -> result tables/figures.

## Reproduction Guide

Run all commands from repository root (`/add/prnp-somatic`).

### 1. Clone and enter repo

```bash
git clone https://github.com/ManfrediC/prnp-somatic.git
cd prnp-somatic
```

### 2. Obtain required resources and package context

Before running workflows, use:

- `doc/reproducibility/tooling_and_reference_provenance.md`

This document lists:

- required reference resources and their source/provenance
- required CLI/R tools and version context
- environment-specific package notes (including junction and ddPCR dependencies)

### 3. Prepare input data in expected locations

Required inputs must be placed as follows:

- ddPCR raw files: `ddPCR/*.csv`
- ddPCR sample metadata: `ddPCR/sample_details.xlsx`
- raw FASTQ files: `fastq`
- resources for GATK workflow: `resources`

See `README` files in the respective directories for details on the required files.

### 4. Create environment for ddPCR + somatic SNV pipeline scripts

```bash
conda --no-plugins create -n prnp-somatic -c conda-forge -c bioconda -y \
  r-base=4.4 r-tidyverse r-openxlsx r-readr r-magrittr r-binom \
  r-dplyr r-stringr r-tidyr r-tibble
conda activate prnp-somatic
```

### 5. Run ddPCR workflow

```bash
bash src/ddPCR/run_ddpcr.sh
```

Expected outputs:

- `results/ddPCR/SNV_data_final.xlsx`
- `results/ddPCR/SNV_pooled_participant.xlsx`
- `results/ddPCR/p0_fallback.csv`

### 6. Run SNV detection pipeline (`src/pipelines`)

Set up pipeline configuration first:

```bash
cp config/preprocessing.env.example config/preprocessing.env
```

Then edit `config/preprocessing.env` for local paths/settings as documented in `config/README.md`.

Run full pipeline in order:

```bash
conda activate prnp-somatic
bash src/pipelines/preflight_preprocessing.sh
DRY_RUN=1 bash src/pipelines/preprocessing.sh
DRY_RUN=0 bash src/pipelines/preprocessing.sh
bash src/pipelines/1_controls_mutect2_no_pon.sh
bash src/pipelines/2_controls_postprocess_no_pon.sh
bash src/pipelines/3_controls_readcount_qc_no_pon.sh
bash src/pipelines/5_controls_variant_qc_no_pon.sh
bash src/pipelines/7_controls_create_pon.sh
bash src/pipelines/8_cjd_dilutions_mutect2_with_pon.sh
bash src/pipelines/9_cjd_dilutions_postprocess_with_pon.sh
bash src/pipelines/10_cjd_dilutions_readcount_qc_with_pon.sh
bash src/pipelines/11_cjd_dilutions_readcount_to_tsv_with_pon.sh
bash src/pipelines/12_cjd_dilutions_variant_qc_with_pon.sh
```

Expected outputs include:

- `results/final_bam/*.bam` (+ indexes)
- `results/mutect2_controls_no_pon/variant_qc/*`
- `results/mutect2_controls_pon/panel_of_normals/CJD_controls_PoN.vcf.gz`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/cjd/*`
- `results/mutect2_cjd_dilutions_with_pon/variant_qc/dilutions/*`

If preprocessing/Mutect outputs already exist and you only need final Stage-12 tables:

```bash
conda activate prnp-somatic
bash src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh
```

### 7. Create environment for junction workflow

```bash
conda --no-plugins create -n prnp-junctions -c conda-forge -c bioconda -y \
  r-base=4.4 bioconductor-genomicfeatures bioconductor-txdbmaker \
  bioconductor-genomicranges bioconductor-genomeinfodb bioconductor-rsamtools \
  bioconductor-biostrings bioconductor-rtracklayer bioconductor-genomicalignments \
  r-dplyr
conda activate prnp-junctions
```

### 8. Run junction workflow

```bash
TMPDIR=/tmp TEMP=/tmp TMP=/tmp bash src/junctions/run_junctions.sh
```

Expected outputs:

- `results/junctions/junction_counts/prnp_junction_counts.tsv`
- `results/junctions/junction_counts/prnp_junction_summary.tsv`

### 9. Optional: regenerate manuscript artifacts

```bash
conda activate prnp-somatic
Rscript manuscript/run_all.R
```

### 10. Verify final outputs

```bash
bash bin/verify_output_checksums.sh --mode check
```

Expected manifest/checksum files:

- `doc/reproducibility/final_outputs_manifest.tsv`
- `doc/reproducibility/final_outputs.sha256`

Detailed workflow docs:

- `src/ddPCR/README.md`
- `src/junctions/README.md`
- `src/pipelines/README.md`
- `manuscript/README.md`

Makefile helper targets are documented in:

- `doc/MAKEFILE.md`

## Controls Variant Calling (No PoN)

For controls-only somatic variant calling without a panel of normals:

1. Run Stage 1 Mutect2: `src/pipelines/1_controls_mutect2_no_pon.sh`
2. Run Stages 2-7 post-processing: `src/pipelines/2_controls_postprocess_no_pon.sh`
3. Run readcount metrics: `src/pipelines/3_controls_readcount_qc_no_pon.sh`
4. Run variant extraction + QC: `src/pipelines/5_controls_variant_qc_no_pon.sh` (calls `src/pipelines/6_controls_variant_table_qc_no_pon.R`)
5. Final controls QC tables are written to: `results/mutect2_controls_no_pon/variant_qc/`
6. Create controls PoN: `src/pipelines/7_controls_create_pon.sh` (writes to `results/mutect2_controls_pon/panel_of_normals/`)
7. Run CJD + dilution Mutect2 with PoN: `src/pipelines/8_cjd_dilutions_mutect2_with_pon.sh`
8. Run CJD + dilution post-processing with PoN: `src/pipelines/9_cjd_dilutions_postprocess_with_pon.sh`
9. Run CJD + dilution readcount collection: `src/pipelines/10_cjd_dilutions_readcount_qc_with_pon.sh`
10. Parse CJD + dilution readcounts to metrics TSV: `src/pipelines/11_cjd_dilutions_readcount_to_tsv_with_pon.sh`
11. Run CJD + dilution variant table + QC: `src/pipelines/12_cjd_dilutions_variant_qc_with_pon.sh` (calls `src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R`)

Pipeline details and outputs are documented in `src/pipelines/README.md`.

### Funcotator resources

The controls post-processing pipeline requires a local Funcotator datasource tree.

- Config keys in `config/preprocessing.env`:
- `REF_FASTA`
- `FUNCOTATOR_DS`
- `GNOMAD_AF_VCF`
- `MANUAL_POP_FREQ_TSV`
- Recommended repo-relative paths:
- `resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38`
- `resources/somatic-hg38_af-only-gnomad.hg38.vcf.gz`

Reference FASTA (`REF_FASTA`) and Funcotator datasources (`FUNCOTATOR_DS`) are separate requirements.
`gnomAD_exome` and `gnomAD_genome` are intentionally excluded from the active Funcotator datasource tree to avoid requester-pays access. They are archived under `resources/backup/funcotator_excluded_datasources/`.

## Repository structure
- `bin/` command-line entrypoints / wrappers
- `src/` reusable code (Python/R)
- `config/` configuration templates
- `resources/` small static artefacts tracked in git (BEDs, schemas). Larger resources must be downloaded manually
- `authoritative_files/` manifests and sequencing-metrics validation utilities
- `results/` outputs *not tracked* (date-stamped runs)
- `env/` container/environment definitions
- `doc/` notes and documentation

## Raw Data Placement (Ignored)

For the workflow families in current reproducibility scope:

- ddPCR raw exports: `ddPCR/*.csv` (ignored via `ddPCR/.gitignore`)
- ddPCR metadata sheet: `ddPCR/sample_details.xlsx` (ignored via `ddPCR/.gitignore`)
- junction BAM inputs: `results/final_bam/*.bam` (ignored via repo `.gitignore` patterns for BAM/results)
- pipeline run intermediates: `runs/**` (ignored via repo `.gitignore`)
- pipeline/junction outputs: `results/**` (ignored except tracked placeholders)

Generated runtime index sidecars are also ignored (for example `resources/junctions/*.fa.{amb,ann,bwt,pac,sa}`).

## Conda environment
- Analyses were performed using the `Conda` environment defined in `env/environment.yml`.

## Data availability
Raw data files can be obtained for academic purposes upon reasonable request to the repository owner.
