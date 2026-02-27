# prnp-somatic

This repository contains the reproducible analysis workflows used in the *PRNP* somatic project:

- NGS somatic SNV detection (`src/pipelines`): FASTQ preprocessing, Mutect2 calling, post-processing, QC, and final variant tables.
- ddPCR SNV quantification (`src/ddPCR`): processing of raw droplet exports into long-format and participant-level result tables.
- PRNP exon-exon junction analysis (`src/junctions`): junction-reference construction, read realignment, and junction-count quantification.
- Manuscript artifact generation (`manuscript`): scripts that build figures/tables from outputs in `results/`.

The Reproduction Guide below details required inputs, environment setup and expected outputs.

## Repository structure
- `bin/` command-line entrypoints / wrappers
- `src/` reusable code (Shell/Python/R)
- `config/` configuration templates
- `resources/` small static artefacts tracked in git (BEDs, schemas). Larger resources must be downloaded manually
- `authoritative_files/` manifests and sequencing-metrics validation utilities
- `results/` outputs *not tracked* (date-stamped runs)
- `env/` container/environment definitions
- `doc/` notes and documentation
- `manuscript/` manuscript artefact generation

## Root files
- `README.md`: top-level project overview and run guide
- `Makefile`: convenience wrappers for running workflows and checks
- `CITATION.cff`: citation metadata for this repository
- `LICENCE`: repository software licence (MIT)
- `.gitignore`: rules for untracked raw data, run outputs and generated artefacts
- `.gitattributes`: git text/binary handling rules for selected tracked files

## Reproduction Guide

Run all commands from the repository root (your cloned `prnp-somatic` directory).

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

Funcotator resource acquisition, layout and datasource details are documented in:

- `resources/README.md`

### 3. Prepare input data in expected locations

Required inputs must be placed as follows:

- ddPCR raw files: `ddPCR/*.csv`
- ddPCR sample metadata: `ddPCR/sample_details.xlsx`
- raw FASTQ files: `fastq`
- resources for GATK workflow: `resources`

See `README` files in the respective directories for details on the required files.

### 4. Create environment for ddPCR scripts

```bash
conda env create -f env/ddpcr.environment.yml
conda activate prnp-somatic-ddpcr
```

### 5. Create environment for somatic SNV pipeline scripts

```bash
conda --no-plugins create -n prnp-somatic -c conda-forge -c bioconda -y \
  python=3.10 pip samtools=1.20 bcftools=1.20 htslib=1.20 \
  bedtools=2.31.1 fastqc=0.12.1 multiqc=1.33 openjdk=17 gatk4=4.5.0.0
conda activate prnp-somatic
```

### 6. Run ddPCR workflow

```bash
bash src/ddPCR/run_ddpcr.sh
```

Expected outputs:

- `results/ddPCR/SNV_data_final.xlsx`
- `results/ddPCR/SNV_pooled_participant.xlsx`

### 7. Run SNV detection pipeline (`src/pipelines`)

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

### 8. Create environment for junction workflow

```bash
conda --no-plugins create -n prnp-junctions -c conda-forge -c bioconda -y \
  r-base=4.4 bioconductor-genomicfeatures bioconductor-txdbmaker \
  bioconductor-genomicranges bioconductor-genomeinfodb bioconductor-rsamtools \
  bioconductor-biostrings bioconductor-rtracklayer bioconductor-genomicalignments \
  r-dplyr
conda activate prnp-junctions
```

### 9. Run junction workflow

```bash
TMPDIR=/tmp TEMP=/tmp TMP=/tmp bash src/junctions/run_junctions.sh
```

Expected outputs:

- `results/junctions/junction_counts/prnp_junction_counts.tsv`
- `results/junctions/junction_counts/prnp_junction_summary.tsv`

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

### 11. Makefile commands (optional wrappers)

As an alternative to executing the pipelines individually, it's possible to run each part of the project using `Makefile`.

Show all available commands:

```bash
make help
```

Run main workflows:

- `make ddpcr` (requires active env: `prnp-somatic-ddpcr`)
- `make snv` (requires active env: `prnp-somatic`)
- `make junctions` (requires active env: `prnp-junctions`)
- `make all` (runs ddPCR + SNV + junctions via `conda run`; expects envs `prnp-somatic-ddpcr`, `prnp-somatic`, and `prnp-junctions`)

Run integrity checks:

- `make check` (requires active env: `prnp-somatic`; runs resource checks + final output checksum check)
- `make verify_resources` (checks `resources/SHA256SUMS.txt`)

Toolchain and QC helpers:

- `make versions` (quick local tool-version report)
- `make toolchain_lock` (writes `doc/tool_versions.lock.txt`)
- `make qc_validate` (validate manifest inputs)
- `make qc_metrics` (compute sequencing metrics TSV under `results/qc/<QC_RUN>/`)
- `make clean_qc QC_RUN=<label>` (remove one QC run directory)
- `make print_qc_paths` (print resolved QC file paths)
- `make preprocessing_preflight` (wrapper for `src/pipelines/preflight_preprocessing.sh`)
- `make preprocessing_dry` (wrapper for `DRY_RUN=1 src/pipelines/preprocessing.sh`)
- `make preprocessing_run` (wrapper for `DRY_RUN=0 src/pipelines/preprocessing.sh`)

Target implementations are in the repository root `Makefile`.

## Raw Data Placement (Git-Ignored)

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
Raw data can be obtained for academic purposes upon reasonable request to the repository owner.
