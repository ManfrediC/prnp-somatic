# prnp-somatic

This repository contains the reproducible analysis workflows used in the *PRNP* somatic project:

- NGS somatic SNV detection (`src/pipelines`): FASTQ preprocessing, Mutect2 calling, post-processing, QC, and final variant tables.
- ddPCR SNV quantification (`src/ddpcr`): processing of raw droplet exports into long-format and participant-level result tables.
- PRNP exon-exon junction analysis (`src/junctions`): junction-reference construction, read realignment, and junction-count quantification.
- PRNP octapeptide repeat region analysis (`src/repeats`): repeat-aware screening of OPRI/OPRD candidates from existing BAMs.
- Manuscript artifact generation (`manuscript`): scripts that build figures/tables from outputs in `results/`.

The Reproduction Guide below details required inputs, environment setup and expected outputs.

## Repository structure
- `src/` reusable code (Shell/Python/R)
- `config/` configuration templates and sample manifests
- `reproducibility/` final-output manifests, checksums, provenance notes and repository maintenance checks
- `resources/` small static artefacts tracked in git (BEDs, schemas). Larger resources must be downloaded manually
- `results/` outputs *not tracked* (date-stamped runs)
- `env/` container/environment definitions
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

- `reproducibility/tooling_and_reference_provenance.md`

This document lists:

- required reference resources and their source/provenance
- required CLI/R tools and version context
- environment-specific package notes (including junction and ddPCR dependencies)

The remaining reproducibility manifests and checks are summarised in:

- `reproducibility/README.md`

Funcotator resource acquisition, layout and datasource details are documented in:

- `resources/README.md`

### 3. Prepare input data in expected locations

Required inputs must be placed as follows:

- ddPCR reviewed-gate raw exports: `raw/ddpcr/{ddpcr_archive,csv_export,archive_contents,layout_xlsx,manifests}/`
- ddPCR original-gate provenance exports: `raw/ddpcr_original/`
- ddPCR sample metadata: `raw/ddpcr/sample_details.xlsx`
- raw FASTQ files: `raw/fastq`
- resources for GATK workflow: `resources`

See `README` files in the respective directories for details on the required files.

Data can be obtained for academic purposes upon reasonable request to the
repository owner.

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
bash src/ddpcr/run_ddpcr.sh
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

### 9. Create environment for repeat workflow

```bash
conda env create -f env/repeats.environment.yml
conda activate prnp-repeats
```

### 10. Run repeat workflow (`src/repeats`)

Set up repeat configuration first:

```bash
cp config/repeats.env.example config/repeats.env
```

Then edit `config/repeats.env` for local paths/settings as documented in
`config/README.md`.

Main cohort run:

```bash
conda activate prnp-repeats
bash src/repeats/01_run_prnp_orr.sh
```

If you are rerunning into an existing `results/repeats/` directory, set
`ARCHIVE_EXISTING_RUN=1` in `config/repeats.env` first so the current live run
is archived instead of mixed with new outputs.

Expected outputs:

- `results/repeats/sample_calls.tsv`
- `results/repeats/candidate_calls.tsv`
- `results/repeats/cohort_summary.tsv`
- `results/repeats/subclonal_read_support.tsv`
- `results/repeats/gangstr_calls.tsv`
- `results/repeats/somatic_screen.tsv`
- `results/repeats/run_settings.tsv`

Optional one-sample manual mosaic review:

```bash
conda activate prnp-repeats
python3 src/repeats/06_manual_mosaic_prnp_orr.py \
  --bam results/final_bam/CJD1.bam \
  --reference-fasta resources/references/snv/chr2_chr4_chr20.fasta \
  --output-prefix results/repeats/manual/CJD1 \
  --sample-calls-tsv results/repeats/sample_calls.tsv \
  --gangstr-calls-tsv results/repeats/gangstr_calls.tsv
```

This writes `results/repeats/manual/CJD1.reads.tsv` and
`results/repeats/manual/CJD1.summary.tsv`.

For control-first cohort calibration:

```bash
conda activate prnp-repeats
make repeats_manual_controls
make repeats_manual_cjd
```

These write cohort summaries plus per-sample review TSVs under
`results/repeats/manual_cohort/controls/` and
`results/repeats/manual_cohort/cjd/`.

To apply the default background-aware post-processing filter to the CJD cohort:

```bash
conda activate prnp-repeats
make repeats_manual_filter_cjd
```

This writes:

- `results/repeats/manual_cohort/cjd/filtered/default.annotated.tsv`
- `results/repeats/manual_cohort/cjd/filtered/default.candidates.tsv`
- `results/repeats/manual_cohort/cjd/filtered/default.overview.tsv`

### 11. Run junction workflow

```bash
TMPDIR=/tmp TEMP=/tmp TMP=/tmp bash src/junctions/run_junctions.sh
```

Expected outputs:

- `results/junctions/junction_counts/prnp_junction_counts.tsv`
- `results/junctions/junction_counts/prnp_junction_summary.tsv`

### 12. Verify final outputs

```bash
bash reproducibility/verify_output_checksums.sh --mode check
```

Expected manifest/checksum files:

- `reproducibility/final_outputs_manifest.tsv`
- `reproducibility/final_outputs.sha256`

Detailed workflow docs:

- `src/ddpcr/README.md`
- `src/junctions/README.md`
- `src/pipelines/README.md`
- `manuscript/README.md`

### 13. Makefile commands (optional wrappers)

As an alternative to executing the pipelines individually, it's possible to run each part of the project using `Makefile`.

Show all available commands:

```bash
make help
```

Run main workflows:

- `make ddpcr` (requires active env: `prnp-somatic-ddpcr`)
- `make snv` (requires active env: `prnp-somatic`; runs the Stage-12 publication-path wrapper)
- `make repeats` (requires active env: `prnp-repeats`)
- `make repeats_manual_controls` (requires active env: `prnp-repeats`)
- `make repeats_manual_cjd` (requires active env: `prnp-repeats`)
- `make repeats_manual_filter_cjd` (requires active env: `prnp-repeats`)
- `make junctions` (requires active env: `prnp-junctions`)
- `make all` (runs ddPCR + SNV Stage-12 wrapper + junctions via `conda run`; expects envs `prnp-somatic-ddpcr`, `prnp-somatic`, and `prnp-junctions`)

Run integrity checks:

- `make check` (requires active env: `prnp-somatic`; runs resource checks + final output checksum check)
- `make verify_resources` (checks `resources/SHA256SUMS.txt`)

Toolchain and QC helpers:

- `make versions` (quick local tool-version report)
- `make toolchain_lock` (writes `reproducibility/tool_versions.lock.txt`)
- `make qc_metrics` (compute sequencing metrics TSV under `results/sequencing_qc/`)
- `make dna_quality` (build DNA-quality evidence table under `results/dna_quality/`)
- `make clean_qc` (remove the canonical sequencing QC output directory)
- `make print_qc_paths` (print resolved QC file paths)
- `make preprocessing_preflight` (wrapper for `src/pipelines/preflight_preprocessing.sh`)
- `make preprocessing_dry` (wrapper for `DRY_RUN=1 src/pipelines/preprocessing.sh`)
- `make preprocessing_run` (wrapper for `DRY_RUN=0 src/pipelines/preprocessing.sh`)

Target implementations are in the repository root `Makefile`.

## Raw Data Placement (Git-Ignored)

Raw data and sequencing pipeline outputs are gitignored, to account for GitHub storage constraints.

- DNA quality raw inputs: `raw/dna_quality/`
- ddPCR raw exports: `raw/ddpcr/` (ignored via repo `.gitignore` except `README.md`)
- ddPCR metadata sheet: `raw/ddpcr/sample_details.xlsx` (ignored via repo `.gitignore`)
- junction BAM inputs: `results/final_bam/*.bam` (ignored via repo `.gitignore` patterns for BAM/results)
- pipeline run intermediates: `runs/**` (ignored via repo `.gitignore`)
- pipeline/junction outputs: `results/**` (ignored except tracked placeholders)

Generated runtime index sidecars are also ignored (for example `resources/junctions/*.fa.{amb,ann,bwt,pac,sa}`).

## Data availability

Sequencing FASTQ files have been uploaded to NCBI SRA.

- BioProject: `PRJNA1484136`
  (<http://www.ncbi.nlm.nih.gov/bioproject/1484136>)
- SRA submission: `SUB16298068`
- BioSample accessions: `SAMN61220231` through `SAMN61220265` for the 35
  submitted samples
- Release date: 2027-07-30, or with release of linked data
- Per-sample accession and FASTQ filename map:
  `reproducibility/sra_accessions.tsv`
- SRA experiment and run accessions: pending; these will be added to the
  accession map once assigned
