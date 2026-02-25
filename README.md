# prnp-somatic

PRNP somatic variant analysis:
FASTQ -> BAM (GATK best practices) -> Mutect2 -> QC -> result tables/figures.

## Reproducible Workflow Entrypoints (repo root)

Run each workflow from the repository root with one command:

- ddPCR SNV dataframe: `bash src/ddPCR/run_ddpcr.sh`
- PRNP exon-exon junctions: `TMPDIR=/tmp bash src/junctions/run_junctions.sh`
- Stage-12 CJD+dilutions QC table regeneration: `bash src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh`
- Manuscript figures/tables: `Rscript manuscript/run_all.R`

Detailed inputs and expected outputs for each workflow are documented in:

- `src/ddPCR/README.md`
- `src/junctions/README.md`
- `src/pipelines/README.md`
- `manuscript/README.md`

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
- `resources/` small static artefacts tracked in git (BEDs, schemas)
- `authoritative_files/` manifests and sequencing-metrics validation utilities
- `data/` input data *not tracked* (see `data/README.md`)
- `results/` outputs *not tracked* (date-stamped runs)
- `env/` container/environment definitions
- `doc/` notes and documentation
- `tests/` smoke tests

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
