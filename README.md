# prnp-somatic

PRNP somatic variant analysis:
FASTQ -> BAM (GATK best practices) -> Mutect2 -> QC -> result tables/figures.

## Quickstart
1. Prepare inputs (FASTQ + manifest)
2. Configure paths in `config/config.local.yaml` (not committed)
3. Run: `make all`

## Controls Variant Calling (No PoN)

For controls-only somatic variant calling without a panel of normals:

1. Run Stage 1 Mutect2: `src/pipelines/mutect2_controls_no_pon.sh`
2. Run Stages 2-6 post-processing: `src/pipelines/mutect2_controls_postprocess_no_pon.sh`

Pipeline details and outputs are documented in `src/pipelines/README.md`.

## Repository structure
- `bin/` command-line entrypoints / wrappers
- `src/` reusable code (Python/R)
- `config/` configuration templates
- `resources/` small static artefacts tracked in git (BEDs, schemas)
- `data/` input data *not tracked* (see `data/README.md`)
- `results/` outputs *not tracked* (date-stamped runs)
- `env/` container/environment definitions
- `doc/` notes and documentation
- `tests/` smoke tests

## Conda environment
- Analyses were performed using the `Conda` environment defined in `env/environment.yml`.

## Data availability
FASTQs are hosted on xxx.
