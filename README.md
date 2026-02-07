# prnp-somatic

PRNP somatic variant analysis:
FASTQ -> BAM (GATK best practices) -> Mutect2 -> QC -> result tables/figures.

## Quickstart
1. Prepare inputs (FASTQ + manifest)
2. Configure paths in `config/config.local.yaml` (not committed)
3. Run: `make all`

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

## Data availability
FASTQs are hosted on xxx.
