# Makefile Usage

This repository includes a root `Makefile` for QC and preprocessing utilities.
Run all commands from the repository root (your cloned `prnp-somatic` directory).

## Prerequisites

- `make` installed.
- Conda environment `prnp-somatic` activated for targets that depend on `check_conda`.

```bash
conda activate prnp-somatic
```

## Quick commands

- Show available commands:

```bash
make help
```

- Print tool versions:

```bash
make versions
```

- Lock versions to file:

```bash
make toolchain_lock
```

- Show resolved QC paths:

```bash
make print_qc_paths
```

## QC targets

- Validate manifests and write a log:

```bash
make qc_validate
```

- Compute sequencing metrics TSV (runs `qc_validate` first):

```bash
make qc_metrics
```

- Use a custom run label:

```bash
make qc_metrics QC_RUN=2026-02-25_review
```

- Remove one QC run folder:

```bash
make clean_qc QC_RUN=2026-02-25_review
```

## Resource checksum verification

```bash
make verify_resources
```

## Preprocessing helper targets

- Preflight checks:

```bash
make preprocessing_preflight
```

- Dry run:

```bash
make preprocessing_dry
```

- Real run:

```bash
make preprocessing_run
```

## Variables you can override

- `QC_RUN` (default: `latest`)
- `RESULTS_DIR` (default: `results`)
- `AUTH_DIR` (default: `authoritative_files`)

Example:

```bash
make qc_metrics RESULTS_DIR=results QC_RUN=my_test
```

## Notes

- Reproducible workflow entrypoints are still the wrapper scripts in `src/`:
  - `bash src/ddPCR/run_ddpcr.sh`
  - `TMPDIR=/tmp bash src/junctions/run_junctions.sh`
  - `bash src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh`
- The modular make fragments included by root `Makefile` are documented in `mk/README.md`.
