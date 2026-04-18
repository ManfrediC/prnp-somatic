# Environment specifications

This directory contains Conda environment specifications used by the reproducible workflows in this repository.

## Files

- `environment.yml`:
  - baseline `prnp-somatic` environment specification for the SNV pipeline and shared tooling.

- `environment.lock.yml`:
  - resolved lockfile for `environment.yml`.

- `junctions.environment.yml`:
  - dedicated `prnp-junctions` environment specification for the exon-junction workflow (`src/junctions`).

- `ddpcr.environment.yml`:
  - dedicated `prnp-somatic-ddpcr` environment specification for the ddPCR workflow (`src/ddPCR`).

- `repeats.environment.yml`:
  - dedicated `prnp-repeats` environment specification for ExpansionHunter, GangSTR, and cohort summarization in the PRNP ORR workflow (`src/repeats`).

- `reviewer.environment.yml`:
  - dedicated `prnp-reviewer` environment specification for REViewer in the PRNP ORR workflow (`src/repeats`).

## Usage

Run from repository root.

Create baseline environment:

```bash
conda env create -f env/environment.yml
conda activate prnp-somatic
```

Create baseline environment from lockfile:

```bash
conda env create -f env/environment.lock.yml
conda activate prnp-somatic
```

Create junction workflow environment:

```bash
conda env create -f env/junctions.environment.yml
conda activate prnp-junctions
```

Create ddPCR workflow environment:

```bash
conda env create -f env/ddpcr.environment.yml
conda activate prnp-somatic-ddpcr
```

Create repeat workflow environment:

```bash
conda env create -f env/repeats.environment.yml
conda activate prnp-repeats
```

Create REViewer environment:

```bash
conda env create -f env/reviewer.environment.yml
```

## Notes

- Environment creation can vary slightly across Conda versions and channels.
