# Environment specifications

This directory contains Conda environment specifications used by the reproducible workflows in this repository.

## Files

- `environment.yml`:
  - baseline `prnp-somatic` environment specification for the SNV pipeline and shared tooling.

- `environment.lock.yml`:
  - resolved lockfile for `environment.yml`.

- `junctions.environment.yml`:
  - dedicated `prnp-junctions` environment specification for the exon-junction workflow (`src/junctions`).

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

## Notes

- Environment creation can vary slightly across Conda versions and channels.
