# Manuscript build scaffold

This folder contains scripts and configuration to generate manuscript figures and tables from repository outputs.

## Layout

- `run_all.R`: entrypoint to generate all manuscript artifacts
- `scripts/`: figure/table scripts
- `config/`: manuscript-specific constants and mappings
- `figures/`: generated figure files (gitignored)
- `tables/`: generated table files (gitignored)

## Inputs

Scripts are expected to read from repository outputs (for example, `results/` and `ddPCR/`) using repo-relative paths.

## Run

```bash
Rscript manuscript/run_all.R
```

## Policy

Generated manuscript artifacts are not tracked in git by default.
