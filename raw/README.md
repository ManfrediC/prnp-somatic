# Raw data

This directory is the repository-local home for raw input data. Raw data are not
tracked in git.

Expected layout:

- `raw/fastq/`: raw FASTQ files, organised by sequencing batch.
- `raw/ddpcr/`: reviewed-gate ddPCR raw database and `sample_details.xlsx`.
- `raw/ddpcr_original/`: original-gate ddPCR export provenance.
- `raw/dna_quality/`: DNA-quality source files for Tapestation, ddPCR quantity
  records, sample manifests, and related metadata.

Derived outputs remain outside this directory, under `runs/` and `results/`.
