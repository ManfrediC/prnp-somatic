# bin utilities

This directory contains small maintenance utilities used for reproducibility checks and repository inventory.

## What This Directory Is For

- verify final-output integrity against a declared manifest
- regenerate the active script inventory used for documentation/curation

## Script Guide

### `archive_large_files.sh`

Purpose:

- inventory large files that are good archive candidates
- compress archive-worthy files in place to reclaim disk space
- restore archived files back to their original names when a workflow needs them

Default scope:

- scans `runs`, `results`, `resources`, and `fastq`
- targets uncompressed large-file extensions: `bam`, `bai`, `sam`, `fa`, `fasta`, `fna`, `gtf`, `vcf`, `tsv`
- skips files that are already compressed, such as `*.fastq.gz`, `*.vcf.gz`, `*.tar.gz`, `*.cram`

Usage:

```bash
bash bin/archive_large_files.sh inventory
bash bin/archive_large_files.sh compress --dirs runs,results --exts bam,bai,sam
bash bin/archive_large_files.sh restore --dirs runs,results
```

Notes:

- default minimum size is `100 MiB`
- default archive format is `xz` because it typically shrinks BAMs and text references better than plain gzip
- restore operates on `*.xz` and `*.gz` files in the selected directories
- use `--dry-run` first to review what would change

Suggested use in this repository:

- archive intermediate BAMs under `runs/preprocessing`
- archive final BAMs under `results/final_bam` only if you do not need them immediately
- archive large uncompressed references under `resources/` when they are not in active use
- do not expect meaningful savings from files already ending in `.gz`

### `verify_output_checksums.sh`

Purpose:

- create or verify SHA256 checksums for outputs declared in `doc/reproducibility/final_outputs_manifest.tsv`

Default files:

- manifest: `doc/reproducibility/final_outputs_manifest.tsv`
- checksums: `doc/reproducibility/final_outputs.sha256`

Modes:

- `--mode write`
  - recomputes checksums and writes checksum file
- `--mode check`
  - recomputes checksums and compares with checksum file

Usage:

```bash
bash bin/verify_output_checksums.sh --mode write
bash bin/verify_output_checksums.sh --mode check
```

Optional arguments:

- `--manifest <path>`
- `--checksums <path>`

Failure behavior:

- exits non-zero if any `required=yes` output in the manifest is missing
- exits non-zero in check mode if checksum diff is non-empty

### `make_inventory.py`

Purpose:

- regenerate inventory of active (non-legacy) scripts under `src/`

Output:

- `doc/inventory.tsv`

Usage:

```bash
python3 bin/make_inventory.py
```

Current scope:

- includes script-like files in `src/`
- excludes `src/legacy/`

## Notes

- run these commands from repository root
- rerun `make_inventory.py` before final repo freeze/commit to capture latest script set
