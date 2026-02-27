# bin utilities

This directory contains small maintenance utilities used for reproducibility checks and repository inventory.

## What This Directory Is For

- verify final-output integrity against a declared manifest
- regenerate the active script inventory used for documentation/curation

## Script Guide

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
