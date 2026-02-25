# bin utilities

Helper scripts for repository maintenance and reproducibility checks.

## Scripts

- `verify_output_checksums.sh`
  - Purpose: generate or verify SHA256 checksums for expected final outputs listed in `doc/reproducibility/final_outputs_manifest.tsv`.
  - Modes:
    - `--mode write`: writes checksum file (default path: `doc/reproducibility/final_outputs.sha256`)
    - `--mode check`: verifies current outputs against checksum file
  - Usage:
    - `bash bin/verify_output_checksums.sh --mode write`
    - `bash bin/verify_output_checksums.sh --mode check`
  - Optional arguments:
    - `--manifest <path>`
    - `--checksums <path>`

- `make_inventory.py`
  - Purpose: create a TSV inventory of script files under `src/legacy/`.
  - Output: `doc/inventory.tsv`
  - Usage:
    - `python3 bin/make_inventory.py`

## Notes

- Run commands from repository root.
- `verify_output_checksums.sh` fails if any `required=yes` entry in the manifest is missing.
