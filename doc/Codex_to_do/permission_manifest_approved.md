# Permission Manifest Used (Approved)

Approved by user on 2026-02-25 ("ok go"), constrained to the preflight scope.

## Write scope used

- `src/ddPCR/*`
- `src/junctions/*`
- `src/pipelines/*` (including temporary AAF-filter override in `12_cjd_dilutions_variant_table_qc_with_pon.R`)
- `doc/Codex_to_do/*`
- `doc/reproducibility/*`
- `bin/*`
- `README.md`

No destructive git commands were used. No raw data were moved/deleted.

## Command scope used

- Read/discovery: `ls`, `find`, `rg`, `sed`, `cat`, `git status`
- Dependency checks: `command -v ...`, `Rscript -e ...`
- Workflow run:
  - `bash src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh`
- Verification:
  - `bash bin/verify_output_checksums.sh --mode write`
  - `bash bin/verify_output_checksums.sh --mode check`

## Network/install policy used

- No dependency installations were performed.
- No network-dependent install commands were executed.

## Escalation

- No escalated command execution was required for this pass.
