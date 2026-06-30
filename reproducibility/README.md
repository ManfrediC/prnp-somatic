# Reproducibility

This directory contains the small, tracked files used to check and document the
repository's reproducibility state. It is for manifests, checksums, provenance
notes and finalisation utilities, not primary analysis code.

## Key Files

- `tooling_and_reference_provenance.md`: required tools, reference resources
  and package-context notes.
- `final_outputs_manifest.tsv`: final outputs expected from the workflows.
- `final_outputs.sha256`: checksums for declared final outputs.
- `sra_accessions.tsv`: submitted SRA/BioSample accession map.
- `data_availability.md`: reviewer-facing data availability summary.
- `inventory.tsv`: active script inventory under `src/`.
- `finalization_checklist.md`: short pre-freeze checklist.

## Utilities

Run commands from the repository root.

Regenerate the active script inventory:

```bash
python3 reproducibility/make_inventory.py
```

Verify final-output checksums:

```bash
bash reproducibility/verify_output_checksums.sh --mode check
```

Refresh checksums only after intentionally changing final outputs:

```bash
bash reproducibility/verify_output_checksums.sh --mode write
```
