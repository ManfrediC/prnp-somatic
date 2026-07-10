# Reproducibility

This directory contains the small, tracked files used to check and document the
repository's reproducibility state. It is for manifests, checksums and
provenance notes, not primary analysis code or private reviewer-response
material.

## Key Files

- `tooling_and_reference_provenance.md`: required tools, reference resources
  and package-context notes.
- `final_outputs_manifest.tsv`: final outputs expected from the workflows.
- `final_outputs.sha256`: checksums for declared final outputs.
- `sra_accessions.tsv`: processed SRA study/BioSample/run accession map.
- `data_availability.md`: reviewer-facing data availability summary.
- `inventory.tsv`: maintainer-side active script inventory under `src/`.
- `make_inventory.py`: optional helper for regenerating `inventory.tsv`.

## Utilities

Run commands from the repository root.

Verify final-output checksums:

```bash
bash reproducibility/verify_output_checksums.sh --mode check
```

Refresh checksums only after intentionally changing final outputs:

```bash
bash reproducibility/verify_output_checksums.sh --mode write
```

Regenerate the maintainer inventory when active scripts change:

```bash
python3 reproducibility/make_inventory.py
```
