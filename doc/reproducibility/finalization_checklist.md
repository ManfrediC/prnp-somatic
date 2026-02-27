# Finalization Checklist

Use this checklist immediately before final reproducibility commit/tag.

## Inventory

- Regenerate active script inventory:

```bash
python3 bin/make_inventory.py
```

- Review and stage inventory update:
  - `doc/inventory.tsv`

## Output verification

- Recompute checksums (if any outputs changed):

```bash
bash bin/verify_output_checksums.sh --mode write
```

- Verify checksums:

```bash
bash bin/verify_output_checksums.sh --mode check
```

## Documentation consistency

- Ensure root `README.md` run steps match current script names and order.
- Ensure workflow READMEs (`src/ddPCR`, `src/junctions`, `src/pipelines`) reflect current inputs/outputs.
- Ensure provenance/version notes in `doc/reproducibility/tooling_and_reference_provenance.md` are up to date.
