# Finalization Checklist

Use this checklist immediately before final reproducibility commit/tag.

## Inventory

- Regenerate active script inventory:

```bash
python3 doc/reproducibility/make_inventory.py
```

- Review and stage inventory update:
  - `doc/inventory.tsv`

## Output verification

- Recompute checksums (if any outputs changed):

```bash
bash doc/reproducibility/verify_output_checksums.sh --mode write
```

- Verify checksums:

```bash
bash doc/reproducibility/verify_output_checksums.sh --mode check
```

## Documentation consistency

- Ensure root `README.md` run steps match current script names and order.
- Ensure workflow READMEs (`src/ddpcr`, `src/junctions`, `src/pipelines`) reflect current inputs/outputs.
- In `src/pipelines/README.md`, ensure stage output lists match real script outputs (and do not list per-stage log files unless scripts actually emit them).
- Ensure provenance/version notes in `doc/reproducibility/tooling_and_reference_provenance.md` are up to date.

Optional version sanity check (against reproducibility envs):

```bash
conda run -n prnp-somatic python --version
conda run -n prnp-somatic samtools --version | sed -n '1p'
conda run -n prnp-somatic bedtools --version
conda run -n prnp-somatic tabix --version | sed -n '1p'
conda run -n prnp-somatic gatk --version | sed -n '1,4p'
conda run -n prnp-junctions Rscript --version 2>&1 | sed -n '1p'
conda run -n prnp-somatic-ddpcr Rscript --version 2>&1 | sed -n '1p'
```
