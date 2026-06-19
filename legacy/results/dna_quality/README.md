# DNA Quality Results

Derived outputs from the DNA-quality proxy workflow are written under `results/dna_quality/<RUN_ID>/`.

Current maintained command:

```bash
make dna_quality
```

Default raw inputs are read from the gitignored repo-local `raw/dna_quality/` tree:

- `raw/dna_quality/sureselect`
- `raw/dna_quality/ddpcr`
- `raw/dna_quality/samples`

Key outputs:

- `sample_quality_evidence_table.tsv`
