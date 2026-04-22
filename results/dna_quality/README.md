# DNA Quality Results

Derived outputs from the DNA-quality proxy workflow are written under `results/dna_quality/<RUN_ID>/`.

Current development command:

```bash
bash src/dna_quality/01_run_dna_quality_analysis.sh --output-run latest
```

Key outputs:

- `file_inventory.tsv`
- `library_qc.tsv`
- `prep_metadata.tsv`
- `input_dna_quantity.tsv`
- `sample_aliases.tsv`
- `sample_quality_master.tsv`
- `dna_quality_scorecard.tsv`
- `analysis_summary.tsv`
- `report.md`

The `latest/` directory is convenient for iteration. For archival reruns, use a timestamped `--output-run`.
