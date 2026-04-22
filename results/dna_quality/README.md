# DNA Quality Results

Derived outputs from the DNA-quality proxy workflow are written under `results/dna_quality/<RUN_ID>/`.

Current development command:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File src\dna_quality\run_dna_quality_analysis.ps1 -OutputRun latest
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

The `latest/` directory is convenient for iteration. For archival reruns, use a timestamped `-OutputRun`.
