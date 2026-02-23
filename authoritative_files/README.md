# Authoritative Files

This directory contains the manifest and validation utilities used for sequencing-metrics generation.

- `manifest.tsv`: canonical sample manifest (`sample_id`, `group`, `batch`, `input_dir`).
- `validate_manifest.sh`: validates manifest rows and resolves BAM/BAI/Picard metrics paths.
- `manifest_qc.tsv`: validation output from `validate_manifest.sh`.
- `compute_sequencing_metrics.py`: computes per-sample sequencing metrics from the manifest.
- `sequencing_metrics_per_sample.schema.tsv`: output schema used by the metrics script.

## Current manifest convention

- `manifest.tsv` currently points `input_dir` to `results/final_bam` for all samples.
- BAMs are expected as `results/final_bam/<sample_id>.bam` with matching `.bai`.

## Refresh `manifest_qc.tsv`

Run:

```bash
cd authoritative_files
./validate_manifest.sh manifest.tsv manifest_qc.tsv
```
