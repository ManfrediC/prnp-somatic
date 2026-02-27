# Authoritative Files

This directory defines the canonical sample manifest and the QC/metrics utilities used by the
sequencing workflow. If someone needs to understand "which samples are expected" and "how
sample-level QC metrics are generated", this is the reference location.

## What This Directory Is For

- define the canonical sample inventory (`manifest.tsv`)
- validate whether expected BAM/BAI/metrics files exist (`validate_manifest.sh`)
- produce a validated manifest with resolved paths (`manifest_qc.tsv`)
- compute per-sample sequencing QC metrics in a stable schema (`compute_sequencing_metrics.py`)

## File Guide

- `manifest.tsv`
  - canonical sample table with columns:
    - `sample_id`
    - `group` (`dilution`, `CJD`, `control`)
    - `batch`
    - `input_dir`

- `validate_manifest.sh`
  - validates manifest structure and allowed group values
  - resolves BAM paths using supported naming styles:
    - `<sample_id>.bam` (`short`)
    - `<sample_id>.bwa.picard.markedDup.recal.bam` (`long`)
  - verifies `.bai` index presence
  - requires Picard metrics only for `long` BAM style
  - writes `manifest_qc.tsv`

- `manifest_qc.tsv`
  - generated validation report with resolved paths and per-row error flags
  - intended for quick inspection and fail-fast checks in automation

- `compute_sequencing_metrics.py`
  - computes read/depth/on-target metrics from BAM + BED inputs
  - emits TSV to stdout in locked schema order
  - uses:
    - `manifest.tsv`
    - `sequencing_metrics_per_sample.schema.tsv`
    - BED resources from `resources/`

- `sequencing_metrics_per_sample.schema.tsv`
  - defines output column order for metrics emission
  - keeps downstream parsing stable

## Typical Usage

Run from repository root:

```bash
bash authoritative_files/validate_manifest.sh \
  authoritative_files/manifest.tsv \
  authoritative_files/manifest_qc.tsv

python3 authoritative_files/compute_sequencing_metrics.py \
  > results/qc/latest/sequencing_metrics_per_sample.tsv
```

## Notes

- Current manifest convention points `input_dir` to `results/final_bam` for active samples.
- The validator fails when required files are missing, so it can be used as a preflight gate.
