# Pipelines

This folder contains executable pipeline entry points. For preprocessing (FASTQ ? final BAM), use:

- `preflight_preprocessing.sh` (checks only; fail fast)
- `preprocessing.sh` (runs the pipeline; no preflight checks inside)

## Preprocessing overview

**Inputs**
- Sample sheet: `config/preprocessing_samples.tsv` (committed)
- Local config: `config/preprocessing.env` (gitignored)
- Raw FASTQs: `fastq/<batch_id>/...` (gitignored; large)
- Reference resources: `resources/` (committed; verified with SHA256)

**Outputs**
- Intermediates and logs: `runs/preprocessing/<batch_id>/<sample_id>/` (gitignored; large)
- Published final BAMs: `results/final_bam/<sample_id>.bam` + `.bam.bai` (gitignored; large)

## Quick start

1) Activate the Conda environment (recommended baseline):

```conda activate prnp-somatic```

2) Run "preflight" check (checks tools, references, sample sheet, FASTQ presence):

`src/pipelines/preflight_preprocessing.sh`

3) Preview what would run (no filesystem changes):

`DRY_RUN=1 src/pipelines/preprocessing.sh`

4) Run the pipeline:

`DRY_RUN=0 src/pipelines/preprocessing.sh`

## Batch selection

Both scripts contain a `BATCHES=(...)` block. Enable/disable batches by commenting lines, e.g.

```BATCHES=(
  CJD_16_samples
  # CJD_8_samples
)```

The sample sheet is filtered to these selected batch IDs.

## Behaviour and reproducibility

### Resume
The pipeline is designed to be re-runnable:

* If the expected output for a stage already exists and is non-empty, that stage is skipped.
* If `FORCE=1`, existing per-sample outputs are removed and re-created.

Example:

```bash
FORCE=1 DRY_RUN=0 src/pipelines/preprocessing.sh
```

### Logs and metadata

For each processed sample:

* Logs are written to: `runs/preprocessing/<batch>/<sample>/logs/`
* Minimal run metadata is written to: `runs/preprocessing/<batch>/<sample>/RUN_META.txt`

### Configuration file

The pipeline reads `config/preprocessing.env`. A template is provided as:

* `config/preprocessing.env.example` (committed)

Create your local config (gitignored):

```bash
cp -n config/preprocessing.env.example config/preprocessing.env
```

Key variables include:

* `THREADS`, `JAVA_MEM_GB`
* `REF_FASTA`, `REF_DICT`, `DBSNP_VCF`, `MILLS_VCF`, `ADAPTERS_FA`
* `RUNS_DIR`, `FINAL_BAM_DIR`, `SAMPLES_TSV`

## Troubleshooting

* If preflight fails with “missing FASTQ”: check `config/preprocessing_samples.tsv` paths and that FASTQs are present under `fastq/`.
* If a run stops mid-sample: re-run `preprocessing.sh` (it will skip completed stages).
* If you suspect a stale output is being reused: run with `FORCE=1` for that sample’s batch.
* Tool availability: preflight expects `bwa`, `samtools`, `gatk`, and (once enabled) also `trimmomatic` and `picard`.

````

---

## 2) Add this section to the repo root `README.md`

Insert somewhere near “Getting started”:

```markdown
## Preprocessing (FASTQ ? final BAM)

This repository can reproduce the preprocessing pipeline used in the PRNP project.

Directory conventions:
- Raw FASTQs: `fastq/` (not committed; large)
- Intermediates/logs: `runs/` (not committed; large)
- Final BAMs: `results/final_bam/` (not committed; large)
- Reference assets: `resources/` (committed; checksummed)

Workflow:
```bash
conda activate prnp-somatic
src/pipelines/preflight_preprocessing.sh
DRY_RUN=1 src/pipelines/preprocessing.sh
DRY_RUN=0 src/pipelines/preprocessing.sh
````

Configuration:

* Template: `config/preprocessing.env.example`
* Local (gitignored): `config/preprocessing.env`
* Sample sheet (committed): `config/preprocessing_samples.tsv`

Batch selection is controlled by toggling comment lines in the `BATCHES=(...)` block inside the pipeline scripts.