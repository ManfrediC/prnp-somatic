# Config

Committed:
- `preprocessing_samples.tsv` — sample manifest for preprocessing (batch_id, sample_id, R1, R2)
- `preprocessing.env.example` — template for local configuration

Not committed (local machine only):
- `preprocessing.env` — local paths and settings (ignored by `.gitignore`)

Rationale: `preprocessing.env` may include machine-specific paths and runtime settings; large data directories (`fastq/`, `runs/`, `results/final_bam/`) are intentionally not committed.