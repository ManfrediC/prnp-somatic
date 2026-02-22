# Config

Committed:
- `preprocessing_samples.tsv` - sample manifest for preprocessing (`batch_id`, `sample_id`, `R1`, `R2`)
- `preprocessing.env.example` - template for local configuration

Not committed (local machine only):
- `preprocessing.env` - local paths and settings (ignored by `.gitignore`)

Rationale: `preprocessing.env` may include machine-specific paths and runtime settings; large data directories (`fastq/`, `runs/`, `results/final_bam/`) are intentionally not committed.

## Funcotator-related settings

For controls post-processing and annotation, set these in `preprocessing.env`:

- `FUNCOTATOR_DS`
- `REF_FASTA`

Recommended `FUNCOTATOR_DS` value:

- `resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38`
