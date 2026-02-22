# Config

Committed:
- `preprocessing_samples.tsv` - sample manifest for preprocessing (`batch_id`, `sample_id`, `R1`, `R2`)
- `preprocessing.env.example` - template for local configuration

Not committed (local machine only):
- `preprocessing.env` - local paths and settings (ignored by `.gitignore`)

Rationale: `preprocessing.env` may include machine-specific paths and runtime settings; large data directories (`fastq/`, `runs/`, `results/final_bam/`) are intentionally not committed.

## Controls post-processing settings

For controls post-processing and annotation, set these in `preprocessing.env`:

- `FUNCOTATOR_DS`
- `REF_FASTA`
- `GNOMAD_AF_VCF`

Recommended `FUNCOTATOR_DS` value:

- `resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38`

Recommended `GNOMAD_AF_VCF` value:

- `resources/somatic-hg38_af-only-gnomad.hg38.vcf.gz`

The corresponding index (`.tbi` or `.csi`) must also be present.

For controls post-processing, `gnomAD_exome` and `gnomAD_genome` should remain excluded from
the active datasource tree to avoid requester-pays access requirements.
Archived copies are kept under:

- `resources/backup/funcotator_excluded_datasources/`
