# Config

Configuration files are repo-relative unless an absolute path is supplied.

## Tracked files

- `config.example.yaml`: generic configuration template.
- `config.local.yaml`: repository-local configuration snapshot.
- `preprocessing_samples.tsv`: preprocessing sample manifest (`batch_id`, `sample_id`, `r1`, `r2`, optional `display_label`).
- `fastq_filename_manifest.tsv`: reviewer-facing FASTQ aliases mapped back to original sequencer filenames and hashes.
- `dna_quality_sample_alias_overrides.tsv`: alias overrides for the DNA-quality evidence table.
- `preprocessing.env.example`: SNV/preprocessing and Mutect2 configuration template.
- `junctions.env.example`: exon-exon junction workflow configuration template.
- `repeats.env.example`: PRNP ORR repeat workflow configuration template.

## Local files

These are gitignored and machine-specific:

- `preprocessing.env`
- `junctions.env`
- `repeats.env`

Create them from the matching `.example` file only when local overrides are
needed:

```bash
cp config/preprocessing.env.example config/preprocessing.env
cp config/junctions.env.example config/junctions.env
cp config/repeats.env.example config/repeats.env
```

## Common paths

- SNV reference FASTA: `resources/chr2_chr4_chr20.fasta`
- Funcotator datasource root: `resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38`
- gnomAD AF-only VCF: `resources/somatic-hg38_af-only-gnomad.hg38.vcf.gz`
- manual population-frequency table: `resources/annotations/manual_population_freq.tsv`
- repeat BAM directory: `results/final_bam`
- repeat output directory: `results/repeats`

Use the example files for the full set of supported variables and defaults.
