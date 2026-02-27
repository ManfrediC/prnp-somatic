# Resources

This directory contains reference files used by the reproducible workflows.

Use `resources/INDEX.tsv` as the canonical inventory of resource files, workflow scope, provenance notes, and tracking status.

## Policy

- Keep script paths stable under `resources/`.
- Treat files in `resources/junctions/` (`prnp_junctions.fa`, `PRNP.pad1kb.hg38.bed`) as committed reference artefacts for the junction workflow.
- Treat BWA index sidecars for `resources/junctions/prnp_junctions.fa` (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) as runtime-generated files.
- Configure pipeline scripts via `config/preprocessing.env` and `config/junctions.env` when local overrides are needed.

## Canonical paths used by current workflows

- SNV pipeline subset reference FASTA: `resources/chr2_chr4_chr20.fasta`
- SNV pipeline dbSNP: `resources/dbsnp_146.hg38.vcf.gz`
- SNV pipeline Mills indels: `resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`
- SNV pipeline gnomAD AF-only: `resources/somatic-hg38_af-only-gnomad.hg38.vcf.gz`
- SNV pipeline manual annotation table: `resources/annotations/manual_population_freq.tsv`
- SNV pipeline Funcotator datasource root: `resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38`
- Junction workflow GTF: `resources/Homo_sapiens.GRCh38.110.gtf.gz`
- Junction workflow reference FASTA: `resources/hg38.fa`
- Junction workflow generated reference files: `resources/junctions/prnp_junctions.fa`, `resources/junctions/PRNP.pad1kb.hg38.bed`

## Canonical file forms

- Keep `hg38.fa` as the canonical junction FASTA path.
- Use `somatic-hg38_1000g_pon.hg38.vcf.gz` (plus index) for PoN resource VCF usage.

## Integrity checks

From repository root:

```bash
make verify_resources
```

This validates checksums listed in `resources/SHA256SUMS.txt`.
