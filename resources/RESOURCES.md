
This directory contains *small, redistributable* reference files vendored into the repository to enable
reproducible runs without relying on machine-specific paths.

## What is vendored here (small files)

Typical examples:
- BED files (targets, PRNP coding region)
- Interval lists
- Index/dictionary sidecars (e.g. `.fai`, `.dict`) when small enough

## Large resources (NOT stored in Git)

The following are required for some pipelines but are too large to commit directly to Git.
Obtain them separately and place them in a local, user-specific path (configured via
the active env config files under `config/`, for example `preprocessing.env` and `junctions.env`)
or use Git LFS if you explicitly choose to version them.

### Reference FASTA
- Expected file(s): `chr2_chr4_chr20.fasta` and BWA index sidecars (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`), plus `.fai` and `.dict`
- Purpose: alignment / variant calling reference
- How to obtain:
  - Provide the exact source of your FASTA (e.g., hg38 build, subset creation procedure).
  - If you subset the genome, document the command(s) used (samtools faidx / seqkit / custom script).
- Recommended verification:
  - Record SHA-256 checksums for the FASTA and indices.

### gnomAD �af-only� resource VCF
- Expected file(s): `somatic-hg38_af-only-gnomad.hg38.vcf.gz` (+ `.tbi` and possibly `.idx`)
- Purpose: Mutect2 germline resource
- How to obtain:
  - Document the upstream source and version (URL or release tag) and any processing you performed.

### Panel of Normals (PoN)
- Expected file(s): `somatic-hg38_1000g_pon.hg38.vcf.gz` (+ `.tbi`) or your custom PoN
- Purpose: Mutect2 PoN filtering
- How to obtain:
  - If using Broad-provided PoN: document which one (hg38) and source.
  - If custom: document how it was generated (inputs, commands, date).

### dbSNP VCF
- Expected file(s): `dbsnp_146.hg38.vcf.gz` (+ `.tbi`)
- Purpose: annotation / known sites (depending on your pipeline)
- How to obtain:
  - Document source (e.g., NCBI dbSNP build number) and genome build (hg38).

### Mills and 1000G indels VCF
- Expected file(s): `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` (+ `.tbi`)
- Purpose: known indels (often for BQSR / filtering depending on workflow)
- How to obtain:
  - Document source and version.

### GTF
- Expected file(s): `Homo_sapiens.GRCh38.110.gtf.gz`
- Purpose: annotation / plotting / transcript coordinates
- How to obtain:
  - Document source (e.g., Ensembl release 110) and genome build.

## Suggested local layout (outside Git)

Example (not prescriptive):
- `data/ref/` (local, not committed)
  - large FASTA + indices
  - large VCFs + indices
  - large GTFs

Point your `config/*.env` values to these local paths.

## Checksums

For reproducibility, record SHA-256 checksums for all large resources.
Example command:

    sha256sum <file> > <file>.sha256
