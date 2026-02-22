
This directory contains reference resources and interval definitions used by the PRNP somatic pipeline.

**Policy**
- Small, hand-authored files (e.g. BEDs, interval lists, this README) are versioned in Git.
- Large reference datasets (FASTA/VCF/GTF and their large indices) are **not** versioned in Git and must be obtained separately.

Where possible, the pipeline should refer to resources via **relative paths** under `resources/`.

---

## Files tracked in this repository (small)

These files are expected to be present after cloning the repo:

### Capture / target definitions
- `targets.bed`  
  Capture target intervals (BED).
- `capture_targets.interval_list`  
  Target intervals in Picard/GATK interval_list format.
- `prnp_coding.bed`  
  BED defining the PRNP coding region used for coverage summaries.

### Minimal manifests / example inputs
- `dataset.tsv`  

### Small reference sidecar files
These may be tracked because they are small and sometimes helpful for validation, but are not sufficient alone without the corresponding large reference:
- `hg38.dict`
- `hg38.fa.fai`
- `chr2_chr4_chr20.dict`
- `chr2_chr4_chr20.fasta.fai`
- `chr2_chr4_chr20.fasta.amb`
- `chr2_chr4_chr20.fasta.ann`

---

## Large files not tracked (must be obtained separately)

The following are typically required for full execution, but are **not** stored in Git due to size/licensing:

### Reference genomes
- `hg38.fa` (GRCh38 reference FASTA) and index files  
- `chr2_chr4_chr20.fasta` (subset FASTA used in some workflows) and BWA indices  
  e.g. `.bwt`, `.pac`, `.sa` (and optionally `.amb/.ann/.fai`)

### Annotation / population resources (from the GATK resource bundle)
- `Homo_sapiens.GRCh38.110.gtf.gz` (Ensembl GTF; version 110)
- `dbsnp_146.hg38.vcf.gz` (+ `.tbi`)
- `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` (+ `.tbi`)
- `somatic-hg38_af-only-gnomad.hg38.vcf.gz` (+ `.tbi` and/or `.idx`)
- `somatic-hg38_1000g_pon.hg38.vcf.gz` (+ `.tbi`) (panel of normals; if used)

### Funcotator datasource bundle (required for annotation stages)
- `funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz` (archive)
- `funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38/` (extracted datasource root)

The pipeline expects `FUNCOTATOR_DS` to point to:

- `resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38`

Typical datasource subfolders include:

- `clinvar/`
- `dbsnp/`
- `gencode/`
- `gnomAD_exome/`
- `gnomAD_genome/`
- `hgnc/`

---

## How to obtain large resources

This repository intentionally does not bundle large reference files. Large reference files can be obtained from the GATK resource bundle, Ensembl and gnomAD.

After downloading, place the files into `resources/` using the filenames listed above so that relative paths remain stable.

---

## Integrity checks (recommended)

After placing resources, record and verify checksums:

```bash
cd resources
sha256sum * > SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt
