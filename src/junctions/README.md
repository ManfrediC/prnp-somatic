# PRNP Exon-Exon Junction Scripts

This directory contains the reproducible exon-exon junction workflow scripts used in this repository.
Paths are repo-relative and resolve from repository root.

## Scripts

- `01_build_prnp_junction_fasta.R`
- `02_make_prnp_bed.R`
- `03_process_bam.sh`
- `04_count_prnp_junctions.R`

## Inputs

Defaults are repo-relative and resolved from script location.

- `resources/Homo_sapiens.GRCh38.110.gtf.gz`
- `resources/hg38.fa`
- `results/final_bam/*.bam` (for `03_process_bam.sh`)

## Command

From repository root (single entrypoint):

```bash
conda activate prnp-junctions
TMPDIR=/tmp TEMP=/tmp TMP=/tmp bash src/junctions/run_junctions.sh
```

Stepwise commands:

```bash
Rscript src/junctions/01_build_prnp_junction_fasta.R
Rscript src/junctions/02_make_prnp_bed.R
bash src/junctions/03_process_bam.sh
Rscript src/junctions/04_count_prnp_junctions.R
```

## Outputs

- `resources/junctions/prnp_junctions.fa`
- `resources/junctions/PRNP.pad1kb.hg38.bed`
- `results/junctions/junction_align/`
- `results/junctions/junction_counts/`

Runtime-generated BWA index sidecars for the FASTA are created in `resources/junctions/` (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) and should stay untracked.

Final expected count outputs:

- `results/junctions/junction_counts/prnp_junction_counts.tsv`
- `results/junctions/junction_counts/prnp_junction_summary.tsv`

These expected output paths are also listed in:

- `doc/reproducibility/final_outputs_manifest.tsv`

## Config Overrides

You can override defaults via environment variables:

- `SAMPLEDIR`
- `BED`
- `JUNC_FA`
- `OUTDIR`
- `THREADS`
- `PRNP_JUNCTION_GTF`
- `PRNP_JUNCTION_ALIGN_DIR`
- `PRNP_JUNCTION_COUNT_DIR`

Example:

```bash
cp config/junctions.env.example config/junctions.env
bash src/junctions/run_junctions.sh
```

## R Package Check

The junction R scripts require:

- `GenomicFeatures`
- `GenomicRanges`
- `GenomeInfoDb`
- `Rsamtools`
- `Biostrings`
- `rtracklayer`
- `GenomicAlignments`
- `dplyr`
- `txdbmaker`

Quick check:

```bash
Rscript -e 'pkgs <- c("GenomicFeatures","GenomicRanges","GenomeInfoDb","Rsamtools","Biostrings","rtracklayer","GenomicAlignments","dplyr","txdbmaker"); print(pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)])'
```
