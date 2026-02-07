#!/bin/bash
set -euo pipefail

# on 2025-04-21: No-PoN pipeline with MQ, strand-bias & orientation metrics
VCF_DIR="/add/results/no_PoN/annot"
SELECT_DIR="/add/results/no_PoN/select_variants"
OUTPUT_DIR="/add/results/no_PoN/variant_tables"

mkdir -p "$SELECT_DIR" "$OUTPUT_DIR"

for vcf in "$VCF_DIR"/*.vcf.gz; do
  sample=$(basename "$vcf" .func.af.vcf.gz)
  echo "Processing sample: $sample"

  # 1) Select only PASS variants
  gatk SelectVariants \
    -V "$vcf" \
    --exclude-filtered \
    -O "$SELECT_DIR/$sample.PASS.vcf.gz"

  # 2) Dump to TSV with desired metrics:
  # INFO-level: CHROM, POS, REF, ALT, FILTER,
  #             QUAL, MQ, QD,
  #             AF, GNOMAD_AF, FS,
  #             BaseQualityRankSumTest, MappingQualityRankSumTest,
  #             AS_SB_TABLE, STRANDQ, FUNCOTATION
  # FORMAT-level: GT, DP, AD, F1R2, F2R1, SB
  gatk VariantsToTable \
    -V "$SELECT_DIR/$sample.PASS.vcf.gz" \
    -O "$OUTPUT_DIR/$sample.noPoN.tsv" \
    -F CHROM \
    -F POS \
    -F REF \
    -F ALT \
    -F FILTER \
    -F QUAL \
    -F MQ \
    -F QD \
    -F AF \
    -F GNOMAD_AF \
    -F FS \
    -F BaseQualityRankSumTest \
    -F MappingQualityRankSumTest \
    -F AS_SB_TABLE \
    -F STRANDQ \
    -F FUNCOTATION \
    -GF GT \
    -GF DP \
    -GF AD \
    -GF F1R2 \
    -GF F2R1 \
    -GF SB

done
