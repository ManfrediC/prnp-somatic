#!/bin/bash

# Directory with filtered VCF files

# on 2025-04-19: PoN pipeline results
VCF_DIR="/add/results/annotation_output/funcotator_PoN/"
OUTPUT_DIR="/add/results/variant_tables/samples/"

mkdir -p "$OUTPUT_DIR"

# Loop over each VCF file
for vcf in "$VCF_DIR"/*.vcf.gz; do
  sample=$(basename "$vcf" | sed 's/\.funcotated\.vcf\.gz$//')
  echo "Processing sample: $sample from file: $vcf"
  
  gatk VariantsToTable \
    --variant "$vcf" \
    --output "${OUTPUT_DIR}/${sample}.tsv"

done



# old version

#VCF_DIR="/add/results/2025-04-01_filtered_prnp_vcf/samples"
#OUTPUT_DIR="/add/results/2025-04-05_combined_variant_tables/samples"

#mkdir -p "$OUTPUT_DIR"

# Loop over each VCF file
#for vcf in "$VCF_DIR"/*.vcf; do
#  sample=$(basename "$vcf" | sed 's/_filtered_prnp.vcf//')
#  echo "Processing sample: $sample from file: $vcf"
  
#  gatk VariantsToTable \
#    --variant "$vcf" \
#    --output "${OUTPUT_DIR}/${sample}.tsv"

#done
