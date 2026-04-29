#!/bin/bash

# Directory containing the VCF files
directory="/add/results/mutect_output/2025-02-08_dilutions_bwa"

# Path to the reference genome
reference_genome="/home/mcarta/databases/hg38.fa"

# Loop over each .vcf.gz file in the directory
for file in "$directory"/*.vcf.gz; do
  # Extract the base name of the file without the extension
  base_name=$(basename "$file" .vcf.gz)

  # Define the output file name
  output_file="${directory}/${base_name}_filtered.vcf.gz"

  # Run FilterMutectCalls with the reference genome
  gatk FilterMutectCalls \
    -R "$reference_genome" \
    -V "$file" \
    -O "$output_file"
done
