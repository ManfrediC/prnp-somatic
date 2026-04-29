#!/bin/bash

# Base directory containing the subdirectories with VCF files
base_directory="/add/results/mutect_output"

# Path to the reference genome
reference_genome="/home/mcarta/databases/hg38.fa"

# List of directories to process
directories=(
  "2021-02-03_sequencing_of_dilutions"
  "2021-10-19_first_CJD_seq"
  "2023-06-18_CJD_16_samples"
  "2023-06-18_CJD_8_samples"
  "2025-02-08_dilutions_bwa"
  "2025-02-08_samples_bwa"
)

# Loop over each directory
for dir in "${directories[@]}"; do
  # Full path to the current directory
  current_directory="$base_directory/$dir"

  # Loop over each .vcf.gz file in the current directory
  for file in "$current_directory"/*.vcf.gz; do
    # Extract the base name of the file without the extension
    base_name=$(basename "$file" .vcf.gz)

    # Define the output file name
    output_file="${current_directory}/${base_name}_filtered.vcf.gz"

    # Run FilterMutectCalls with the reference genome
    gatk FilterMutectCalls \
      -R "$reference_genome" \
      -V "$file" \
      -O "$output_file"
  done
done
