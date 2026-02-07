#!/bin/bash
# date: 2025-04-12
# Pre-Processing Checks:
# This script verifies that normalized VCFs are properly formatted, compressed,
# and indexed; then runs quality metric checks using bcftools stats.

##########################################
# Directories containing the normalized VCF files
##########################################
vcf_samples="/add/results/annotation_output/1_normalisation/samples/"
vcf_dilutions="/add/results/annotation_output/1_normalisation/dilutions/"

# Set up directories to store quality metrics reports
stats_samples="/add/results/annotation_output/1_normalisation/samples_stats/"
stats_dilutions="/add/results/annotation_output/1_normalisation/dilutions_stats/"

# Create stats output directories if they don't exist
mkdir -p "$stats_samples"
mkdir -p "$stats_dilutions"

##########################################
# Pre-Processing Checks Loop
##########################################

for type in samples dilutions; do
    if [ "$type" == "samples" ]; then
        vcf_dir="$vcf_samples"
        stats_dir="$stats_samples"
    else
        vcf_dir="$vcf_dilutions"
        stats_dir="$stats_dilutions"
    fi

    # Loop over each normalized VCF file in the specified directory
    for vcf in "$vcf_dir"*.vcf.gz; do
        base=$(basename "$vcf" .vcf.gz)
        echo "Processing $vcf for quality checks ..."

        # Check if the VCF index (.tbi) exists; if not, create it
        if [ ! -f "${vcf}.tbi" ]; then
            echo "Indexing $vcf ..."
            bcftools index "$vcf"
        fi

        # Optionally, inspect the header for formatting correctness
        echo "Inspecting header for $vcf ..."
        bcftools view -h "$vcf" | head -n 10

        # Generate basic quality statistics with bcftools stats
        stats_file="${stats_dir}${base}.stats"
        echo "Generating statistics for $vcf ..."
        bcftools stats "$vcf" > "$stats_file"

        echo "Completed processing for $vcf; stats saved to $stats_file"
    done
done

echo "Pre-processing checks completed."
