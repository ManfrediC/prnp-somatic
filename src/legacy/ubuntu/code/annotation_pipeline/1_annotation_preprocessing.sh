#!/bin/bash
# date: 2025-04-12
# This pipeline serves to prepare annotation of variants found with the BWA ("aggressive")
# preprocessing pipeline, followed by Mutect2 and FilterMutectCalls

# Databases
dbDir="/home/mcarta/databases"
fastaFile="$dbDir/chr2_chr4_chr20.fasta"

##########################################
# Input Directories for filtered VCF files
##########################################
input_samples="/add/results/mutect_filtered_collection/samples/"
input_dilutions="/add/results/mutect_filtered_collection/dilutions/"

##########################################
# Output Directories for normalized VCF files
##########################################
out_samples="/add/results/annotation_output/1_normalisation/samples/"
out_dilutions="/add/results/annotation_output/1_normalisation/dilutions/"

# Create output directories if they don't exist
mkdir -p "$out_samples"
mkdir -p "$out_dilutions"

##########################################
# Normalisation with bcftools
##########################################

# Loop over both sample types
for type in samples dilutions; do
    # Set directories based on type
    if [ "$type" == "samples" ]; then
        input_dir="$input_samples"
        output_dir="$out_samples"
    else
        input_dir="$input_dilutions"
        output_dir="$out_dilutions"
    fi

    # Loop through each filtered VCF file in the current directory
    for vcf in "$input_dir"*filtered*.vcf.gz; do
        # Obtain the base filename without extensions for output naming
        base=$(basename "$vcf" .vcf.gz)
        
        echo "Processing $vcf ..."

        # Step 1: Decompose Multi-allelic Sites with bcftools norm -m -any
        decomposed_vcf="${output_dir}${base}.decomposed.vcf.gz"
        bcftools norm -m -any -Oz -o "$decomposed_vcf" "$vcf"

        # Step 2: Left-align Indels using the reference genome
        normalized_vcf="${output_dir}${base}.normalized.vcf.gz"
        bcftools norm -f "$fastaFile" -Oz -o "$normalized_vcf" "$decomposed_vcf"

        # Index the normalized VCF file for fast access in downstream analyses
        bcftools index "$normalized_vcf"
        
        echo "Finished processing $vcf"
    done
done
