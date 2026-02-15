#!/bin/bash
# date: 2025-04-12
# This script runs GATK Funcotator on normalized VCF files from the dilutions dataset.
# It uses a custom reference (chr2, chr4, chr20) and the re-organized Funcotator data bundle
# located at /add/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/ (which now contains
# a subdirectory named "hg38" with all the data sources).
# Annotated VCFs are saved in a dedicated output directory.

# Set paths for the reference and data sources:
dbDir="/home/mcarta/databases"
fastaFile="$dbDir/chr2_chr4_chr20.fasta"

# The data sources path now points to the re-organized bundle.
funcotatorData="/add/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38"

##########################################
# Input and Output Directories for Dilutions
##########################################
vcf_dilutions="/add/results/annotation_output/1_normalisation/dilutions/"
out_dilutions="/add/results/annotation_output/2_funcotator/dilutions/"

# Create the output directory if it doesn't exist
mkdir -p "$out_dilutions"

##########################################
# Run Funcotator on each normalized VCF file
##########################################
for vcf in "$vcf_dilutions"*.vcf.gz; do
    base=$(basename "$vcf" .vcf.gz)
    echo "Annotating $vcf with Funcotator..."
    
    # Define output filename with ".funcotated.vcf.gz" suffix
    output_file="${out_dilutions}${base}.funcotated.vcf.gz"
    
    gatk Funcotator \
        --variant "$vcf" \
        --reference "$fastaFile" \
        --ref-version hg38 \
        --data-sources-path "$funcotatorData" \
        --output "$output_file" \
        --output-file-format VCF

    echo "Annotation complete for: $output_file"
done

echo "Funcotator annotation completed for all dilutions files."
