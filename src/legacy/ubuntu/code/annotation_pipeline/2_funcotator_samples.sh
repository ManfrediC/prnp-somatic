#!/bin/bash
# date: 2025-04-12
# This script runs GATK Funcotator on normalized VCF files from the samples dataset.
# It uses a custom reference (chr2, chr4, chr20) and the re-organized Funcotator data bundle
# located at /add/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38.
# Annotated VCF files are written to the designated output directory.

# Set up the reference and data sources:
dbDir="/home/mcarta/databases"
fastaFile="$dbDir/chr2_chr4_chr20.fasta"
funcotatorData="/add/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38"

# Define input and output directories for the samples:
vcf_samples="/add/results/annotation_output/1_normalisation/samples/"
out_samples="/add/results/annotation_output/2_funcotator/samples/"

# Create the output directory if it does not already exist
mkdir -p "$out_samples"

# Process each normalized VCF file in the samples directory
for vcf in "$vcf_samples"*.vcf.gz; do
    base=$(basename "$vcf" .vcf.gz)
    echo "Annotating $vcf with Funcotator..."
    
    # Define the output filename, appending ".funcotated.vcf.gz" to the original base name
    output_file="${out_samples}${base}.funcotated.vcf.gz"
    
    gatk Funcotator \
        --variant "$vcf" \
        --reference "$fastaFile" \
        --ref-version hg38 \
        --data-sources-path "$funcotatorData" \
        --output "$output_file" \
        --output-file-format VCF
    
    echo "Annotation complete for: $output_file"
done

echo "Funcotator annotation completed for all samples files."
