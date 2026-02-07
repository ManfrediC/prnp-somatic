#!/bin/bash
# date: 2025-04-12
# This pipeline uses GATK Funcotator to functionally annotate normalized VCF files.
# The normalized files have been prepared by decomposing multi-allelic sites and left-aligning indels.

# Databases and Reference
dbDir="/home/mcarta/databases"
fastaFile="$dbDir/chr2_chr4_chr20.fasta"
funcotatorData="/home/mcarta/funcotator_data"  # Path to Funcotator data bundle

##########################################
# Input Directories (normalized VCFs)
##########################################
vcf_samples="/add/results/annotation_output/1_normalisation/samples/"
vcf_dilutions="/add/results/annotation_output/1_normalisation/dilutions/"

##########################################
# Output Directories for Funcotator annotated VCF files
##########################################
out_samples="/add/results/annotation_output/2_funcotator/samples/"
out_dilutions="/add/results/annotation_output/2_funcotator/dilutions/"

# Create output directories if they do not exist
mkdir -p "$out_samples"
mkdir -p "$out_dilutions"

##########################################
# Funcotator Annotation with GATK
##########################################
for type in samples dilutions; do
    if [ "$type" == "samples" ]; then
        input_dir="$vcf_samples"
        output_dir="$out_samples"
    else
        input_dir="$vcf_dilutions"
        output_dir="$out_dilutions"
    fi

    # Loop through each normalized VCF file
    for vcf in "$input_dir"*.vcf.gz; do
        base=$(basename "$vcf" .vcf.gz)
        echo "Annotating $vcf with Funcotator ..."

        # Define the output filename; appending .funcotated.vcf.gz to the basename
        output_file="${output_dir}${base}.funcotated.vcf.gz"

        gatk Funcotator \
            -R "$fastaFile" \
            -V "$vcf" \
            -O "$output_file" \
            --data-sources-path "$funcotatorData" \
            --ref-version custom \
            --output-file-format VCF

        echo "Annotation complete: $output_file"
    done
done

echo "Funcotator annotation completed for all files."
