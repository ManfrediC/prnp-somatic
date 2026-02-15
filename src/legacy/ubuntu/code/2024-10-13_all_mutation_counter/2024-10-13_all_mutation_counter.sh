#!/bin/bash

# Define input files
samplesTSV="/add/code/mutect_force_variant/samples.tsv"
fastaFile="/home/mcarta/databases/hg38.fa"
variantList="/add/variant_VCF/variant_list_short.tsv"

# Loop through each variant in the variant list
while IFS=$'\t' read -r variant dbSNP position REF ALT; do
    # Skip header line if present
    if [[ "$variant" == "variant" && "$dbSNP" == "dbSNP id" ]]; then
        continue
    fi

    # Define output directory and mutation-specific VCF file
    outputDir="/add/results/2024-10-13_countAll/${variant}"
    mutationVCF="/add/variant_VCF/${variant}.vcf.gz"

    # Create the output directory if it doesn't exist
    mkdir -p "$outputDir"

    # Loop: read sampleName and inputFile from TSV for each mutation
    while IFS=$'\t' read -r sampleName inputFile; do
        # Skip header line if present
        if [[ "$sampleName" == "sampleName" && "$inputFile" == "inputFile" ]]; then
            continue
        fi

        # Perform Mutect2 analysis for each sample and mutation
        gatk Mutect2 \
            -R "$fastaFile" \
            --alleles "$mutationVCF" \
            -L "$mutationVCF" \
            -I "$inputFile" \
            -O "${outputDir}/${sampleName}_${variant}_frequency.vcf"
    done < "$samplesTSV"

done < "$variantList"
