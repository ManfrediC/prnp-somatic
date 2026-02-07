#!/bin/bash

# Define variables
samplesTSV="/add/results/2024-03-04_force_mutect/samples.tsv"
fastaFile="/home/mcarta/databases/hg38.fa"
mutationName="E200K"
outputDir="/add/results/2024-03-04_force_mutect/VCFoutput/"
mutationVCF="/add/variant_VCF/E200K.vcf.gz"

# Loop: read sampleName and inputFile from TSV
while IFS=$'\t' read -r sampleName inputFile; do
    # Skip header line if present
    if [[ "$sampleName" == "sampleName" && "$inputFile" == "inputFile" ]]; then
        continue
    fi

    # Perform Mutect2 analysis for each sample
    gatk Mutect2 \
        -R "$fastaFile" \
        --alleles "$mutationVCF" \
        -L "$mutationVCF" \
        -I "$inputFile" \
        -O "${outputDir}/${sampleName}_${mutationName}_frequency.vcf"
done < "$samplesTSV"
