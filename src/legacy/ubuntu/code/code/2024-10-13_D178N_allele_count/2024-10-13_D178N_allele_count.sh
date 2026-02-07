#!/bin/bash


# here I use the "hack" version of the D178N VCF (E200K VCF modified with the D178N variant) 
# this file contains some formatting (what exactly??) that causes the output VCF to be created as desired



# Define variables
samplesTSV="/add/code/mutect_force_variant/samples.tsv"
fastaFile="/home/mcarta/databases/hg38.fa"
mutationName="D178N"
outputDir="/add/results/2024-10-13_count_D178N"
mutationVCF="/add/variant_VCF/D178N.vcf.gz"

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
