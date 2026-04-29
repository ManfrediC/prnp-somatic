#!/bin/bash

# this script takes the samples listed in samples.tsv and looks for the variant specified by $mutationVCF (e.g. E200K)

# Define variables
samplesTSV="/add/results/2024-03-04_force_mutect/samples.tsv"
fastaFile="/home/mcarta/databases/hg38.fa"
mutationName="D178N"
outputDir="/add/results/2024-03-04_force_mutect/2024-03-06_D178N/"
mutationVCF="/add/results/2024-03-05_ref_VCFs_for_mutations/D178N.vcf.gz"

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
