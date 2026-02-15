#!/bin/bash

# Define variables
samplesTSV="/add/code/mutect_force_variant/samples.tsv"
fastaFile="/home/mcarta/databases/hg38.fa"
mutationName="M129V"
outputDir="/add/results/2024-03-04_force_mutect/2024-03-08_M129V/"
#outputDir="/add/results/2024-03-04_force_mutect/2024-03-08_M129V_test/"

mutationVCF="/add/variant_VCF/M129V.vcf.gz"
#mutationVCF="/add/variant_VCF/CJD6_M129V_frequency.vcf.gz"


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