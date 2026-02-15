#!/bin/bash

# this script is used to run the force-call Mutect2 workflow on the dilution samples (output of Bwa pipeline)

set -e

# Define essential inputs.
fastaFile="/home/mcarta/databases/hg38.fa"
variantList="/add/variant_VCF/variant_list_short.tsv"

# Input directory with dilution samples (aggressive pipeline).
inputDir="/add/seq_data/2025-02-07_dilutions_bwa/"

# Base output directory for forced calling results.
outputBase="/add/results/2025-03-15_bwa_force_call_dilutions/"

# Loop through each variant in the variant list.
while IFS=$'\t' read -r variant dbSNP position REF ALT; do
    # Skip header line if present.
    if [[ "$variant" == "variant" && "$dbSNP" == "dbSNP id" ]]; then
        continue
    fi

    # Define output directory for this variant and the mutation-specific VCF file.
    outputDir="${outputBase}/${variant}"
    mutationVCF="/add/variant_VCF/${variant}.vcf.gz"

    # Create the output directory if it doesn't exist.
    mkdir -p "$outputDir"

    echo "------------------------------------------------------------"
    echo "Processing variant: $variant (Position: $position, REF: $REF, ALT: $ALT)"
    
    # Loop over each BAM file matching the expected suffix in the dilution input directory.
    for bamFile in "$inputDir"/*bwa.picard.markedDup.recal.bam; do
        # Extract the sample name: this will remove the directory and the trailing ".bwa.picard.markedDup.recal.bam".
        sampleName=$(basename "$bamFile" ".bwa.picard.markedDup.recal.bam")
        echo "   Processing sample: $sampleName"
        
        # Run forced calling with Mutect2 for this variant on the sample's BAM file.
        gatk Mutect2 \
            -R "$fastaFile" \
            --alleles "$mutationVCF" \
            -L "$mutationVCF" \
            -I "$bamFile" \
            -O "${outputDir}/${sampleName}_bwa_${variant}_frequency.vcf"
    done

done < "$variantList"

echo "Forced Mutect2 processing completed for all dilution variants and samples."
