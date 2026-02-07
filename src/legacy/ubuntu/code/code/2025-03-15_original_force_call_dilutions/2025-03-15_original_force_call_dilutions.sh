#!/bin/bash
set -e

# Define essential inputs.
fastaFile="/home/mcarta/databases/hg38.fa"
variantList="/add/variant_VCF/variant_list_short.tsv"

# Input directory with original pipeline dilution BAM files.
inputDir="/add/seq_data/2021-02-03_sequencing_of_dilutions/"

# Base output directory for forced calling results.
outputBase="/add/results/2025-03-15_original_force_call_dilutions/"

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
    
    # Loop over each BAM file matching the expected suffix in the input directory.
    for bamFile in "$inputDir"/*picard.markedDup.recal.bam; do
        # Extract the base filename (removing the trailing .picard.markedDup.recal.bam)
        base=$(basename "$bamFile" ".picard.markedDup.recal.bam")
        # The base filename is like: 20210203.0-o23928_1_1-NA100_undil_D01
        # Extract the sample name: take the part after the last dash, then use awk to combine the first two underscore fields.
        sample_full=$(echo "$base" | awk -F '-' '{print $NF}')  # e.g., NA100_undil_D01
        sample=$(echo "$sample_full" | awk -F '_' '{print $1"_"$2}')  # yields "NA100_undil"
        echo "   Processing sample: $sample"
        
        # Run forced calling with Mutect2 for this variant on the sample's BAM.
        gatk Mutect2 \
            -R "$fastaFile" \
            --alleles "$mutationVCF" \
            -L "$mutationVCF" \
            -I "$bamFile" \
            -O "${outputDir}/${sample}_orig_${variant}_frequency.vcf"
    done

done < "$variantList"

echo "Forced Mutect2 processing completed for all dilution variants and samples."
