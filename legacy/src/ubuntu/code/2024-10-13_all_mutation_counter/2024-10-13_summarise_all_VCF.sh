#!/bin/bash

# Define the directories for input, output, and merged files
inputBaseDir="/add/results/2024-10-13_countAll_test"
outputBaseDir="/add/results/2024-10-13_VCF_summary"
mergedOutputDir="/add/results/2024-10-13_merged_VCF"
allSummariesDir="${outputBaseDir}/all_summaries"
variantList="/add/variant_VCF/variant_list_short.tsv"

# Ensure the output directories exist
mkdir -p "$outputBaseDir"
mkdir -p "$mergedOutputDir"
mkdir -p "$allSummariesDir"

# Read each variant from the variant list (skipping the header)
tail -n +2 "$variantList" | while IFS=$'\t' read -r variant dbSNP position REF ALT; do
    # Define input directory for this variant
    inputDir="${inputBaseDir}/${variant}"
    
    # Check if the input directory exists
    if [[ ! -d "$inputDir" ]]; then
        echo "Directory $inputDir does not exist. Skipping..."
        continue
    fi

    # Define the output directory for this variant
    outputDir="${outputBaseDir}/${variant}"
    mkdir -p "$outputDir"

    # Compress and index the VCF files for this variant
    for file in "${inputDir}"/*.vcf; do
        if [[ -f "$file" ]]; then  # Check if file exists
            bgzip -c "$file" > "${outputDir}/${file##*/}.gz"
            tabix -p vcf "${outputDir}/${file##*/}.gz"
        fi
    done

    # Merge the compressed VCF files for this variant into a single merged VCF file in the merged directory
    mergedVCF="${mergedOutputDir}/merged_${variant}.vcf"
    bcftools merge "${outputDir}"/*.vcf.gz > "$mergedVCF"

    # Extract and format variant information, saving it to a summary file in the all_summaries directory
    summaryFile="${allSummariesDir}/${variant}_summary.tsv"
    {
        printf 'Sample\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tFORMAT-DP\tAD-REF\tAD-ALT\n';
        bcftools query --format "[%SAMPLE ]\t%CHROM\t%POS0\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT ]\t[%DP ]\t[%AD ] " "$mergedVCF" |
        awk -F'\t' -v OFS="\t" '
        { 
            split($1, samples, " "); 
            split($9, genotypes, " "); 
            split($10, dps, " ");
            split($11, ads, " "); 
            for(s in samples){ 
               split(ads[s], sampleAd, ",");  
               print samples[s], $2, $3, $4, $5, $6, $7, $8, genotypes[s], dps[s], sampleAd[1], sampleAd[2]
            }
        }';
    } | column -s$'\t' -t > "$summaryFile"

    echo "Processed and summarized ${variant} VCF files. Merged file saved as ${mergedVCF}. Summary file saved as ${summaryFile}"
done
