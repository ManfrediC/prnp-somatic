#!/bin/bash
set -e

# Define directories for input, output, merged files, and the variant list.
# (These directories reflect the output from the force-call workflow for dilution samples.)
inputBaseDir="/add/results/2025-03-15_bwa_force_call_dilutions"
outputBaseDir="/add/results/2025-03-15_bwa_forced_calls_VCF_summary_dilutions"
mergedOutputDir="/add/results/2025-03-15_bwa_forced_calls_merged_VCF_dilutions"
allSummariesDir="${outputBaseDir}/all_summaries"
variantList="/add/variant_VCF/variant_list_short.tsv"

# Ensure the output directories exist.
mkdir -p "$outputBaseDir"
mkdir -p "$mergedOutputDir"
mkdir -p "$allSummariesDir"

# Loop through each variant in the variant list (skipping the header).
tail -n +2 "$variantList" | while IFS=$'\t' read -r variant dbSNP position REF ALT; do
    # Define the input directory for this variant (e.g., /add/results/2025-03-15_bwa_force_call_dilutions/E200K)
    inputDir="${inputBaseDir}/${variant}"
    
    # Check if the input directory exists; if not, skip this variant.
    if [[ ! -d "$inputDir" ]]; then
        echo "Directory $inputDir does not exist. Skipping variant $variant..."
        continue
    fi

    echo "Processing variant: $variant (Position: $position, REF: $REF, ALT: $ALT)"
    
    # Define and create the output directory for this variant.
    variantOutputDir="${outputBaseDir}/${variant}"
    mkdir -p "$variantOutputDir"

    # Compress and index the VCF files for this variant.
    for file in "$inputDir"/*.vcf; do
        if [[ -f "$file" ]]; then
            outFile="${variantOutputDir}/$(basename "$file")"
            bgzip -c "$file" > "${outFile}.gz"
            tabix -p vcf "${outFile}.gz"
        fi
    done

    # Merge the compressed VCF files for this variant into a single merged VCF file.
    mergedVCF="${mergedOutputDir}/merged_${variant}.vcf"
    bcftools merge "${variantOutputDir}"/*.vcf.gz > "$mergedVCF"

    # Extract and format variant information into a summary TSV.
    summaryFile="${allSummariesDir}/${variant}_summary.tsv"
    {
        printf 'Sample\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tFORMAT-DP\tAD-REF\tAD-ALT\n'
        bcftools query --format "[%SAMPLE ]\t%CHROM\t%POS0\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT ]\t[%DP ]\t[%AD ] " "$mergedVCF" |
        awk -F'\t' -v OFS="\t" '
        {
            split($1, samples, " ");
            split($9, genotypes, " ");
            split($10, dps, " ");
            split($11, ads, " ");
            for(i in samples){
                split(ads[i], sampleAd, ",");
                print samples[i], $2, $3, $4, $5, $6, $7, $8, genotypes[i], dps[i], sampleAd[1], sampleAd[2]
            }
        }'
    } | column -s $'\t' -t > "$summaryFile"

    echo "Processed variant $variant. Merged VCF: $mergedVCF. Summary: $summaryFile"
done

echo "Forced-calling summary processing completed for all dilution variants."
