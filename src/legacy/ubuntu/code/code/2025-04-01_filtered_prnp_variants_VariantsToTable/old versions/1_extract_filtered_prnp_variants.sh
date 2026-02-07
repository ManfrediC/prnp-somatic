#!/bin/bash
# Define directories

### For dilutions (disabled)
# INPUT_DIR="/add/results/mutect_output/2025-02-08_dilutions_bwa"
# OUTPUT_DIR="/add/results/2025-04-01_filtered_prnp_vcf/dilutions"

### For samples (enabled)
INPUT_DIR="/add/results/mutect_filtered_collection/samples"
OUTPUT_DIR="/add/results/2025-04-01_filtered_prnp_vcf/samples"

TMP_DIR="/tmp"  # Change this if you prefer another temporary location

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all *_filtered.vcf.gz files in the input directory
for file in "$INPUT_DIR"/*_filtered.vcf.gz; do
    # Extract the sample name by removing the directory and the suffix
    sample=$(basename "$file" | sed 's/\.Mutect_filtered\.vcf\.gz//')
    echo "Processing sample: $sample"
    
    # Define a temporary file name for the unzipped VCF
    temp_vcf="${TMP_DIR}/${sample}_temp.vcf"

    # Unzip the file to a temporary file
    gunzip -c "$file" > "$temp_vcf"

    # Filter the VCF:
    # - Always include header lines (starting with "#") to preserve correct VCF formatting.
    # - Include only variant lines on "chr20" with positions between 4686134 and 4701605 (PRNP protein coding) that have FILTER field equal to "PASS".
    awk -v OFS="\t" '{
        if($0 ~ /^#/) {
            print $0
        } else if($1=="chr20" && $2>=4686134 && $2<=4701605 && $7=="PASS") {
            print $0
        }
    }' "$temp_vcf" > "$OUTPUT_DIR/${sample}_filtered_prnp.vcf"

    echo "Filtered VCF saved as: $OUTPUT_DIR/${sample}_filtered_prnp.vcf"

    # Automatically delete the temporary file to save disk space
    rm "$temp_vcf"
    echo "Deleted temporary file: ${temp_vcf}"
done
