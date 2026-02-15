#!/bin/bash
#
# This script copies dilution BAM and BAI files from the source directory:
#   /add/seq_data/2021-02-03_sequencing_of_dilutions/
# It processes only files ending in "picard.markedDup.recal.bam" or "picard.markedDup.recal.bai".
#
# An array of sample names is defined:
#    A100_1to2, NA99A1_1to5, NA100_1to10, NA995A05_undil, NA99A1_undil, NA100_undil
#
# For each sample, the script searches for matching files in the source directory.
# If a matching file is found, it is copied to the destination directory
# (/add/MosaicForecast/input) and renamed as:
#    <Sample>_original.bam   or   <Sample>_original.bai
#
# Adjust the source and destination directories if needed.

set -euo pipefail
shopt -s nullglob

# Define source and destination directories.
SRC_DIR="/add/seq_data/2021-02-03_sequencing_of_dilutions"
DEST_DIR="/add/MosaicForecast/input"
mkdir -p "$DEST_DIR"

# Define the list of sample names.
samples=("A100_1to2" "NA99A1_1to5" "NA100_1to10" "NA995A05_undil" "NA99A1_undil" "NA100_undil")

# Loop over each sample name.
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"
    
    # Search for the BAM file for this sample.
    bam_files=("$SRC_DIR"/*"${sample}"*.picard.markedDup.recal.bam)
    if [ ${#bam_files[@]} -eq 0 ]; then
        echo "  No BAM file found for sample $sample"
    else
        if [ ${#bam_files[@]} -gt 1 ]; then
            echo "  Multiple BAM files found for sample $sample, using the first one: ${bam_files[0]}"
        fi
        cp "${bam_files[0]}" "$DEST_DIR/${sample}_original.bam"
        echo "  Copied $(basename "${bam_files[0]}") as ${sample}_original.bam"
    fi

    # Search for the BAI file for this sample.
    bai_files=("$SRC_DIR"/*"${sample}"*.picard.markedDup.recal.bai)
    if [ ${#bai_files[@]} -eq 0 ]; then
        echo "  No BAI file found for sample $sample"
    else
        if [ ${#bai_files[@]} -gt 1 ]; then
            echo "  Multiple BAI files found for sample $sample, using the first one: ${bai_files[0]}"
        fi
        cp "${bai_files[0]}" "$DEST_DIR/${sample}_original.bai"
        echo "  Copied $(basename "${bai_files[0]}") as ${sample}_original.bai"
    fi

    echo ""
done

echo "All required files have been copied and renamed successfully."
