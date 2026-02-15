#!/bin/bash
# move_original_files_for_MF.sh
#
# This script copies sample files from the specified source directories into the
# MosaicForecast input directory (/add/MosaicForecast/input), but only those files
# that end in "picard.markedDup.recal.bam" or "picard.markedDup.recal.bai".
#
# When copying, the file is renamed so that the sample identifier is extracted from
# the original filename and the new name is:
#   <SampleID>_original.bam   OR   <SampleID>_original.bai
#
# Example:
#   "20211027.0-o26424_1_1-CJD1_fr.picard.markedDup.recal.bam" will be copied as "CJD1_original.bam"
#   "312092_11-CJD28_F01_S6_R1_001.fastq.gz.picard.markedDup.recal.bai" will be copied as "CJD28_original.bai"
#
# Adjust the source directories as needed.

set -euo pipefail

# Enable nullglob so that if no files match a glob, the loop will simply skip.
shopt -s nullglob

# Destination directory for MosaicForecast input files
DEST_DIR="/add/MosaicForecast/input"
mkdir -p "$DEST_DIR"

# Define the source directories where sample files are located
SRC_DIRS=(
    "/add/seq_data/2021-10-19_first_CJD_seq"
    "/add/seq_data/2023-06-18_CJD_16_samples"
    "/add/seq_data/2023-06-18_CJD_8_samples"
)

echo "===== Processing BAM Files (picard.markedDup.recal.bam) ====="
for src_dir in "${SRC_DIRS[@]}"; do
    echo "Processing source directory: $src_dir"
    # Loop only over files that end with picard.markedDup.recal.bam
    for file in "$src_dir"/*picard.markedDup.recal.bam; do
        # If no file matches the glob, continue.
        if [ ! -e "$file" ]; then
            continue
        fi

        filename=$(basename "$file")
        # Extract sample identifier: looks for either "CJD" or "Ctrl" followed by digits.
        sample=$(echo "$filename" | grep -oE '(CJD|Ctrl)[0-9]+')
        if [ -z "$sample" ]; then
            echo "Warning: Could not extract sample ID from '$filename'. Skipping."
            continue
        fi

        new_filename="${sample}_original.bam"
        echo "Copying $filename to $new_filename"
        cp "$file" "$DEST_DIR/$new_filename"
    done
done

echo "===== Processing BAI Files (picard.markedDup.recal.bai) ====="
for src_dir in "${SRC_DIRS[@]}"; do
    echo "Processing source directory: $src_dir"
    # Loop only over files that end with picard.markedDup.recal.bai
    for file in "$src_dir"/*picard.markedDup.recal.bai; do
        if [ ! -e "$file" ]; then
            continue
        fi

        filename=$(basename "$file")
        sample=$(echo "$filename" | grep -oE '(CJD|Ctrl)[0-9]+')
        if [ -z "$sample" ]; then
            echo "Warning: Could not extract sample ID from '$filename'. Skipping."
            continue
        fi

        new_filename="${sample}_original.bai"
        echo "Copying $filename to $new_filename"
        cp "$file" "$DEST_DIR/$new_filename"
    done
done

echo "All files have been copied and renamed successfully."
