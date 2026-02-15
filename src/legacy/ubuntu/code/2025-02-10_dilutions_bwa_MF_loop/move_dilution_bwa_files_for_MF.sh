#!/bin/bash
#
# This script copies dilution BAM and BAI files from the source directory:
#   /add/seq_data/2025-02-07_dilutions_bwa/
# but only those files that end in "picard.markedDup.recal.bam" or "picard.markedDup.recal.bai".
#
# It extracts the sample name from each file by removing the suffix ".bwa.picard.markedDup.recal.bam"
# (or ".bwa.picard.markedDup.recal.bai") and then renames the file to:
#   <Sample>_bwa.bam   or   <Sample>_bwa.bai
#
# For example:
#   "NA995A05_undil.bwa.picard.markedDup.recal.bam" becomes "NA995A05_undil_bwa.bam"
#   "NA995A05_undil.bwa.picard.markedDup.recal.bai" becomes "NA995A05_undil_bwa.bai"
#
# Adjust the source and destination directories if needed.

set -euo pipefail
shopt -s nullglob

# Source directory containing the bwa files
SRC_DIR="/add/seq_data/2025-02-07_dilutions_bwa"

# Destination directory for MosaicForecast input files
DEST_DIR="/add/MosaicForecast/input"
mkdir -p "$DEST_DIR"

echo "===== Processing BAM Files ====="
for file in "$SRC_DIR"/*picard.markedDup.recal.bam; do
    # Skip if no matching file is found.
    [ -e "$file" ] || continue
    filename=$(basename "$file")
    # Remove the suffix ".bwa.picard.markedDup.recal.bam" to extract the sample name.
    sample=$(echo "$filename" | sed 's/\.bwa\.picard\.markedDup\.recal\.bam$//')
    new_filename="${sample}_bwa.bam"
    echo "Copying $filename to $new_filename"
    cp "$file" "$DEST_DIR/$new_filename"
done

echo "===== Processing BAI Files ====="
for file in "$SRC_DIR"/*picard.markedDup.recal.bai; do
    [ -e "$file" ] || continue
    filename=$(basename "$file")
    sample=$(echo "$filename" | sed 's/\.bwa\.picard\.markedDup\.recal\.bai$//')
    new_filename="${sample}_bwa.bai"
    echo "Copying $filename to $new_filename"
    cp "$file" "$DEST_DIR/$new_filename"
done

echo "All required files have been copied and renamed successfully."
