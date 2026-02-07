#!/bin/bash

# this script moves files to the MosaicForecast input container, to prepare for the analysis

# here: bwa samples are used (modify if needed)

set -e

# Define source and destination directories.
srcDir="/add/seq_data/2025-02-07_samples_bwa/"
dstDir="/add/MosaicForecast/input/"

# Create the destination directory if it doesn't exist.
mkdir -p "$dstDir"

# Loop over each BAM file matching the pattern in the source directory.
for bamFile in "$srcDir"/*bwa.picard.markedDup.recal.bam; do
    # Extract the sample name by removing the trailing suffix.
    sample=$(basename "$bamFile" ".bwa.picard.markedDup.recal.bam")
    echo "Copying sample: $sample"
    
    # Copy the BAM file and rename it to sample_bwa.bam.
    cp "$bamFile" "$dstDir/${sample}_bwa.bam"
    
    # The corresponding index file should have the same prefix and a .bai extension.
    baiFile="${bamFile%.bam}.bai"
    if [ -f "$baiFile" ]; then
        cp "$baiFile" "$dstDir/${sample}_bwa.bam.bai"
    else
        echo "Warning: Index file for $bamFile not found!"
    fi
done

echo "All files copied and renamed to $dstDir."
