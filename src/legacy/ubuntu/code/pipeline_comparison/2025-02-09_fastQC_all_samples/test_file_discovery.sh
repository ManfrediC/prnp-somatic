#!/bin/bash

# this script tests whether the files we need for the FastQC loop can be found

set -e

# Define the input directories.
inputDirs=(
  "/add/seq_data/2021-10-19_first_CJD_seq/"
  "/add/seq_data/2023-06-18_CJD_16_samples/"
  "/add/seq_data/2023-06-18_CJD_8_samples/"
)

echo "Testing file discovery and sample name extraction..."
echo "-----------------------------------------"

# Loop over each input directory.
for dir in "${inputDirs[@]}"; do
    echo "Directory: $dir"
    # Loop over all files matching the pattern *_R1_001.fastq.gz in the directory.
    for read1 in $dir*_R1_001.fastq.gz; do
        # If no file is found, skip to next directory.
        [ -e "$read1" ] || continue

        # Get the base filename.
        base=$(basename "$read1")
        # Use sed to extract the sample name.
        # For example, "312091_5-Ctrl1_G01_S6_R1_001.fastq.gz" yields "Ctrl1".
        sample=$(echo "$base" | sed -E 's/.*-([^_]+)_.*/\1/')
        # Determine the corresponding R2 file by replacing _R1_001.fastq.gz with _R2_001.fastq.gz.
        read2="${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

        # Print the discovered information.
        echo "Sample: $sample"
        echo "R1: $read1"
        echo "R2: $read2"
        echo "----------------"
    done
done

echo "File discovery test completed."
