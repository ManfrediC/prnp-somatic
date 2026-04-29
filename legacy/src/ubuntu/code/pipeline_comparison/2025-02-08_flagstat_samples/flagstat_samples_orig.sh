#!/bin/bash
# This script runs samtools flagstat on original pipeline BAM files from three input directories:
#   1. /add/seq_data/2023-06-18_CJD_8_samples/
#   2. /add/seq_data/2023-06-18_CJD_16_samples/
#   3. /add/seq_data/2021-10-19_first_CJD_seq/
#
# Only files ending in "picard.markedDup.recal.bam" are processed.
# The sample name is extracted using a regex that captures the text between the dash (-)
# and the first underscore (_) in the filename.
#
# Output summaries are written to:
#   /add/code/pipeline_comparison/2025-02-08_flagstat_samples/

# Set output directory and create it if needed
outputdir="/add/results/pipeline_comparison/2025-02-08_flagstat_samples/"
mkdir -p "$outputdir"

# Define the input directories
inputDirs=(
  "/add/seq_data/2023-06-18_CJD_8_samples/"
  "/add/seq_data/2023-06-18_CJD_16_samples/"
  "/add/seq_data/2021-10-19_first_CJD_seq/"
)

# Loop over each input directory
for inputdir in "${inputDirs[@]}"; do
    # Process only files ending with "picard.markedDup.recal.bam"
    for bamFile in "$inputdir"*picard.markedDup.recal.bam; do
        # Get the base filename (without the directory path)
        base=$(basename "$bamFile")
        # Extract the sample name using sed.
        # This regex removes everything up to the dash and then captures the text until the first underscore.
        sample=$(echo "$base" | sed -E 's/.*-([^_]+)_.*/\1/')
        echo "Processing sample: $sample from file: $bamFile"
        samtools flagstat "$bamFile" > "${outputdir}${sample}_orig.flagstat.txt"
    done
done

echo "Flagstat processing completed for all original pipeline samples."
