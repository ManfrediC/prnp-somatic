#!/bin/bash
# This script runs samtools flagstat on each BAM file in the dilutions pipeline and writes the results
# to an output directory with filenames based on the sample names.

# Directories
inputdir="/add/seq_data/2025-02-07_dilutions_bwa/"
outputdir="/add/results/pipeline_comparison/2025-02-08_flagstat_dilutions/"

# Create the output directory if it doesn't exist
mkdir -p "$outputdir"

# Define the samples and corresponding BAM files.
samples=("NA100_undil" "NA100_1to10" "NA99A1_undil" "A100_1to2" "NA99A1_1to5" "NA995A05_undil" "NA100_1to2")
bamFiles=(
  "NA100_undil.bwa.picard.markedDup.recal.bam"
  "NA100_1to10.bwa.picard.markedDup.recal.bam"
  "NA99A1_undil.bwa.picard.markedDup.recal.bam"
  "A100_1to2.bwa.picard.markedDup.recal.bam"
  "NA99A1_1to5.bwa.picard.markedDup.recal.bam"
  "NA995A05_undil.bwa.picard.markedDup.recal.bam"
  "NA100_1to2.bwa.picard.markedDup.recal.bam"
)

# Loop over the arrays and run samtools flagstat
for i in "${!samples[@]}"; do
    sampleName="${samples[$i]}"
    bamFile="${bamFiles[$i]}"
    echo "Processing sample: $sampleName"
    samtools flagstat "${inputdir}${bamFile}" > "${outputdir}${sampleName}_bwa.flagstat.txt"
done

echo "Flagstat processing completed for all samples."
