#!/bin/bash
# This script runs samtools flagstat on each BAM file in the dilutions pipeline and writes the results
# to an output directory with filenames based on the sample names.

# Directories
inputdir="/add/seq_data/2021-02-03_sequencing_of_dilutions/"
outputdir="/add/results/pipeline_comparison/2025-02-08_flagstat_dilutions/"

# Create the output directory if it doesn't exist
mkdir -p "$outputdir"

# Define the samples and corresponding BAM files.
samples=("NA100_undil" "NA100_1to10" "NA99A1_undil" "A100_1to2" "NA99A1_1to5" "NA995A05_undil" "NA100_1to2")
bamFiles=(
  "20210203.0-o23928_1_1-NA100_undil_D01.picard.markedDup.recal.bam"
  "20210203.0-o23928_1_4-NA100_1to10_B02.picard.markedDup.recal.bam"
  "20210203.0-o23928_1_2-NA99A1_undil_D02.picard.markedDup.recal.bam"
  "20210203.0-o23928_1_5-A100_1to2_F01.picard.markedDup.recal.bam"
  "20210203.0-o23928_1_3-NA99A1_1to5_E01.picard.markedDup.recal.bam"
  "20210203.0-o23928_1_7-NA995A05_undil_G01.picard.markedDup.recal.bam"
  "20210203.0-o23928_1_5-A100_1to2_F01.picard.markedDup.recal.bam"
)

# Loop over the arrays and run samtools flagstat
for i in "${!samples[@]}"; do
    sampleName="${samples[$i]}"
    bamFile="${bamFiles[$i]}"
    echo "Processing sample: $sampleName"
    samtools flagstat "${inputdir}${bamFile}" > "${outputdir}${sampleName}_orig.flagstat.txt"
done

echo "Flagstat processing completed for all samples."
