#!/bin/bash
# flagstat_samples_bwa.sh
# This script runs samtools flagstat on each BAM file in the samples_bwa directory and
# writes the output to an output directory, with filenames that include the sample name.
#
# Input directory and output directory:
inputdir="/add/seq_data/2025-02-07_samples_bwa/"
outputdir="/add/results/pipeline_comparison/2025-02-08_flagstat_samples/"

# Create the output directory if it doesn't exist
mkdir -p "$outputdir"

# Define arrays for sample names and their corresponding BAM filenames.
samples=("CJD31" "CJD22" "CJD8" "CJD12" "CJD13" "CJD7" "CJD4" "CJD26" "CJD19" "CJD2" "CJD24" "CJD18" "Ctrl2" "Ctrl1" "CJD14" "CJD28" "CJD30" "Ctrl5" "Ctrl7" "CJD15" "Ctrl3" "CJD21" "CJD1" "CJD3" "CJD23" "CJD6" "CJD27" "CJD11" "CJD29" "CJD25" "Ctrl4" "CJD5")

bamFiles=(
  "CJD31.bwa.picard.markedDup.recal.bam"
  "CJD22.bwa.picard.markedDup.recal.bam"
  "CJD8.bwa.picard.markedDup.recal.bam"
  "CJD12.bwa.picard.markedDup.recal.bam"
  "CJD13.bwa.picard.markedDup.recal.bam"
  "CJD7.bwa.picard.markedDup.recal.bam"
  "CJD4.bwa.picard.markedDup.recal.bam"
  "CJD26.bwa.picard.markedDup.recal.bam"
  "CJD19.bwa.picard.markedDup.recal.bam"
  "CJD2.bwa.picard.markedDup.recal.bam"
  "CJD24.bwa.picard.markedDup.recal.bam"
  "CJD18.bwa.picard.markedDup.recal.bam"
  "Ctrl2.bwa.picard.markedDup.recal.bam"
  "Ctrl1.bwa.picard.markedDup.recal.bam"
  "CJD14.bwa.picard.markedDup.recal.bam"
  "CJD28.bwa.picard.markedDup.recal.bam"
  "CJD30.bwa.picard.markedDup.recal.bam"
  "Ctrl5.bwa.picard.markedDup.recal.bam"
  "Ctrl7.bwa.picard.markedDup.recal.bam"
  "CJD15.bwa.picard.markedDup.recal.bam"
  "Ctrl3.bwa.picard.markedDup.recal.bam"
  "CJD21.bwa.picard.markedDup.recal.bam"
  "CJD1.bwa.picard.markedDup.recal.bam"
  "CJD3.bwa.picard.markedDup.recal.bam"
  "CJD23.bwa.picard.markedDup.recal.bam"
  "CJD6.bwa.picard.markedDup.recal.bam"
  "CJD27.bwa.picard.markedDup.recal.bam"
  "CJD11.bwa.picard.markedDup.recal.bam"
  "CJD29.bwa.picard.markedDup.recal.bam"
  "CJD25.bwa.picard.markedDup.recal.bam"
  "Ctrl4.bwa.picard.markedDup.recal.bam"
  "CJD5.bwa.picard.markedDup.recal.bam"
)

# Loop over each sample and run samtools flagstat
for i in "${!samples[@]}"; do
    sampleName="${samples[$i]}"
    bamFile="${bamFiles[$i]}"
    echo "Processing sample: $sampleName"
    samtools flagstat "${inputdir}${bamFile}" > "${outputdir}${sampleName}_bwa.flagstat.txt"
done

echo "Flagstat processing completed for all samples."
