#!/bin/bash
# 2025-02-08_mutect_dilutions_bwa.sh
# This script runs GATK Mutect2 on dilution sample BAM files from the aggressive pipeline.
#
# Expected input files in /add/seq_data/2025-02-07_dilutions_bwa/:
#   NA100_undil.bwa.picard.markedDup.recal.bam
#   NA100_1to10.bwa.picard.markedDup.recal.bam
#   NA99A1_undil.bwa.picard.markedDup.recal.bam
#   A100_1to2.bwa.picard.markedDup.recal.bam
#   NA99A1_1to5.bwa.picard.markedDup.recal.bam
#   NA995A05_undil.bwa.picard.markedDup.recal.bam
#   NA100_1to2.bwa.picard.markedDup.recal.bam
#
# Output files will be written to:
#   /add/results/mutect_output/2025-02-08_dilutions_bwa/
#
# Resources
dbDir="/home/mcarta/databases"
fastaFile="$dbDir/chr2_chr4_chr20.fasta"
germlineResource="$dbDir/somatic-hg38_af-only-gnomad.hg38.vcf"

# Directories for input and output
inputDir="/add/seq_data/2025-02-07_dilutions_bwa"
outDir="/add/results/mutect_output/2025-02-08_dilutions_bwa"
mkdir -p "$outDir"

# Define arrays for the dilution samples and their corresponding BAM files.
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

# Loop over each dilution sample and run Mutect2.
for i in "${!samples[@]}"; do
    sample="${samples[$i]}"
    bamFile="${bamFiles[$i]}"
    inputFile="$inputDir/$bamFile"
    echo "-----------------------------------------"
    echo "Processing dilution sample: $sample"
    echo "Input file: $inputFile"
    
    gatk Mutect2 \
       -R "$fastaFile" \
       -I "$inputFile" \
       --germline-resource "$germlineResource" \
       -O "$outDir/${sample}.Mutect.vcf.gz"
done

echo "Mutect2 processing completed for all dilution samples."
