#!/bin/bash
# runMutect2_samples_bwa.sh
# This script runs GATK Mutect2 on aggressive pipeline output BAM files located in:
#   /add/seq_data/2025-02-07_samples_bwa/
#
# The BAM files are expected to be named in the format:
#   <SampleName>.bwa.picard.markedDup.recal.bam
# For example:
#   CJD31.bwa.picard.markedDup.recal.bam
#   Ctrl2.bwa.picard.markedDup.recal.bam
#
# The sample name is extracted by removing the suffix ".bwa.picard.markedDup.recal.bam".
# Mutect2 outputs VCF files into the output directory:
#   /add/results/mutect_output/2025-02-08_samples_bwa/

# Resources
dbDir="/home/mcarta/databases"
fastaFile="$dbDir/chr2_chr4_chr20.fasta"
germlineResource="$dbDir/somatic-hg38_af-only-gnomad.hg38.vcf"

# Directories
inputDir="/add/seq_data/2025-02-07_samples_bwa"
outDir="/add/results/mutect_output/2025-02-08_samples_bwa"
mkdir -p "$outDir"

# Loop over each BAM file matching the expected pattern
for bamFile in "$inputDir"/*.bwa.picard.markedDup.recal.bam; do
    # Extract the sample name by removing the known suffix from the filename.
    sample=$(basename "$bamFile" ".bwa.picard.markedDup.recal.bam")
    echo "-----------------------------------------"
    echo "Processing sample: $sample"
    echo "Input file: $bamFile"
    
    # Run Mutect2
    gatk Mutect2 \
       -R "$fastaFile" \
       -I "$bamFile" \
       --germline-resource "$germlineResource" \
       -O "$outDir/${sample}.Mutect.vcf.gz"
done

echo "Mutect2 processing completed for all samples."
