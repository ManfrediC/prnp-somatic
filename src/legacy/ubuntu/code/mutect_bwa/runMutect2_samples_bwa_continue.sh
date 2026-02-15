#!/bin/bash
# runMutect2_samples_bwa_continue.sh
# This script processes only those aggressive pipeline BAM files (in /add/seq_data/2025-02-07_samples_bwa/)
# for which a corresponding Mutect2 VCF file (in /add/results/mutect_output/2025-02-08_samples_bwa/) does not exist.
#
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
    outputVCF="$outDir/${sample}.Mutect.vcf.gz"
    
    if [ -f "$outputVCF" ]; then
        echo "Skipping sample: $sample (already processed)"
    else
        echo "-----------------------------------------"
        echo "Processing sample: $sample"
        echo "Input file: $bamFile"
        
        gatk Mutect2 \
           -R "$fastaFile" \
           -I "$bamFile" \
           --germline-resource "$germlineResource" \
           -O "$outputVCF"
    fi
done

echo "Mutect2 processing completed for all missing samples."
