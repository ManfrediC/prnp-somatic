#!/bin/bash

###this script forces Mutect2 to call D178N

#database
fastaFile="/home/mcarta/databases/hg38.fa"

#define mutation
mutationName="D178N"
positionOfInterest="4699752"
baseOfInterest="chr20:4699752-4699752"

#VCF file that contains this defined mutation of interest (once it works: set up a new directory)
mutationVCF="/add/results/2024-03-05_ref_VCFs_for_mutations/D178Nnew.vcf.gz"

#define sample to be processed
sampleName="CJD1"
#inputFile="/add/seq_data/2023-06-18_CJD_8_samples/312091_5-Ctrl1_G01_S6_R1_001.fastq.gz.picard.markedDup.recal.bam"
inputFile="/add/seq_data/2021-10-19_first_CJD_seq/20211027.0-o26424_1_1-CJD1_fr.picard.markedDup.recal.bam"

#where to save the output VCF files (one for every sample)
outputDir="/add/results/2024-03-04_force_mutect/2024-03-08_D178N/"

# Step 1: VCF file containing the variant of interest is required -> force Mutect to call it

gatk Mutect2 \
  -R "$fastaFile" \
  --alleles "$mutationVCF"  \
  -L "$mutationVCF" \
  -I "$inputFile" \
  -O "${outputDir}/${sampleName}_${mutationName}_frequency.vcf"