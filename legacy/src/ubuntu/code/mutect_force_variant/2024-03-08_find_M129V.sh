#!/bin/bash

###this script forces Mutect2 to call D178N

#database
fastaFile="/home/mcarta/databases/hg38.fa"

#define mutation
mutationName="M129V"
positionOfInterest="4699605"
baseOfInterest="chr20:4699605-4699605"

#VCF file that contains this defined mutation of interest
mutationVCF="/add/variant_VCF/M129V.vcf.gz"

#define sample to be processed
sampleName="CJD6"
inputFile="/add/seq_data/2021-10-19_first_CJD_seq/20211027.0-o26424_1_3-CJD6_fr.picard.markedDup.recal.bam"

#where to save the output VCF files (one for every sample)
outputDir="/add/results/2024-03-04_force_mutect/2024-03-08_M129V/"

# Step 1: VCF file containing the variant of interest is required -> force Mutect to call it

gatk Mutect2 \
  -R "$fastaFile" \
  --alleles "$mutationVCF"  \
  -L "$mutationVCF" \
  -I "$inputFile" \
  -O "${outputDir}/${sampleName}_${mutationName}_frequency.vcf"