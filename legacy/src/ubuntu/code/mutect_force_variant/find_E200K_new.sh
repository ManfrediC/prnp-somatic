
#based on https://bioinformatics.stackexchange.com/questions/22237/how-do-i-quantify-a-specific-somatic-variant

fastaFile="/home/mcarta/databases/hg38.fa"

mutationName="E200K"
positionOfInterest="4699818"
baseOfInterest="chr20:4699818-4699818"


sampleName="CJD30"
inputFile="/add/2023-06-18_CJD_16_samples/312092_13-CJD30_G01_S16_R1_001.fastq.gz.picard.markedDup.recal.bam"

outputDir="/add/2024-03-04_force_mutect/output/"

mutationVCF="/add/2024-03-04_force_mutect/E200K.vcf.gz"


# Step 1: VCF file containing the variant of interest is required -> force Mutect to call it
# the issue is that the allelic depth parameter doesn't contain all reads (only about 100  instead of 400)

fastaFile="/home/mcarta/databases/hg38.fa"
sampleName="CJD30"
mutationName="E200K"

gatk Mutect2 \
  -R $fastaFile \
  --alleles $mutationVCF  \
  -L $mutationVCF \
  -I $inputFile \
  -O ${outputDir}/${sampleName}_${mutationName}_frequency.vcf 
