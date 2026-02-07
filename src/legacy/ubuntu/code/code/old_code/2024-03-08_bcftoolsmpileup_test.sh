
### define the variables
fastaFile="/home/mcarta/databases/hg38.fa"

# mutation of interest
mutationName="D178N"
positionOfInterest="chr20:4699752"
baseOfInterest="chr20:4699752-4699752"

# sample to use
sampleName="Ctrl1"
inputFile="/add/seq_data/2023-06-18_CJD_8_samples/312091_5-Ctrl1_G01_S6_R1_001.fastq.gz.picard.markedDup.recal.bam"

# where to save the output
outputDir="/add/results/2024-03-05_ref_VCFs_for_mutations"



### Step 1: Extract reads overlapping the position of interest (this step works correctly) + indexing

# name of output file
extractedReads="$outputDir/${sampleName}_${mutationName}_reads.bam"

# extract reads
samtools view -b -h "$inputFile" "$baseOfInterest" > "$extractedReads"

# index BAM file -> .bam.bai
samtools index $extractedReads

# Note: I have converted the bam to sam and manually checked it, it does overlap the site of interest -> Step 1 works



### Step 2: Filter reads supporting the alternate allele


# PROBLEM: bcftools will only produce an output if the variant is actually present. In this case, it isn't.


# name of output file
variantVCF="$outputDir/testD178N.vcf.gz"

#bcftools mpileup -Ou \
#    -f $fastaFile \
#    --annotate FORMAT/AD \
#    -d 1000000 \
#    -r chr20:4699752 \
#    "$extractedReads" \
#    -a "FORMAT/AD" \
#    -o $variantVCF.mpileup




#bcftools mpileup -Ou \
#    -f $fastaFile \
#    --annotate FORMAT/AD \
#   -d 1000000 \
#    -r chr20:4699752 \
#    "$extractedReads" \
#    -o $variantVCF.mpileup
    



#bcftools mpileup -Ou \
#    -f $fastaFile \
#   --annotate FORMAT/AD \
#    -d 1000000 \
#    -r chr20:4699752 \
#    "$extractedReads" | \
#    bcftools call -mv -Oz \
#    -o $variantVCF \
#    -c vcf


#bcftools mpileup -Ou \
#    -f $fastaFile \
#    --annotate FORMAT/AD \
#    -d 1000000 \
#    -r chr20:4699752 \
#    "$extractedReads" \
#    -o $variantVCF.mpileup

#echo $variantVCF

#gunzip -k "$variantVCF"
