
#the directories have since changed

fastaFile="/home/mcarta/databases/hg38.fa"

mutationName="E200K"
positionOfInterest="4699818"
baseOfInterest="chr20:4699818-4699818"


sampleName="CJD30"
inputFile="/add/2023-06-18_CJD_16_samples/312092_13-CJD30_G01_S16_R1_001.fastq.gz.picard.markedDup.recal.bam"

outputDir="/add/2024-03-03_count_path_variants/"

# Step 1: Extract reads overlapping the position of interest
extractedReads="$outputDir/${sampleName}_${mutationName}_reads.bam"

samtools view -b -h "$inputFile" "$baseOfInterest" > "$extractedReads"


# Step 2: Filter reads supporting the alternate allele

variantVCF="$outputDir/${sampleName}_${mutationName}_variants.vcf.gz"

# output maximum allele depth for reference (REF) and alternate (ALT). Number of reads to consider 1'000'000 (otherwise it defaults to 250).
bcftools mpileup -Ou -f $fastaFile --annotate FORMAT/AD -d 1000000 "$extractedReads" | bcftools call -mv -Oz -o $variantVCF
gunzip -k "$variantVCF"



