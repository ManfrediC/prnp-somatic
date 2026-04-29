#this generates a VCF with Mutect calls for D178N

fastaFile="/home/mcarta/databases/hg38.fa"

mutationName="D178N"
positionOfInterest="4699752"
baseOfInterest="chr20:4699752-4699752"


sampleName="Ctrl1"
inputFile="/add/seq_data/2023-06-18_CJD_8_samples/312091_5-Ctrl1_G01_S6_R1_001.fastq.gz.picard.markedDup.recal.bam"

outputDir="/add/results/2024-03-05_generate_mutation_VCFs/"

# Step 1: Extract reads overlapping the position of interest
extractedReads="$outputDir/${sampleName}_${mutationName}_reads.bam"

samtools view -b -h "$inputFile" "$baseOfInterest" > "$extractedReads"


# Step 2: Filter reads supporting the alternate allele

variantVCF="$outputDir/${sampleName}_${mutationName}_variants.vcf.gz"

# output maximum allele depth for reference (REF) and alternate (ALT). Number of reads to consider 1'000'000 (otherwise it defaults to 250).
bcftools mpileup -Ou -f $fastaFile --annotate FORMAT/AD -d 1000000 "$extractedReads" | bcftools call -mv -Oz -o $variantVCF
gunzip -k "$variantVCF"



