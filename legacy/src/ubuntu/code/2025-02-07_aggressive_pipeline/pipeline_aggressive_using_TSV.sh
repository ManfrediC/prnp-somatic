#!/bin/bash
# Usage: bash trial_pipeline.sh /path/to/input_directory
# This script processes one sample from dataset.tsv using aggressive trimming and bwa-mem

# Check for input directory parameter
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 input_directory"
    exit 1
fi

# Set working directory
inputDir=$1
cd $inputDir

# Parse the dataset file (assumes dataset.tsv exists in the directory)
# Skip header line using grep -v on the column name
sampleNames=(dum $(cut -f 1 dataset.tsv | grep -v Name))
case1=(dum $(cut -f 2 dataset.tsv | grep -v Read1))
case2=(dum $(cut -f 3 dataset.tsv | grep -v Read2))

# For testing, choose the first sample (after the dummy element)
sampleIndex=1
fileName=$(basename ${case1[$sampleIndex]} _R1.fastq.gz)
sampleName=${sampleNames[$sampleIndex]}
echo "Processing sample: $sampleName"

# Define resources and parameters
cores=4
dbDir="/home/mcarta/databases"
fastaFile="$dbDir/chr2_chr4_chr20.fasta"
adapterFile="/home/mcarta/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"

# Aggressive trimming parameters: increase quality thresholds (LEADING, TRAILING, SLIDINGWINDOW)
trimmString="ILLUMINACLIP:$adapterFile:1:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:25 AVGQUAL:30 MINLEN:80"
read1=${case1[$sampleIndex]}
read2=${case2[$sampleIndex]}

echo "Starting aggressive trimming with Trimmomatic..."
trimmomatic PE -threads $cores -phred33 $read1 $read2 \
    $fileName.R1.trimmed.fastq  $fileName.R1.unpaired \
    $fileName.R2.trimmed.fastq $fileName.R2.unpaired $trimmString

echo "Trimming completed. Starting alignment with bwa-mem..."

# Alignment using bwa-mem (using the trimmed paired reads)
bwa mem -t $cores $fastaFile $fileName.R1.trimmed.fastq $fileName.R2.trimmed.fastq > $fileName.bwa.sam

echo "Alignment completed. Converting SAM to sorted BAM..."

# Convert SAM to BAM, sort, and index using samtools
samtools view -hbS $fileName.bwa.sam | samtools sort -o $fileName.bwa.sorted.bam
samtools index $fileName.bwa.sorted.bam

echo "BAM file generated and indexed. Adding read groups with Picard..."

# Add read groups using Picard
picard AddOrReplaceReadGroups I=$fileName.bwa.sorted.bam \
    O=$fileName.bwa.picard.bam SORT_ORDER=coordinate RGID=$fileName \
    RGLB=Paired_end RGPL=illumina RGSM=$sampleName RGPU=illumina \
    USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

echo "Read groups added. Marking duplicates with Picard..."

# Mark duplicates
picard MarkDuplicates I=$fileName.bwa.picard.bam \
    O=$fileName.bwa.picard.markedDup.bam M=$fileName.bwa.picard.markedDup.metrics \
    USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

echo "Duplicates marked. Starting GATK BaseRecalibrator..."

# GATK BaseRecalibrator (ensure you have the known sites files in your databases directory)
gatk --java-options "-Xmx8g" BaseRecalibrator \
    -I $fileName.bwa.picard.markedDup.bam \
    -R $fastaFile \
    --known-sites $dbDir/dbsnp_146.hg38.vcf.gz \
    --known-sites $dbDir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O $fileName.bwa.recal_data.table

echo "Base recalibration table generated. Applying BQSR..."

# Apply Base Quality Score Recalibration
gatk --java-options "-Xmx8g" ApplyBQSR \
    -R $fastaFile \
    -I $fileName.bwa.picard.markedDup.bam \
    --bqsr-recal-file $fileName.bwa.recal_data.table \
    -O $fileName.bwa.picard.markedDup.recal.bam

echo "BQSR applied. Pipeline trial for sample $sampleName is complete."
