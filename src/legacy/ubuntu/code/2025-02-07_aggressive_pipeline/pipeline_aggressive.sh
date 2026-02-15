#!/bin/bash
# Usage:
#   bash aggressive_pipeline_two_samples.sh
#
# This script processes two paired-end samples:
#   Sample "CJD14":
#     R1: /add/seq_data/2023-06-18_CJD_16_samples/312092_04-CJD14_B02_S15_R1_001.fastq.gz
#     R2: /add/seq_data/2023-06-18_CJD_16_samples/312092_04-CJD14_B02_S15_R2_001.fastq.gz
#
#   Sample "CJD21":
#     R1: /add/seq_data/2023-06-18_CJD_16_samples/312092_08-CJD21_D02_S1_R1_001.fastq.gz
#     R2: /add/seq_data/2023-06-18_CJD_16_samples/312092_08-CJD21_D02_S1_R2_001.fastq.gz
#
# All outputs will be saved into:
#   /add/seq_data/2025-02-07_CJD_aggressive_test/

# Define output directory and create it if needed
outDir="/add/seq_data/2025-02-07_CJD_aggressive_test/"
mkdir -p $outDir

# Define input file paths for each sample
# Arrays for sample names, R1 files, and R2 files
samples=("CJD14" "CJD21")
r1Files=(
  "/add/seq_data/2023-06-18_CJD_16_samples/312092_04-CJD14_B02_S15_R1_001.fastq.gz"
  "/add/seq_data/2023-06-18_CJD_16_samples/312092_08-CJD21_D02_S1_R1_001.fastq.gz"
)
r2Files=(
  "/add/seq_data/2023-06-18_CJD_16_samples/312092_04-CJD14_B02_S15_R2_001.fastq.gz"
  "/add/seq_data/2023-06-18_CJD_16_samples/312092_08-CJD21_D02_S1_R2_001.fastq.gz"
)

# Define common parameters and resources
cores=4
dbDir="/home/mcarta/databases"
fastaFile="$dbDir/chr2_chr4_chr20.fasta"
adapterFile="/home/mcarta/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"
# Aggressive trimming parameters: increased LEADING and TRAILING thresholds and a more stringent SLIDINGWINDOW.
trimmString="ILLUMINACLIP:${adapterFile}:1:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:25 AVGQUAL:30 MINLEN:80"

# Loop over each sample and process
for i in "${!samples[@]}"; do
    sample="${samples[$i]}"
    read1="${r1Files[$i]}"
    read2="${r2Files[$i]}"
    
    # Use sample name as file prefix
    filePrefix="$sample"
    echo "Processing sample: $sample"
    
    ##############################
    # 1. Aggressive Trimming Step
    ##############################
    echo "Starting aggressive trimming for $sample..."
    trimmomatic PE -threads $cores -phred33 \
        "$read1" "$read2" \
        "$outDir/${filePrefix}.R1.trimmed.fastq" "$outDir/${filePrefix}.R1.unpaired.fastq" \
        "$outDir/${filePrefix}.R2.trimmed.fastq" "$outDir/${filePrefix}.R2.unpaired.fastq" \
        $trimmString

    ##############################
    # 2. Alignment Using bwa-mem
    ##############################
    echo "Starting alignment with bwa-mem for $sample..."
    bwa mem -t $cores $fastaFile \
        "$outDir/${filePrefix}.R1.trimmed.fastq" \
        "$outDir/${filePrefix}.R2.trimmed.fastq" > "$outDir/${filePrefix}.bwa.sam"

    ###########################################
    # 3. Convert SAM to Sorted BAM and Index
    ###########################################
    echo "Converting SAM to sorted BAM for $sample..."
    samtools view -hbS "$outDir/${filePrefix}.bwa.sam" | samtools sort -o "$outDir/${filePrefix}.bwa.sorted.bam"
    samtools index "$outDir/${filePrefix}.bwa.sorted.bam"

    #################################################
    # 4. Add Read Groups Using Picard (Picard tools)
    #################################################
    echo "Adding read groups for $sample..."
    picard AddOrReplaceReadGroups \
        I="$outDir/${filePrefix}.bwa.sorted.bam" \
        O="$outDir/${filePrefix}.bwa.picard.bam" \
        SORT_ORDER=coordinate \
        RGID="$filePrefix" \
        RGLB=Paired_end \
        RGPL=illumina \
        RGSM="$sample" \
        RGPU=illumina \
        USE_JDK_DEFLATER=true \
        USE_JDK_INFLATER=true

    #################################################
    # 5. Mark Duplicates Using Picard
    #################################################
    echo "Marking duplicates for $sample..."
    picard MarkDuplicates \
        I="$outDir/${filePrefix}.bwa.picard.bam" \
        O="$outDir/${filePrefix}.bwa.picard.markedDup.bam" \
        M="$outDir/${filePrefix}.bwa.picard.markedDup.metrics" \
        USE_JDK_DEFLATER=true \
        USE_JDK_INFLATER=true

    #################################################
    # 6. GATK BaseRecalibrator and Apply BQSR
    #################################################
    echo "Starting GATK BaseRecalibrator for $sample..."
    gatk --java-options "-Xmx8g" BaseRecalibrator \
        -I "$outDir/${filePrefix}.bwa.picard.markedDup.bam" \
        -R $fastaFile \
        --known-sites "$dbDir/dbsnp_146.hg38.vcf.gz" \
        --known-sites "$dbDir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
        -O "$outDir/${filePrefix}.bwa.recal_data.table"

    echo "Applying Base Quality Score Recalibration for $sample..."
    gatk --java-options "-Xmx8g" ApplyBQSR \
        -R $fastaFile \
        -I "$outDir/${filePrefix}.bwa.picard.markedDup.bam" \
        --bqsr-recal-file "$outDir/${filePrefix}.bwa.recal_data.table" \
        -O "$outDir/${filePrefix}.bwa.picard.markedDup.recal.bam"

    echo "Sample $sample processing complete. Output files are in $outDir."
    echo "-------------------------------------------"
done

echo "All samples processed."
