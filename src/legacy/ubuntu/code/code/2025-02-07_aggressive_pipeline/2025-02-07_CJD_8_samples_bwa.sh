#!/bin/bash
# Usage:
#   bash pipeline_CJD8_samples.sh
#
# This script processes paired-end samples from:
# /add/seq_data/2023-06-18_CJD_8_samples/
#
# The input files are:
#   312091_1-CJD3_E01_S1_R1_001.fastq.gz
#   312091_1-CJD3_E01_S1_R2_001.fastq.gz
#   312091_2-CJD4_E02_S2_R1_001.fastq.gz
#   312091_2-CJD4_E02_S2_R2_001.fastq.gz
#   312091_3-CJD7_F01_S8_R1_001.fastq.gz
#   312091_3-CJD7_F01_S8_R2_001.fastq.gz
#   312091_4-CJD11_F02_S4_R1_001.fastq.gz
#   312091_4-CJD11_F02_S4_R2_001.fastq.gz
#   312091_5-Ctrl1_G01_S6_R1_001.fastq.gz
#   312091_5-Ctrl1_G01_S6_R2_001.fastq.gz
#   312091_6-Ctrl2_G02_S7_R1_001.fastq.gz
#   312091_6-Ctrl2_G02_S7_R2_001.fastq.gz
#   312091_7-Ctrl3_H01_S5_R1_001.fastq.gz
#   312091_7-Ctrl3_H01_S5_R2_001.fastq.gz
#   312091_8-Ctrl5_H02_S3_R1_001.fastq.gz
#   312091_8-Ctrl5_H02_S3_R2_001.fastq.gz
#
# The sample name is extracted from the filename using a regex.
# For example, from "312091_1-CJD3_E01_S1_R1_001.fastq.gz" the sample name "CJD3" is obtained.
#
# Output files are written to:
#   /add/seq_data/2025-02-05_samples_bwa/

# Define directories
inDir="/add/seq_data/2023-06-18_CJD_8_samples"
outDir="/add/seq_data/2025-02-05_samples_bwa"
mkdir -p "$outDir"

# Define arrays for R1 and R2 files (order matters: first file of R1 pairs with first file of R2)
r1Files=(
  "$inDir/312091_1-CJD3_E01_S1_R1_001.fastq.gz"
  "$inDir/312091_2-CJD4_E02_S2_R1_001.fastq.gz"
  "$inDir/312091_3-CJD7_F01_S8_R1_001.fastq.gz"
  "$inDir/312091_4-CJD11_F02_S4_R1_001.fastq.gz"
  "$inDir/312091_5-Ctrl1_G01_S6_R1_001.fastq.gz"
  "$inDir/312091_6-Ctrl2_G02_S7_R1_001.fastq.gz"
  "$inDir/312091_7-Ctrl3_H01_S5_R1_001.fastq.gz"
  "$inDir/312091_8-Ctrl5_H02_S3_R1_001.fastq.gz"
)
r2Files=(
  "$inDir/312091_1-CJD3_E01_S1_R2_001.fastq.gz"
  "$inDir/312091_2-CJD4_E02_S2_R2_001.fastq.gz"
  "$inDir/312091_3-CJD7_F01_S8_R2_001.fastq.gz"
  "$inDir/312091_4-CJD11_F02_S4_R2_001.fastq.gz"
  "$inDir/312091_5-Ctrl1_G01_S6_R2_001.fastq.gz"
  "$inDir/312091_6-Ctrl2_G02_S7_R2_001.fastq.gz"
  "$inDir/312091_7-Ctrl3_H01_S5_R2_001.fastq.gz"
  "$inDir/312091_8-Ctrl5_H02_S3_R2_001.fastq.gz"
)

# Automatically extract sample names from the R1 filenames using sed.
# The pattern removes everything up to the dash and extracts the text before the next underscore.
samples=()
for file in "${r1Files[@]}"; do
  base=$(basename "$file")
  sample=$(echo "$base" | sed -E 's/.*-([^_]+)_.*/\1/')
  samples+=("$sample")
done

# Print extracted sample names for verification
echo "Extracted sample names:"
printf "%s\n" "${samples[@]}"

# Define common parameters and resource files
cores=4
dbDir="/home/mcarta/databases"
fastaFile="$dbDir/chr2_chr4_chr20.fasta"
adapterFile="/home/mcarta/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"
# Aggressive trimming parameters:
trimmString="ILLUMINACLIP:${adapterFile}:1:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:25 AVGQUAL:30 MINLEN:80"

# Process each sample in a loop
for i in "${!samples[@]}"; do
    sample="${samples[$i]}"
    r1="${r1Files[$i]}"
    r2="${r2Files[$i]}"
    filePrefix="$sample"
    echo "-----------------------------------------"
    echo "Processing sample: $sample"
    
    ##############################
    # 1. Aggressive Trimming Step
    ##############################
    echo "Starting aggressive trimming for $sample..."
    trimmomatic PE -threads $cores -phred33 \
      "$r1" "$r2" \
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
    # 4. Add Read Groups Using Picard
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
done

echo "-----------------------------------------"
echo "All samples processed."
