#!/bin/bash
#
# This script processes paired-end samples from:
# /add/seq_data/2021-02-03_sequencing_of_dilutions/
#
# Original samples and corresponding files:
#
# 1. Sample: NA100_undil
#    R1: 20210203.0-o23928_1_1-NA100_undil_D01_R1.fastq.gz
#    R2: 20210203.0-o23928_1_1-NA100_undil_D01_R2.fastq.gz
#
# 2. Sample: NA99A1_undil
#    R1: 20210203.0-o23928_1_2-NA99A1_undil_D02_R1.fastq.gz
#    R2: 20210203.0-o23928_1_2-NA99A1_undil_D02_R2.fastq.gz
#
# 3. Sample: NA99A1_1to5
#    R1: 20210203.0-o23928_1_3-NA99A1_1to5_E01_R1.fastq.gz
#    R2: 20210203.0-o23928_1_3-NA99A1_1to5_E01_R2.fastq.gz
#
# 4. Sample: NA100_1to10
#    R1: 20210203.0-o23928_1_4-NA100_1to10_B02_R1.fastq.gz
#    R2: 20210203.0-o23928_1_4-NA100_1to10_B02_R2.fastq.gz
#
# 5. Sample: A100_1to2
#    R1: 20210203.0-o23928_1_5-A100_1to2_F01_R1.fastq.gz
#    R2: 20210203.0-o23928_1_5-A100_1to2_F01_R2.fastq.gz
#
# 6. Sample: NA100_1to2
#    R1: 20210203.0-o23928_1_6-NA100_1to2_F02_R1.fastq.gz
#    R2: 20210203.0-o23928_1_6-NA100_1to2_F02_R2.fastq.gz
#
# 7. Sample: NA995A05_undil
#    R1: 20210203.0-o23928_1_7-NA995A05_undil_G01_R1.fastq.gz
#    R2: 20210203.0-o23928_1_7-NA995A05_undil_G01_R2.fastq.gz
#
# Output directory for all files:
#   /add/seq_data/2025-02-05_dilutions_bwa/

# Define directories
inDir="/add/seq_data/2021-02-03_sequencing_of_dilutions"
outDir="/add/seq_data/2025-02-05_dilutions_bwa/"
mkdir -p "$outDir"

# Define sample arrays
samples=("NA100_undil" "NA99A1_undil" "NA99A1_1to5" "NA100_1to10" "A100_1to2" "NA100_1to2" "NA995A05_undil")
r1Files=(
  "$inDir/20210203.0-o23928_1_1-NA100_undil_D01_R1.fastq.gz"
  "$inDir/20210203.0-o23928_1_2-NA99A1_undil_D02_R1.fastq.gz"
  "$inDir/20210203.0-o23928_1_3-NA99A1_1to5_E01_R1.fastq.gz"
  "$inDir/20210203.0-o23928_1_4-NA100_1to10_B02_R1.fastq.gz"
  "$inDir/20210203.0-o23928_1_5-A100_1to2_F01_R1.fastq.gz"
  "$inDir/20210203.0-o23928_1_6-NA100_1to2_F02_R1.fastq.gz"
  "$inDir/20210203.0-o23928_1_7-NA995A05_undil_G01_R1.fastq.gz"
)
r2Files=(
  "$inDir/20210203.0-o23928_1_1-NA100_undil_D01_R2.fastq.gz"
  "$inDir/20210203.0-o23928_1_2-NA99A1_undil_D02_R2.fastq.gz"
  "$inDir/20210203.0-o23928_1_3-NA99A1_1to5_E01_R2.fastq.gz"
  "$inDir/20210203.0-o23928_1_4-NA100_1to10_B02_R2.fastq.gz"
  "$inDir/20210203.0-o23928_1_5-A100_1to2_F01_R2.fastq.gz"
  "$inDir/20210203.0-o23928_1_6-NA100_1to2_F02_R2.fastq.gz"
  "$inDir/20210203.0-o23928_1_7-NA995A05_undil_G01_R2.fastq.gz"
)

# Common parameters and resource files
cores=4
dbDir="/home/mcarta/databases"
fastaFile="$dbDir/chr2_chr4_chr20.fasta"
adapterFile="/home/mcarta/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"
# Aggressive trimming parameters:
trimmString="ILLUMINACLIP:${adapterFile}:1:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:25 AVGQUAL:30 MINLEN:80"

# Process each sample
for i in "${!samples[@]}"; do
    # Skip the first three samples that do not need to be rerun.
    if [ $i -lt 3 ]; then
        echo "Skipping sample ${samples[$i]}"
        continue
    fi

    sample="${samples[$i]}"
    read1="${r1Files[$i]}"
    read2="${r2Files[$i]}"

    # Use sample name as file prefix for output files
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
    echo "-----------------------------------------"
done

echo "All samples processed."
