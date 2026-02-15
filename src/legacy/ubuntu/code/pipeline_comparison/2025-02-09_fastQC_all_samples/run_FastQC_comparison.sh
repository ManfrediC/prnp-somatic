#!/bin/bash
set -e

# Number of threads for trimming
cores=4

# Adapter file for Trimmomatic
adapterFile="/home/mcarta/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"

# Output directory for all produced files (trimmed FASTQs and FastQC results)
outputDir="/add/results/pipeline_comparison/2025-02-09_FastQC_samples/"
mkdir -p "$outputDir"

# Define the input directories containing raw FASTQ files.
inputDirs=(
  "/add/seq_data/2021-10-19_first_CJD_seq/"
  "/add/seq_data/2023-06-18_CJD_16_samples/"
  "/add/seq_data/2023-06-18_CJD_8_samples/"
)

# Loop over each input directory
for dir in "${inputDirs[@]}"; do
  # Process each R1 file (paired-end files assumed to follow a naming convention with _R1_001.fastq.gz)
  for read1 in $dir*_R1_001.fastq.gz; do
    # If no files match, skip to the next directory.
    [ -e "$read1" ] || continue

    # Extract the base filename.
    base=$(basename "$read1")
    # Use regex to extract the sample name.
    # For example, "312091_5-Ctrl1_G01_S6_R1_001.fastq.gz" yields "Ctrl1".
    sample=$(echo "$base" | sed -E 's/.*-([^_]+)_.*/\1/')
    
    # Determine the corresponding R2 file by replacing _R1_001.fastq.gz with _R2_001.fastq.gz.
    read2="${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
    echo "-----------------------------------------"
    echo "Processing sample: $sample from directory: $dir"
    
    #### Aggressive Trimming (parameters: LEADING:10, TRAILING:10, SLIDINGWINDOW:4:25, AVGQUAL:30, MINLEN:80) ####
    
    # Set output file names for aggressive trimming (paired and unpaired)
    aggro_R1="${outputDir}${sample}_aggro_R1.trimmed.fastq"
    aggro_R2="${outputDir}${sample}_aggro_R2.trimmed.fastq"
    aggro_R1_unpaired="${outputDir}${sample}_aggro_R1.unpaired.fastq"
    aggro_R2_unpaired="${outputDir}${sample}_aggro_R2.unpaired.fastq"
    
    echo "Running aggressive trimming for sample $sample..."
    trimmomatic PE -threads $cores -phred33 "$read1" "$read2" \
       "$aggro_R1" "$aggro_R1_unpaired" "$aggro_R2" "$aggro_R2_unpaired" \
       ILLUMINACLIP:"$adapterFile":1:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:25 AVGQUAL:30 MINLEN:80
       
    echo "Running FastQC on aggressive trimmed files for sample $sample..."
    fastqc -o "$outputDir" "$aggro_R1" "$aggro_R2"
    
    # Remove aggressive trimmed FASTQ files to save disk space.
    rm -f "$aggro_R1" "$aggro_R2" "$aggro_R1_unpaired" "$aggro_R2_unpaired"
    
    #### Original Trimming (parameters: LEADING:5, TRAILING:5, SLIDINGWINDOW:5:20, AVGQUAL:30, HEADCROP:0, MINLEN:80) ####
    
    # Set output file names for original trimming
    orig_R1="${outputDir}${sample}_original_R1.trimmed.fastq"
    orig_R2="${outputDir}${sample}_original_R2.trimmed.fastq"
    orig_R1_unpaired="${outputDir}${sample}_original_R1.unpaired.fastq"
    orig_R2_unpaired="${outputDir}${sample}_original_R2.unpaired.fastq"
    
    echo "Running original trimming for sample $sample..."
    trimmomatic PE -threads $cores -phred33 "$read1" "$read2" \
       "$orig_R1" "$orig_R1_unpaired" "$orig_R2" "$orig_R2_unpaired" \
       ILLUMINACLIP:"$adapterFile":1:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:5:20 AVGQUAL:30 HEADCROP:0 MINLEN:80
       
    echo "Running FastQC on original trimmed files for sample $sample..."
    fastqc -o "$outputDir" "$orig_R1" "$orig_R2"
    
    # Remove original trimmed FASTQ files to save disk space.
    rm -f "$orig_R1" "$orig_R2" "$orig_R1_unpaired" "$orig_R2_unpaired"
    
    echo "Finished processing sample $sample."
  done
done

echo "All samples processed. The output directory ($outputDir) now contains:"
echo "- FastQC reports for every sample, with filenames indicating the sample name and trimming strategy ('aggro' or 'original')."
