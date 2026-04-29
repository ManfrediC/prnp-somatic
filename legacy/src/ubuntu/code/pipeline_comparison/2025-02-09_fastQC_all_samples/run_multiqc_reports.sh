#!/bin/bash
set -e

# Directory where FastQC outputs reside (from both pipelines)
fastqcDir="/add/results/pipeline_comparison/2025-02-09_FastQC_samples/"

# We'll output the MultiQC reports to the same directory
multiqcOutputDir="$fastqcDir"

# Create temporary directories for MultiQC input
mkdir -p "$fastqcDir/multiqc_bwa"
mkdir -p "$fastqcDir/multiqc_orig"

echo "Collecting aggressive (bwa) FastQC files..."
# Symlink FastQC output files (e.g., HTML and ZIP) containing 'aggro' in their names into multiqc_bwa
find "$fastqcDir" -maxdepth 1 -type f -name "*aggro*.html" -exec ln -sf {} "$fastqcDir/multiqc_bwa/" \;
find "$fastqcDir" -maxdepth 1 -type f -name "*aggro*.zip" -exec ln -sf {} "$fastqcDir/multiqc_bwa/" \;

echo "Collecting original FastQC files..."
# Symlink FastQC output files containing 'original' in their names into multiqc_orig
find "$fastqcDir" -maxdepth 1 -type f -name "*original*.html" -exec ln -sf {} "$fastqcDir/multiqc_orig/" \;
find "$fastqcDir" -maxdepth 1 -type f -name "*original*.zip" -exec ln -sf {} "$fastqcDir/multiqc_orig/" \;

echo "Running MultiQC for aggressive (bwa) samples..."
multiqc "$fastqcDir/multiqc_bwa" -o "$multiqcOutputDir" -n multiqc_bwa_report.html

echo "Running MultiQC for original samples..."
multiqc "$fastqcDir/multiqc_orig" -o "$multiqcOutputDir" -n multiqc_original_report.html

# Clean up temporary directories
rm -rf "$fastqcDir/multiqc_bwa" "$fastqcDir/multiqc_orig"

echo "MultiQC reports generated:"
echo "Aggressive (bwa): $multiqcOutputDir/multiqc_bwa_report.html"
echo "Original:         $multiqcOutputDir/multiqc_original_report.html"
