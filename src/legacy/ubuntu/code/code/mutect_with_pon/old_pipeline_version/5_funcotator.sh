#!/usr/bin/env bash
# ------------------------------------------------------------
# Step 5 - Annotate each normalised VCF with Funcotator.
# ------------------------------------------------------------
# REQUIREMENTS
#   - GATK 4.5.0.0 (Funcotator) in $PATH
#   - Funcotator data sources (hg38)
# ------------------------------------------------------------
set -euo pipefail

# ----------------------------------------
# Paths and resources
# ----------------------------------------
dbDir="/home/mcarta/databases"
fastaFile="$dbDir/chr2_chr4_chr20.fasta"
funcotatorData="/add/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38"

# Root results directory from previous steps
rootDir="/add/results/mutect_output/2025-04-18_CJD_with_PoN"

# Input / output dirs
normDir="$rootDir/norm_vcf"
annotDir="/add/results/annotation_output/funcotator_PoN"
mkdir -p "$annotDir"

# ----------------------------------------
# Loop over normalised VCFs
# ----------------------------------------
shopt -s nullglob
norm_vcfs=("$normDir"/*.normalized.vcf.gz)
if [ ${#norm_vcfs[@]} -eq 0 ]; then
  echo "[ERROR] No *.normalized.vcf.gz files found in $normDir" >&2
  exit 1
fi

for vcf in "${norm_vcfs[@]}"; do
    sample=$(basename "$vcf" ".normalized.vcf.gz")
    echo "-----------------------------------------"
    echo "Annotating sample: $sample"
    echo "Input VCF: $vcf"

    output_file="$annotDir/${sample}.funcotated.vcf.gz"

    gatk Funcotator \
        --variant "$vcf" \
        --reference "$fastaFile" \
        --ref-version hg38 \
        --data-sources-path "$funcotatorData" \
        --output "$output_file" \
        --output-file-format VCF

    # Ensure output is bgzipped and indexed (Funcotator may output plain VCF)
    if [[ "$output_file" != *.gz ]]; then
        bgzip -f "$output_file"
        output_file+=".gz"
    fi
    tabix -f -p vcf "$output_file"

done

echo "Funcotator annotation completed. Output written to $annotDir"
