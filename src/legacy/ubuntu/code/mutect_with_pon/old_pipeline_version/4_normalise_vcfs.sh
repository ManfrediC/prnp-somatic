#!/usr/bin/env bash
# ------------------------------------------------------------
# Step 4 - Normalise each filtered VCF:
#          * split multiallelic sites (-m -any)
#          * left-align indels relative to reference
#          * bgzip-compress and index
# ------------------------------------------------------------
# REQUIREMENTS
#   - bcftools =1.10 in $PATH
# ------------------------------------------------------------
set -euo pipefail

# --------------------------------------------------
# Resources and directories
# --------------------------------------------------
# Reference FASTA limited to chr2, chr4, chr20
fastaFile="/home/mcarta/databases/chr2_chr4_chr20.fasta"

rootDir="/add/results/mutect_output/2025-04-18_CJD_with_PoN"
filteredDir="$rootDir/filtered_vcf"   # input from previous step
normDir="$rootDir/norm_vcf"            # output directory

mkdir -p "$normDir"

# --------------------------------------------------
# Loop over each filtered VCF and normalise
# --------------------------------------------------
shopt -s nullglob
filteredVcfs=("$filteredDir"/*.filtered.vcf.gz)

if [ ${#filteredVcfs[@]} -eq 0 ]; then
  echo "[ERROR] No *.filtered.vcf.gz files found in $filteredDir"
  exit 1
fi

echo "[INFO] Normalising ${#filteredVcfs[@]} VCFs"

for vcf in "${filteredVcfs[@]}"; do
  sample=$(basename "$vcf" .filtered.vcf.gz)
  outVCF="$normDir/${sample}.normalized.vcf.gz"

  echo "[INFO] Processing $sample"

  # Split multiallelics and left-align indels
  bcftools norm -m -any -f "$fastaFile" "$vcf" -Oz -o "$outVCF"

  # Index the resulting VCF
  bcftools index -t -f "$outVCF"

done

echo "Normalised VCFs written to $normDir"
