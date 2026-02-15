#!/usr/bin/env bash
# ------------------------------------------------------------
# Step 2 - Apply FilterMutectCalls to each CJD raw VCF using
#          the orientation-bias priors. PoN and germline AFs
#          were already used during Mutect2 calling.
# ------------------------------------------------------------
# REQUIREMENTS
#   - GATK 4.5.0.0 (or later) in $PATH
# ------------------------------------------------------------
set -euo pipefail

# --------------------------------------------------
# Resources and directories (edit if your layout differs)
# --------------------------------------------------
# Reference FASTA limited to chr2, chr4, chr20
fastaFile="/home/mcarta/databases/chr2_chr4_chr20.fasta"

rootDir="/add/results/mutect_output/2025-04-18_CJD_with_PoN"
rawDir="$rootDir"                     # produced by the Mutect2 step
orientationDir="$rootDir/orientation" # produced by LearnReadOrientationModel
filtDir="$rootDir/filtered_vcf"

mkdir -p "$filtDir"

# --------------------------------------------------
# Loop over every raw VCF and apply FilterMutectCalls
# --------------------------------------------------
shopt -s nullglob
rawVcfs=("$rawDir"/*.raw.vcf.gz)

if [ ${#rawVcfs[@]} -eq 0 ]; then
  echo "[ERROR] No *.raw.vcf.gz files found in $rawDir"
  exit 1
fi

echo "[INFO] Running FilterMutectCalls on ${#rawVcfs[@]} samples"

for vcf in "${rawVcfs[@]}"; do
  sample=$(basename "$vcf" .raw.vcf.gz)
  orient="$orientationDir/${sample}_orientation.tar.gz"
  statsFile="$rawDir/${sample}.raw.vcf.gz.stats"

  if [ ! -f "$orient" ]; then
    echo "[WARN] Orientation model missing for $sample — skipping"
    continue
  fi

  if [ ! -f "$statsFile" ]; then
    echo "[WARN] Stats file missing for $sample — skipping"
    continue
  fi

  echo "[INFO] Processing $sample"

  gatk FilterMutectCalls \
      -R "$fastaFile" \
      -V "$vcf" \
      --stats "$statsFile" \
      --orientation-bias-artifact-priors "$orient" \
      -O "$filtDir/${sample}.filtered.vcf.gz"

done

echo "Filtered VCFs written to $filtDir"
