#!/usr/bin/env bash
# ------------------------------------------------------------
# Merge six single-sample Mutect2 VCFs (Ctrl1-5 & Ctrl7) into
# one multi-sample VCF suitable for CreateSomaticPanelOfNormals.
# ------------------------------------------------------------
# REQUIREMENTS
#   • bcftools =1.10 in $PATH (for merge + indexing)
# ------------------------------------------------------------
set -euo pipefail

# ---------------------- Paths -------------------------------
rawDir="/add/results/PoN/controls_tumor_only"   # input VCFs
outDir="/add/results/PoN"                        # output location
mkdir -p "$outDir"

mergedVCF="$outDir/controls_multisample.vcf.gz"

# ------------------ Merge & index ---------------------------
echo "[1/2] Merging control VCFs into multi-sample file …"
bcftools merge -m all -O z -o "$mergedVCF" \
  "$rawDir/Ctrl1_pon_raw.vcf.gz" \
  "$rawDir/Ctrl2_pon_raw.vcf.gz" \
  "$rawDir/Ctrl3_pon_raw.vcf.gz" \
  "$rawDir/Ctrl4_pon_raw.vcf.gz" \
  "$rawDir/Ctrl5_pon_raw.vcf.gz" \
  "$rawDir/Ctrl7_pon_raw.vcf.gz"

echo "[2/2] Indexing merged VCF …"
bcftools index -f "$mergedVCF"

echo "Merged multi-sample VCF ready: $mergedVCF"

# Optional: quick site count
if command -v bcftools &>/dev/null; then
  echo -n "Number of sites in merged VCF: "
  bcftools view -H "$mergedVCF" | wc -l
fi
