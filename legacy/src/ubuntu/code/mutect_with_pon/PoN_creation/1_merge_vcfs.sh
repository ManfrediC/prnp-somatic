#!/usr/bin/env bash
# ------------------------------------------------------------
# Date: 2025-05-04
# Part of variant calling pipeline using a PoN 
# Merge six single-sample filtered Mutect2 VCFs (Ctrl1–5 & Ctrl7)
# into one multi-sample VCF for CreateSomaticPanelOfNormals.
# ------------------------------------------------------------
# REQUIREMENTS
#   • bcftools >=1.10 in $PATH (for merge + indexing)
# ------------------------------------------------------------
set -euo pipefail

# ---------------------- Paths -------------------------------
rawDir="/add/results/no_PoN/filtered"  # input VCFs (filtered with FilterMutectCalls)
outDir="/add/results/PoN/panel_of_normals"              # output location
mkdir -p "$outDir"

mergedVCF="$outDir/controls_multisample.vcf.gz"

# ------------------ Merge & index ---------------------------
echo "[1/2] Merging filtered control VCFs into multi-sample file …"
bcftools merge -m all -O z -o "$mergedVCF" \
  "$rawDir/Ctrl1.filtered.vcf.gz" \
  "$rawDir/Ctrl2.filtered.vcf.gz" \
  "$rawDir/Ctrl3.filtered.vcf.gz" \
  "$rawDir/Ctrl4.filtered.vcf.gz" \
  "$rawDir/Ctrl5.filtered.vcf.gz" \
  "$rawDir/Ctrl7.filtered.vcf.gz"

echo "[2/2] Indexing merged VCF …"
bcftools index -f "$mergedVCF"
gatk IndexFeatureFile -I /add/results/PoN/panel_of_normals/controls_multisample.vcf.gz

echo "Merged multi-sample VCF ready: $mergedVCF"

# Optional: quick site count
if command -v bcftools &>/dev/null; then
  echo -n "Number of sites in merged VCF: "
  bcftools view -H "$mergedVCF" | wc -l
fi
