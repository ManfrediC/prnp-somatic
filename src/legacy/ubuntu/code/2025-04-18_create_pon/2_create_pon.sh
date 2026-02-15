#!/usr/bin/env bash
# ------------------------------------------------------------
# Build the final Panel of Normals (PoN) VCF from the merged
# multi-sample control VCF produced in the previous step.
# Compatible with GATK v4.5.0.0, which requires a single -V.
# ------------------------------------------------------------
# REQUIREMENTS
#   • GATK = 4.5.0.0 in $PATH
# ------------------------------------------------------------
set -euo pipefail

# ---------------------- Paths -------------------------------
dbDir="/home/mcarta/databases"                     # reference & resources
fastaFile="$dbDir/chr2_chr4_chr20.fasta"
mergedVCF="/add/results/PoN/controls_multisample.vcf.gz"  # input (from bcftools)
ponVCF="/add/results/PoN/CJD_controls_PoN.vcf.gz"         # output PoN

# -------------------- Build PoN ----------------------------
echo "[1/1] Running CreateSomaticPanelOfNormals …"

gatk CreateSomaticPanelOfNormals \
    -V "$mergedVCF" \
    -R "$fastaFile" \
    -O "$ponVCF"

echo "PoN created at: $ponVCF"

# quick site count with bcftools
if command -v bcftools &>/dev/null; then
  echo -n "Number of sites in PoN: "
  bcftools view -H "$ponVCF" | wc -l
fi
