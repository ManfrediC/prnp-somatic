#!/usr/bin/env bash
# ------------------------------------------------------------
# Step 1 - Build an orientation-bias model for each CJD sample
#          using the F1R2 data produced by Mutect2.
# ------------------------------------------------------------
# REQUIREMENTS
#   - GATK 4.5.0.0 (or later) in $PATH
# ------------------------------------------------------------
set -euo pipefail

# --------------------------------------------------
# Directories (edit if your layout differs)
# --------------------------------------------------
rootDir="/add/results/mutect_output/2025-04-18_CJD_with_PoN"
rawDir="$rootDir"          # contains *.f1r2.tar.gz files
orientationDir="$rootDir/orientation"

mkdir -p "$orientationDir"

# --------------------------------------------------
# Loop over every F1R2 archive and learn orientation model
# --------------------------------------------------
shopt -s nullglob
f1r2_files=("$rawDir"/*.f1r2.tar.gz)

if [ ${#f1r2_files[@]} -eq 0 ]; then
  echo "[ERROR] No *.f1r2.tar.gz files found in $rawDir"
  exit 1
fi

for f1 in "${f1r2_files[@]}"; do
  sample=$(basename "$f1" .f1r2.tar.gz)
  echo "[INFO] Processing $sample"

  gatk LearnReadOrientationModel \
      -I "$f1" \
      -O "$orientationDir/${sample}_orientation.tar.gz"

done

echo "Orientation-bias models written to $orientationDir"
