#!/usr/bin/env bash
# ------------------------------------------------------------
# Run Mutect2 tumour-only on every CJD BAM **with the new PoN**.
# ------------------------------------------------------------
# REQUIREMENTS
#   • GATK = 4.5.0.0 in $PATH
# ------------------------------------------------------------
set -euo pipefail

# --------------------------------------------------
# Resources
# --------------------------------------------------
dbDir="/home/mcarta/databases"

# Reference FASTA limited to chr2, chr4, chr20
fastaFile="$dbDir/chr2_chr4_chr20.fasta"

# gnomAD germline resource (AF-only)
germlineResource="$dbDir/somatic-hg38_af-only-gnomad.hg38.vcf"

# Interval list for probe-capture targets
intervalList="$dbDir/capture_targets.interval_list"

# Panel of Normals created from Controls
ponVCF="/add/results/PoN/CJD_controls_PoN.vcf.gz"

# --------------------------------------------------
# Directories
# --------------------------------------------------
inputDir="/add/seq_data/2025-02-07_samples_bwa"            # CJD BAMs
outRoot="/add/results/mutect_output/2025-04-18_CJD_with_PoN" # output root
mkdir -p "$outRoot"

# --------------------------------------------------
# Loop over CJD samples
# --------------------------------------------------
shopt -s nullglob
cjdBAMs=("$inputDir"/CJD*.bwa.picard.markedDup.recal.bam)
if (( ${#cjdBAMs[@]} == 0 )); then
  echo "[ERROR] No CJD BAMs found in $inputDir (pattern CJD*.bwa.picard.markedDup.recal.bam)" >&2
  exit 1
fi

for bamFile in "${cjdBAMs[@]}"; do
    sample=$(basename "$bamFile" ".bwa.picard.markedDup.recal.bam")
    echo "-----------------------------------------"
    echo "Running Mutect2 (PoN) on sample: $sample"
    echo "Input BAM: $bamFile"

    gatk Mutect2 \
        -R "$fastaFile" \
        -I "$bamFile" \
        -tumor-sample "$sample" \
        --panel-of-normals "$ponVCF" \
        --germline-resource "$germlineResource" \
        --af-of-alleles-not-in-resource 0.0000025 \
        --intervals "$intervalList" \
        --max-mnp-distance 0 \
        --f1r2-tar-gz "$outRoot/${sample}.f1r2.tar.gz" \
        -O "$outRoot/${sample}.raw.vcf.gz"

done

echo "Mutect2 tumour-only (with PoN) completed for all CJD samples."
