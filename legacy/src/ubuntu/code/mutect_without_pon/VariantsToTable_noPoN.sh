#!/usr/bin/env bash
# ============================================================
# Extract PASS variants to TSV for CJD+Ctrl samples AND dilutions
# Uses GATK SelectVariants and VariantsToTable
# ============================================================
set -euo pipefail

# Directories
VCF_SAMPLES_DIR="/add/results/no_PoN/annot"      # Funcotator outputs for samples
VCF_DIL_DIR="/add/results/no_PoN/annot/dil"     # Funcotator outputs for dilutions
SELECT_SAMPLES_DIR="/add/results/no_PoN/select_variants/samples"
SELECT_DIL_DIR="/add/results/no_PoN/select_variants/dilutions"
TSV_SAMPLES_DIR="/add/results/no_PoN/variant_tables/samples"
TSV_DIL_DIR="/add/results/no_PoN/variant_tables/dilutions"

# Create output directories
mkdir -p "$SELECT_SAMPLES_DIR" "$SELECT_DIL_DIR" "$TSV_SAMPLES_DIR" "$TSV_DIL_DIR"

# Common fields
INFO_FIELDS=(CHROM POS REF ALT FILTER QUAL MQ QD AF GNOMAD_AF FS BaseQualityRankSumTest MappingQualityRankSumTest AS_SB_TABLE STRANDQ FUNCOTATION)
FORMAT_FIELDS=(GT DP AD F1R2 F2R1 SB)

# Function to process one set
process_set() {
  local vcf_dir=$1 select_dir=$2 tsv_dir=$3

  for vcf in "${vcf_dir}"/*.func.vcf.gz; do
    sample=$(basename "$vcf" .func.vcf.gz)
    echo "Processing $sample"

    # 1) Select PASS variants
    out_vcf="${select_dir}/${sample}.PASS.vcf.gz"
    if [[ ! -f "$out_vcf" ]]; then
      gatk SelectVariants \
        -V "$vcf" \
        --exclude-filtered \
        -O "$out_vcf"
      tabix -f -p vcf "$out_vcf"
    else
      echo " [SKIP] PASS VCF exists for $sample"
    fi

    # 2) Export table
    out_tsv="${tsv_dir}/${sample}.noPoN.tsv"
    if [[ ! -f "$out_tsv" ]]; then
      gatk VariantsToTable \
        -V "$out_vcf" \
        -O "$out_tsv" \
        ${INFO_FIELDS[@]/#/-F } \
        ${FORMAT_FIELDS[@]/#/-GF }
    else
      echo " [SKIP] TSV exists for $sample"
    fi
  done
}

# Process samples and dilutions
process_set "$VCF_SAMPLES_DIR" "$SELECT_SAMPLES_DIR" "$TSV_SAMPLES_DIR"
process_set "$VCF_DIL_DIR" "$SELECT_DIL_DIR" "$TSV_DIL_DIR"

echo "Extraction complete. TSVs in $TSV_SAMPLES_DIR and $TSV_DIL_DIR"
