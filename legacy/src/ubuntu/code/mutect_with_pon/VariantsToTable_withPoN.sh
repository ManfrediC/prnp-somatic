#!/usr/bin/env bash
# ============================================================
# Extract PASS variants to TSV for CJD samples processed with the PoN pipeline
# Uses GATK SelectVariants and VariantsToTable
# ============================================================
set -euo pipefail

# Directories
VCF_SAMPLES_DIR="/add/results/with_PoN/annot"              # Funcotator outputs for samples
SELECT_SAMPLES_DIR="/add/results/with_PoN/select_variants/samples"
TSV_SAMPLES_DIR="/add/results/with_PoN/variant_tables/samples"

# Create output directories
mkdir -p "$SELECT_SAMPLES_DIR" "$TSV_SAMPLES_DIR"

# Common fields
INFO_FIELDS=(CHROM POS REF ALT FILTER QUAL MQ QD AF GNOMAD_AF FS BaseQualityRankSumTest MappingQualityRankSumTest AS_SB_TABLE STRANDQ FUNCOTATION)
FORMAT_FIELDS=(GT DP AD F1R2 F2R1 SB)

# Function to process sample VCFs
process_samples() {
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
    out_tsv="${tsv_dir}/${sample}.withPoN.tsv"
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

# Run only for samples
process_samples "$VCF_SAMPLES_DIR" "$SELECT_SAMPLES_DIR" "$TSV_SAMPLES_DIR"

echo "Extraction complete. TSVs saved to $TSV_SAMPLES_DIR"
