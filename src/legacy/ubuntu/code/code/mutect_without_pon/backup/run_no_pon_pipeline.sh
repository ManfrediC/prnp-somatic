#!/usr/bin/env bash
# ============================================================
# Resumable pipeline: Mutect2 (no PoN) -> orientation model
#   -> FilterMutectCalls -> bcftools normalise -> Funcotator
#   CJDs + Ctrls   AND   Dilution series.
# Each stage skips samples whose expected output already exists.
# Missing Mutect2 .stats files are optional: if absent, we run
# FilterMutectCalls without the --stats argument.
# ============================================================
set -euo pipefail
shopt -s nullglob

# ---------------------- resources ---------------------------
dbDir="/home/mcarta/databases"
fasta="$dbDir/chr2_chr4_chr20.fasta"
gnomAD="$dbDir/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
funcDS="/add/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38"

# ---------------------- input sets --------------------------
bamSamples="/add/seq_data/2025-02-07_samples_bwa"      # CJD + Ctrl
bamDilutions="/add/seq_data/2025-02-07_dilutions_bwa"  # dilution series

dilutionNames=( \
  "NA100_undil" "NA100_1to10" "NA99A1_undil" "A100_1to2" \
  "NA99A1_1to5" "NA995A05_undil" "NA100_1to2" )

# ---------------------- output roots ------------------------
rootDir="/add/results/no_PoN"
mkdir -p "$rootDir"{/mutect_raw,/mutect_raw_dil,\
/orientation,/orientation_dil,\
/filtered,/filtered_dil,\
/norm,/norm_dil,\
/annot,/annot/dil}

# ============================================================
# Stage 0 : BAM indexing
# ============================================================
echo "=== Stage 0: indexing BAMs if needed ==="
for bam in "$bamSamples"/CJD*.bwa.picard.markedDup.recal.bam \
           "$bamSamples"/Ctrl*.bwa.picard.markedDup.recal.bam \
           "$bamDilutions"/*.bwa.picard.markedDup.recal.bam; do
  [[ -f "${bam}.bai" ]] && continue
  echo "Indexing $(basename "$bam")"; samtools index -@4 "$bam"
done

# ============================================================
# Stage 1 : Mutect2 tumour-only  (no PoN)
# ============================================================
echo "=== Stage 1: Mutect2 tumour-only (no PoN) ==="

for bam in "$bamSamples"/CJD*.bwa.picard.markedDup.recal.bam \
           "$bamSamples"/Ctrl*.bwa.picard.markedDup.recal.bam; do
  [[ ! -f "$bam" ]] && continue
  sample=$(basename "$bam" .bwa.picard.markedDup.recal.bam)
  outVcf="$rootDir/mutect_raw/${sample}.raw.vcf.gz"
  [[ -f "$outVcf" ]] && { echo "[Stage1] SKIP $sample (raw VCF exists)"; continue; }
  gatk Mutect2 -R "$fasta" -I "$bam" -tumor-sample "$sample" \
      --germline-resource "$gnomAD" \
      --af-of-alleles-not-in-resource 0.0000025 \
      --f1r2-tar-gz "$rootDir/mutect_raw/${sample}.f1r2.tar.gz" \
      -O "$outVcf"
done

# Dilutions
for name in "${dilutionNames[@]}"; do
  bam="$bamDilutions/${name}.bwa.picard.markedDup.recal.bam"
  [[ ! -f "$bam" ]] && { echo "[Stage1] WARN missing $bam"; continue; }
  outVcf="$rootDir/mutect_raw_dil/${name}.raw.vcf.gz"
  [[ -f "$outVcf" ]] && { echo "[Stage1] SKIP $name"; continue; }
  gatk Mutect2 -R "$fasta" -I "$bam" -tumor-sample "$name" \
      --germline-resource "$gnomAD" \
      --af-of-alleles-not-in-resource 0.0000025 \
      --f1r2-tar-gz "$rootDir/mutect_raw_dil/${name}.f1r2.tar.gz" \
      -O "$outVcf"
done

# ============================================================
# Stage 2 : LearnReadOrientationModel
# ============================================================
echo "=== Stage 2: LearnReadOrientationModel ==="
for tar in "$rootDir"/mutect_raw/*.f1r2.tar.gz; do
  sample=$(basename "$tar" .f1r2.tar.gz)
  out="$rootDir/orientation/${sample}_orientation.tar.gz"
  [[ -f "$out" ]] && { echo "[Stage2] SKIP $sample"; continue; }
  gatk LearnReadOrientationModel -I "$tar" -O "$out" --num-em-iterations 50
done
for tar in "$rootDir"/mutect_raw_dil/*.f1r2.tar.gz; do
  sample=$(basename "$tar" .f1r2.tar.gz)
  out="$rootDir/orientation_dil/${sample}_orientation.tar.gz"
  [[ -f "$out" ]] && { echo "[Stage2] SKIP $sample (dilution)"; continue; }
  gatk LearnReadOrientationModel -I "$tar" -O "$out" --num-em-iterations 50
done

# ============================================================
# Stage 3 : FilterMutectCalls (orientation priors, stats optional)
# ============================================================
echo "=== Stage 3: FilterMutectCalls ==="
filter_one() {
  local vcf=$1 orient=$2 out=$3
  local stats="${vcf}.stats"
  if [[ -f "$out" ]]; then echo "[Stage3] SKIP $(basename "$vcf")"; return; fi
  if [[ -f "$stats" ]]; then
    gatk FilterMutectCalls -R "$fasta" -V "$vcf" --orientation-bias-artifact-priors "$orient" --stats "$stats" -O "$out"
  else
    echo "[Stage3] stats missing for $(basename "$vcf"), running without .stats"
    gatk FilterMutectCalls -R "$fasta" -V "$vcf" --orientation-bias-artifact-priors "$orient" -O "$out"
  fi
}
for vcf in "$rootDir"/mutect_raw/*.raw.vcf.gz; do
  sample=$(basename "$vcf" .raw.vcf.gz)
  orient="$rootDir/orientation/${sample}_orientation.tar.gz"
  filter_one "$vcf" "$orient" "$rootDir/filtered/${sample}.filtered.vcf.gz"
done
for vcf in "$rootDir"/mutect_raw_dil/*.raw.vcf.gz; do
  sample=$(basename "$vcf" .raw.vcf.gz)
  orient="$rootDir/orientation_dil/${sample}_orientation.tar.gz"
  filter_one "$vcf" "$orient" "$rootDir/filtered_dil/${sample}.filtered.vcf.gz"
done

# ============================================================
# Stage 4 : bcftools normalise (split + left-align)
# ============================================================
echo "=== Stage 4: bcftools normalise ==="
for vcf in "$rootDir"/filtered/*.filtered.vcf.gz; do
  sample=$(basename "$vcf" .filtered.vcf.gz)
  out="$rootDir/norm/${sample}.norm.vcf.gz"
  [[ -f "$out" ]] && { echo "[Stage4] SKIP $sample"; continue; }
  bcftools norm -m-any -f "$fasta" "$vcf" -Oz -o "$out"; tabix -f -p vcf "$out"
done
for vcf in "$rootDir"/filtered_dil/*.filtered.vcf.gz; do
  sample=$(basename "$vcf" .filtered.vcf.gz)
  out="$rootDir/norm_dil/${sample}.norm.vcf.gz"
  [[ -f "$out" ]] && { echo "[Stage4] SKIP $sample (dilution)"; continue; }
  bcftools norm -m-any -f "$fasta" "$vcf" -Oz -o "$out"; tabix -f -p vcf "$out"
done

# ============================================================
# Stage 5 : Funcotator annotation
# ============================================================
echo "=== Stage 5: Funcotator annotation ==="
annot_one() {
  local vcf=$1 out=$2
  [[ -f "$out" ]] && { echo "[Stage5] SKIP $(basename "$vcf")"; return; }
  gatk Funcotator --variant "$vcf" --reference "$fasta" --data-sources-path "$funcDS" --ref-version hg38 --output "$out" --output-file-format VCF
  tabix -f -p vcf "$out"
}
for vcf in "$rootDir"/norm/*.norm.vcf.gz; do
  sample=$(basename "$vcf" .norm.vcf.gz)
  annot_one "$vcf" "$rootDir/annot/${sample}.func.vcf.gz"
done
for vcf in "$rootDir"/norm_dil/*.norm.vcf.gz; do
  sample=$(basename "$vcf" .norm.vcf.gz)
  annot_one "$vcf" "$rootDir/annot/dil/${sample}.func.vcf.gz"
done

echo "=== Pipeline finished (resume mode) ==="
