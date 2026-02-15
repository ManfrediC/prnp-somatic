#!/usr/bin/env bash
# ============================================================
# Resumable pipeline WITH PoN: Mutect2 -> orientation model
#   -> FilterMutectCalls (PASS only) -> VariantAnnotator
#   -> bcftools normalise -> Funcotator
# Applies to **CJD tumour samples only** (no controls, no dilutions).
# Each stage skips completed outputs.
# ============================================================
set -euo pipefail
shopt -s nullglob

# ---------------------- resources ---------------------------
DB_DIR=/home/mcarta/databases
FASTA=$DB_DIR/chr2_chr4_chr20.fasta
INTERVALS=$DB_DIR/capture_targets.interval_list           # 3-gene target space
GNOMAD=$DB_DIR/somatic-hg38_af-only-gnomad.hg38.vcf.gz

# Panel of Normals built from controls
PON=/add/results/PoN/panel_of_normals/CJD_controls_PoN.vcf.gz
[[ -f ${PON}.idx ]] || gatk IndexFeatureFile -I "$PON"

FUNC_DS=/add/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38

# ---------------------- inputs ------------------------------
BAMS_SAMPLES=/add/seq_data/2025-02-07_samples_bwa   # CJD BAMs only

# ---------------------- output dirs ------------------------
ROOT=/add/results/with_PoN
for d in mutect_raw orientation filtered scores norm annot; do
    mkdir -p "$ROOT/$d"
done

# ------------------------------------------------------------
# Stage 0: BAM indexing
# ------------------------------------------------------------
echo "=== Stage 0: BAM indexing ==="
for bam in "$BAMS_SAMPLES"/CJD*.bwa.picard.markedDup.recal.bam; do
  [[ -f ${bam}.bai ]] && continue
  echo "Indexing $(basename "$bam")"
  samtools index -@4 "$bam"
done

# ------------------------------------------------------------
# Stage 1: Mutect2 (with PoN)
# ------------------------------------------------------------
echo "=== Stage 1: Mutect2 ==="
process_mutect() {
  local bam=$1 sample=$2 outdir=$3
  local outvcf=$ROOT/$outdir/${sample}.raw.vcf.gz
  [[ -f $outvcf ]] && { echo "[Stage1] SKIP $sample"; return; }
  gatk Mutect2 \
       -R "$FASTA" \
       -I "$bam" \
       -tumor-sample "$sample" \
       --panel-of-normals "$PON" \
       --germline-resource "$GNOMAD" \
       --af-of-alleles-not-in-resource 0.0000025 \
       --f1r2-tar-gz "$ROOT/$outdir/${sample}.f1r2.tar.gz" \
       --intervals "$INTERVALS" \
       --annotation RMSMappingQuality \
       -O "$outvcf"
}
for bam in "$BAMS_SAMPLES"/CJD*.bwa.picard.markedDup.recal.bam; do
  sample=$(basename "$bam" .bwa.picard.markedDup.recal.bam)
  process_mutect "$bam" "$sample" mutect_raw
done

# ------------------------------------------------------------
# Stage 2: LearnReadOrientationModel
# ------------------------------------------------------------
echo "=== Stage 2: LearnReadOrientationModel ==="
process_orient() {
  local tar=$1
  local sample
  sample=$(basename "$tar" .f1r2.tar.gz)
  local out=$ROOT/orientation/${sample}_orientation.tar.gz
  [[ -f $out ]] && { echo "[Stage2] SKIP $sample"; return; }
  gatk LearnReadOrientationModel -I "$tar" -O "$out" --num-em-iterations 50
}
for tar in "$ROOT"/mutect_raw/*.f1r2.tar.gz; do
  process_orient "$tar"
done

# ------------------------------------------------------------
# Stage 3: FilterMutectCalls + PASS filtering
# ------------------------------------------------------------
echo "=== Stage 3: FilterMutectCalls (PASS only) ==="
filter_one() {
  local vcf=$1 orient=$2
  local sample
  sample=$(basename "$vcf" .raw.vcf.gz)
  local tmp=$ROOT/filtered/${sample}.tmp.vcf.gz
  local out=$ROOT/filtered/${sample}.filtered.vcf.gz
  [[ -f $out ]] && { echo "[Stage3] SKIP $sample"; return; }
  local stats=${vcf}.stats
  if [[ -f $stats ]]; then
    gatk FilterMutectCalls -R "$FASTA" -V "$vcf" \
         --orientation-bias-artifact-priors "$orient" \
         --stats "$stats" -O "$tmp"
  else
    echo "[Stage3] WARN stats missing for $sample"
    gatk FilterMutectCalls -R "$FASTA" -V "$vcf" \
         --orientation-bias-artifact-priors "$orient" -O "$tmp"
  fi
  bcftools view -f PASS "$tmp" -Oz -o "$out"
  tabix -f -p vcf "$out"
  rm -f "$tmp"
}
for vcf in "$ROOT"/mutect_raw/*.raw.vcf.gz; do
  orient=$ROOT/orientation/$(basename "$vcf" .raw.vcf.gz)_orientation.tar.gz
  filter_one "$vcf" "$orient"
done

# ------------------------------------------------------------
# Stage 4: VariantAnnotator (QUAL, MQ, QD, RankSums)
# ------------------------------------------------------------
echo "=== Stage 4: VariantAnnotator ==="
annotate_scores() {
  local vcf=$1
  local sample
  sample=$(basename "$vcf" .filtered.vcf.gz)
  local out=$ROOT/scores/${sample}.scores.vcf.gz
  [[ -f $out ]] && { echo "[Stage4] SKIP $sample"; return; }
  gatk VariantAnnotator \
       -R "$FASTA" \
       -V "$vcf" \
       --intervals "$INTERVALS" \
       -O "$out" \
       --annotation QualByDepth \
       --annotation MappingQuality \
       --annotation BaseQualityRankSumTest \
       --annotation MappingQualityRankSumTest
  tabix -f -p vcf "$out"
}
for vcf in "$ROOT"/filtered/*.filtered.vcf.gz; do
  annotate_scores "$vcf"
done

# ------------------------------------------------------------
# Stage 5: bcftools normalise (split + left-align)
# ------------------------------------------------------------
echo "=== Stage 5: bcftools normalise ==="
for vcf in "$ROOT"/filtered/*.filtered.vcf.gz; do
  sample=$(basename "$vcf" .filtered.vcf.gz)
  out=$ROOT/norm/${sample}.norm.vcf.gz
  [[ -f $out ]] && { echo "[Stage5] SKIP $sample"; continue; }
  bcftools norm -m-any -f "$FASTA" "$vcf" -Oz -o "$out"
  tabix -f -p vcf "$out"
done

# ------------------------------------------------------------
# Stage 6: Funcotator annotation
# ------------------------------------------------------------
echo "=== Stage 6: Funcotator annotation ==="
funcotate() {
  local vcf=$1
  local sample
  sample=$(basename "$vcf" .norm.vcf.gz)
  local out=$ROOT/annot/${sample}.func.vcf.gz
  [[ -f $out ]] && { echo "[Stage6] SKIP $sample"; return; }
  gatk Funcotator --variant "$vcf" --reference "$FASTA" \
                  --data-sources-path "$FUNC_DS" --ref-version hg38 \
                  --intervals "$INTERVALS" \
                  --output "$out" --output-file-format VCF
  tabix -f -p vcf "$out"
}
for vcf in "$ROOT"/norm/*.norm.vcf.gz; do
  funcotate "$vcf"
done

echo "=== Pipeline finished ==="
