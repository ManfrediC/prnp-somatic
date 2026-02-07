#!/usr/bin/env bash
# ============================================================
# Resumable no-PoN pipeline: Mutect2 -> orientation model
#   -> FilterMutectCalls (PASS only) -> VariantAnnotator
#   -> bcftools normalise -> Funcotator
# Applies to CJD+Ctrl samples AND dilution series.
# Each stage skips completed outputs.
# ============================================================
set -euo pipefail
shopt -s nullglob

# ---------------------- resources ---------------------------
DB_DIR=/home/mcarta/databases
FASTA=$DB_DIR/chr2_chr4_chr20.fasta
INTERVALS=$DB_DIR/capture_targets.interval_list      # restrict to probe targets
GNOMAD=$DB_DIR/somatic-hg38_af-only-gnomad.hg38.vcf.gz
FUNC_DS=/add/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38

# ---------------------- inputs ------------------------------
BAMS_SAMPLES=/add/seq_data/2025-02-07_samples_bwa
BAMS_DIL=/add/seq_data/2025-02-07_dilutions_bwa
DIL_NAMES=(NA100_undil NA100_1to10 NA99A1_undil A100_1to2 NA99A1_1to5 NA995A05_undil NA100_1to2)

# ---------------------- output dirs ------------------------
ROOT=/add/results/no_PoN
for d in mutect_raw mutect_raw_dil orientation orientation_dil filtered filtered_dil scores norm norm_dil annot annot/dil; do
    mkdir -p $ROOT/$d
    [[ $d == scores ]] && mkdir -p $ROOT/scores_dil
done

# ------------------------------------------------------------
# Stage 0: BAM indexing
# ------------------------------------------------------------
echo "=== Stage 0: BAM indexing ==="
for bam in $BAMS_SAMPLES/CJD*.bwa.picard.markedDup.recal.bam \
           $BAMS_SAMPLES/Ctrl*.bwa.picard.markedDup.recal.bam \
           $BAMS_DIL/*.bwa.picard.markedDup.recal.bam; do
  [[ -f ${bam}.bai ]] && continue
  echo "Indexing $(basename $bam)"
  samtools index -@4 $bam
done

# ------------------------------------------------------------
# Stage 1: Mutect2 tumour-only (no PoN)
# ------------------------------------------------------------
echo "=== Stage 1: Mutect2 ==="
process_mutect() {
  local bam=$1 sample=$2 outdir=$3
  local outvcf=$ROOT/$outdir/${sample}.raw.vcf.gz
  [[ -f $outvcf ]] && { echo "[Stage1] SKIP $sample"; return; }
  gatk Mutect2 \
       -R $FASTA \
       -I $bam \
       -tumor-sample $sample \
       --germline-resource $GNOMAD \
       --af-of-alleles-not-in-resource 0.0000025 \
       --f1r2-tar-gz $ROOT/$outdir/${sample}.f1r2.tar.gz \
       --intervals $INTERVALS \
       --annotation RMSMappingQuality \
       -O $outvcf
}
for bam in $BAMS_SAMPLES/CJD*.bwa.picard.markedDup.recal.bam \
           $BAMS_SAMPLES/Ctrl*.bwa.picard.markedDup.recal.bam; do
  [[ ! -f $bam ]] && continue
  sample=$(basename $bam .bwa.picard.markedDup.recal.bam)
  process_mutect $bam $sample mutect_raw
done
for name in "${DIL_NAMES[@]}"; do
  bam=$BAMS_DIL/${name}.bwa.picard.markedDup.recal.bam
  process_mutect $bam $name mutect_raw_dil
done

# ------------------------------------------------------------
# Stage 2: LearnReadOrientationModel
# ------------------------------------------------------------
echo "=== Stage 2: LearnReadOrientationModel ==="
process_orient() {
  local tar=$1 outdir=$2
  local sample=$(basename $tar .f1r2.tar.gz)
  local out=$ROOT/$outdir/${sample}_orientation.tar.gz
  [[ -f $out ]] && { echo "[Stage2] SKIP $sample"; return; }
  gatk LearnReadOrientationModel -I $tar -O $out --num-em-iterations 50
}
for tar in $ROOT/mutect_raw/*.f1r2.tar.gz; do process_orient $tar orientation; done
for tar in $ROOT/mutect_raw_dil/*.f1r2.tar.gz; do process_orient $tar orientation_dil; done

# ------------------------------------------------------------
# Stage 3: FilterMutectCalls + PASS filtering
# ------------------------------------------------------------
echo "=== Stage 3: FilterMutectCalls (PASS only) ==="
filter_one() {
  local vcf=$1 orient=$2 outdir=$3
  local sample=$(basename $vcf .raw.vcf.gz)
  local tmp=$ROOT/$outdir/${sample}.tmp.vcf.gz
  local out=$ROOT/$outdir/${sample}.filtered.vcf.gz
  [[ -f $out ]] && { echo "[Stage3] SKIP $sample"; return; }
  local stats=${vcf}.stats
  if [[ -f $stats ]]; then
    gatk FilterMutectCalls -R $FASTA -V $vcf \
         --orientation-bias-artifact-priors $orient \
         --stats $stats -O $tmp
  else
    echo "[Stage3] WARN stats missing for $sample"
    gatk FilterMutectCalls -R $FASTA -V $vcf \
         --orientation-bias-artifact-priors $orient -O $tmp
  fi
  # keep only PASS variants
  bcftools view -f PASS $tmp -Oz -o $out
  tabix -f -p vcf $out
  rm -f $tmp
}
for vcf in $ROOT/mutect_raw/*.raw.vcf.gz; do
  orient=$ROOT/orientation/$(basename $vcf .raw.vcf.gz)_orientation.tar.gz
  filter_one $vcf $orient filtered
done
for vcf in $ROOT/mutect_raw_dil/*.raw.vcf.gz; do
  orient=$ROOT/orientation_dil/$(basename $vcf .raw.vcf.gz)_orientation.tar.gz
  filter_one $vcf $orient filtered_dil
done

# ------------------------------------------------------------
# Stage 4: VariantAnnotator (QUAL, MQ, QD, BaseQualityRankSumTest, MappingQualityRankSumTest)
# ------------------------------------------------------------
echo "=== Stage 4: VariantAnnotator ==="
annotate_scores() {
  local vcf=$1 outdir=$2
  local sample=$(basename $vcf .filtered.vcf.gz)
  local out=$ROOT/$outdir/${sample}.scores.vcf.gz
  [[ -f $out ]] && { echo "[Stage4] SKIP $sample"; return; }
  gatk VariantAnnotator \
       -R $FASTA \
       -V $vcf \
       -O $out \
       --annotation QualByDepth \
       --annotation MappingQuality \
       --annotation BaseQualityRankSumTest \
       --annotation MappingQualityRankSumTest
  tabix -f -p vcf $out
}
for vcf in $ROOT/filtered/*.filtered.vcf.gz; do annotate_scores $vcf scores; done
for vcf in $ROOT/filtered_dil/*.filtered.vcf.gz; do annotate_scores $vcf scores_dil; done

# ------------------------------------------------------------
# Stage 5: bcftools normalise (split + left-align)
# ------------------------------------------------------------
echo "=== Stage 5: bcftools normalise ==="
for vcf in $ROOT/filtered/*.filtered.vcf.gz; do
  sample=$(basename $vcf .filtered.vcf.gz)
  out=$ROOT/norm/${sample}.norm.vcf.gz
  [[ -f $out ]] && { echo "[Stage5] SKIP $sample"; continue; }
  bcftools norm -m-any -f $FASTA $vcf -Oz -o $out; tabix -f -p vcf $out
done
for vcf in $ROOT/filtered_dil/*.filtered.vcf.gz; do
  sample=$(basename $vcf .filtered.vcf.gz)
  out=$ROOT/norm_dil/${sample}.norm.vcf.gz
  [[ -f $out ]] && { echo "[Stage5] SKIP $sample (dil)"; continue; }
  bcftools norm -m-any -f $FASTA $vcf -Oz -o $out; tabix -f -p vcf $out
done

# ------------------------------------------------------------
# Stage 6: Funcotator annotation
# ------------------------------------------------------------
echo "=== Stage 6: Funcotator annotation ==="
funcotate() {
  local vcf=$1 outdir=$2
  local sample=$(basename $vcf .norm.vcf.gz)
  local out=$ROOT/$outdir/${sample}.func.vcf.gz
  [[ -f $out ]] && { echo "[Stage6] SKIP $sample"; return; }
  gatk Funcotator --variant $vcf --reference $FASTA \
                  --data-sources-path $FUNC_DS --ref-version hg38 \
                  --output $out --output-file-format VCF
  tabix -f -p vcf $out
}
for vcf in $ROOT/norm/*.norm.vcf.gz; do funcotate $vcf annot; done
for vcf in $ROOT/norm_dil/*.norm.vcf.gz; do funcotate $vcf annot/dil; done

echo "=== Pipeline finished ==="
