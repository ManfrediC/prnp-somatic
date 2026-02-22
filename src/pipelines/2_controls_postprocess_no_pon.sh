#!/usr/bin/env bash
# ============================================================
# Controls-only post-processing (no PoN), starting after Mutect2:
#   Stage 2: LearnReadOrientationModel
#   Stage 3: FilterMutectCalls + PASS filtering
#   Stage 4: VariantAnnotator (QUAL/MQ/QD/RankSum annotations)
#   Stage 5: bcftools normalise (split + left-align)
#   Stage 6: Funcotator annotation
#   Stage 7: bcftools annotate gnomAD AF
# ============================================================
set -euo pipefail
shopt -s nullglob

# Keep caller-provided values so ENV_FILE does not silently override them.
CLI_DRY_RUN="${DRY_RUN-}"
CLI_JAVA_MEM_GB="${JAVA_MEM_GB-}"

# -----------------------
# Repo root discovery (works from anywhere)
# -----------------------
find_repo_root() {
  local start d
  start="${1:-}"
  [[ -n "$start" ]] || { echo "ERROR: find_repo_root requires a start directory" >&2; return 2; }
  d="$start"
  while [[ "$d" != "/" ]]; do
    [[ -f "$d/Makefile" ]] && { echo "$d"; return 0; }
    d="$(dirname "$d")"
  done
  echo "ERROR: could not find repo root (Makefile not found)" >&2
  return 1
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(find_repo_root "$(cd "$SCRIPT_DIR/../.." && pwd)")"
cd "$REPO_ROOT"

# -----------------------
# Optional env file
# -----------------------
ENV_FILE_DEFAULT="$REPO_ROOT/config/preprocessing.env"
ENV_FILE="${ENV_FILE:-$ENV_FILE_DEFAULT}"
if [[ -f "$ENV_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$ENV_FILE"
fi

# -----------------------
# Defaults (override via env or ENV_FILE)
# -----------------------
DRY_RUN="${CLI_DRY_RUN:-${DRY_RUN:-0}}"
JAVA_MEM_GB="${CLI_JAVA_MEM_GB:-${JAVA_MEM_GB:-8}}"

# -----------------------
# Path defaults (repo-relative)
# -----------------------
MUTECT2_CONTROLS_OUT_ROOT="${MUTECT2_CONTROLS_OUT_ROOT:-runs/mutect2_controls_no_pon}"
if [[ "$MUTECT2_CONTROLS_OUT_ROOT" == run/* ]]; then
  MUTECT2_CONTROLS_OUT_ROOT="runs/${MUTECT2_CONTROLS_OUT_ROOT#run/}"
fi

OUT_ROOT="${OUT_ROOT:-$MUTECT2_CONTROLS_OUT_ROOT}"

RAW_DIR="${RAW_DIR:-$OUT_ROOT/vcf}"
F1R2_DIR="${F1R2_DIR:-$OUT_ROOT/f1r2}"
ORIENT_DIR="${ORIENT_DIR:-$OUT_ROOT/orientation}"
FILTERED_DIR="${FILTERED_DIR:-$OUT_ROOT/filtered}"
SCORES_DIR="${SCORES_DIR:-$OUT_ROOT/scores}"
NORM_DIR="${NORM_DIR:-$OUT_ROOT/norm}"
ANNOT_DIR="${ANNOT_DIR:-$OUT_ROOT/annot}"
ANNOT_GNOMAD_DIR="${ANNOT_GNOMAD_DIR:-$OUT_ROOT/annot_with_gnomad}"

REF_FASTA="${REF_FASTA:-resources/chr2_chr4_chr20.fasta}"
FUNCOTATOR_DS="${FUNCOTATOR_DS:-resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38}"
GNOMAD_AF_VCF="${GNOMAD_AF_VCF:-resources/somatic-hg38_af-only-gnomad.hg38.vcf.gz}"

to_abs() {
  local p="$1"
  # Allow both repo-relative and absolute paths in ENV overrides.
  if [[ "$p" = /* ]]; then
    echo "$p"
  else
    echo "$REPO_ROOT/$p"
  fi
}

RAW_DIR="$(to_abs "$RAW_DIR")"
F1R2_DIR="$(to_abs "$F1R2_DIR")"
ORIENT_DIR="$(to_abs "$ORIENT_DIR")"
FILTERED_DIR="$(to_abs "$FILTERED_DIR")"
SCORES_DIR="$(to_abs "$SCORES_DIR")"
NORM_DIR="$(to_abs "$NORM_DIR")"
ANNOT_DIR="$(to_abs "$ANNOT_DIR")"
ANNOT_GNOMAD_DIR="$(to_abs "$ANNOT_GNOMAD_DIR")"
REF_FASTA="$(to_abs "$REF_FASTA")"
FUNCOTATOR_DS="$(to_abs "$FUNCOTATOR_DS")"
GNOMAD_AF_VCF="$(to_abs "$GNOMAD_AF_VCF")"

mkdir -p "$ORIENT_DIR" "$FILTERED_DIR" "$SCORES_DIR" "$NORM_DIR" "$ANNOT_DIR" "$ANNOT_GNOMAD_DIR"

# -----------------------
# Tiny helpers
# -----------------------
have() { command -v "$1" >/dev/null 2>&1; }
die() { echo "ERROR: $*" >&2; exit 1; }

run() {
  echo "+ $*"
  if [[ "$DRY_RUN" == "1" ]]; then
    return 0
  fi
  "$@"
}

# -----------------------
# Tool + file checks
# -----------------------
have gatk     || die "gatk not in PATH (did you activate conda?)"
have bcftools || die "bcftools not in PATH (did you activate conda?)"
have tabix    || die "tabix not in PATH (did you activate conda?)"

[[ -s "$REF_FASTA" ]] || die "Missing REF_FASTA: $REF_FASTA"
[[ -s "$REF_FASTA.fai" ]] || die "Missing FASTA index: $REF_FASTA.fai"
[[ -d "$FUNCOTATOR_DS" ]] || die "Missing Funcotator data sources dir: $FUNCOTATOR_DS"
[[ -s "$GNOMAD_AF_VCF" ]] || die "Missing GNOMAD_AF_VCF: $GNOMAD_AF_VCF"
[[ -s "$GNOMAD_AF_VCF.tbi" || -s "$GNOMAD_AF_VCF.csi" ]] || die "Missing index for GNOMAD_AF_VCF: $GNOMAD_AF_VCF.tbi/.csi"

f1r2_tars=( "$F1R2_DIR"/*.f1r2.tar.gz )
raw_vcfs=( "$RAW_DIR"/*.raw.vcf.gz )
[[ "${#f1r2_tars[@]}" -gt 0 ]] || die "No F1R2 tarballs found in: $F1R2_DIR"
[[ "${#raw_vcfs[@]}" -gt 0 ]] || die "No raw VCFs found in: $RAW_DIR"
raw_samples=()
for vcf in "${raw_vcfs[@]}"; do
  # Track Stage-1 sample IDs once so the final completeness check is deterministic.
  # This avoids coupling final checks to whichever intermediate glob happens to exist.
  raw_samples+=( "$(basename "$vcf" .raw.vcf.gz)" )
done

echo "== Controls post-processing (no PoN) =="
echo "Repo root:      $REPO_ROOT"
echo "RAW_DIR:        $RAW_DIR"
echo "F1R2_DIR:       $F1R2_DIR"
echo "ORIENT_DIR:     $ORIENT_DIR"
echo "FILTERED_DIR:   $FILTERED_DIR"
echo "SCORES_DIR:     $SCORES_DIR"
echo "NORM_DIR:       $NORM_DIR"
echo "ANNOT_DIR:      $ANNOT_DIR"
echo "ANNOT_GNOMAD_DIR: $ANNOT_GNOMAD_DIR"
echo "REF_FASTA:      $REF_FASTA"
echo "FUNCOTATOR_DS:  $FUNCOTATOR_DS"
echo "GNOMAD_AF_VCF:  $GNOMAD_AF_VCF"
echo "JAVA_MEM_GB:    $JAVA_MEM_GB"
echo

# ------------------------------------------------------------
# Stage 2: LearnReadOrientationModel
# ------------------------------------------------------------
echo "=== Stage 2: LearnReadOrientationModel ==="
for tar in "${f1r2_tars[@]}"; do
  sample="$(basename "$tar" .f1r2.tar.gz)"
  out="$ORIENT_DIR/${sample}_orientation.tar.gz"
  # Orientation priors are consumed by FilterMutectCalls in Stage 3.
  [[ -s "$out" ]] && { echo "[Stage2] SKIP $sample"; continue; }
  run gatk --java-options "-Xmx${JAVA_MEM_GB}g" LearnReadOrientationModel \
    -I "$tar" \
    -O "$out" \
    --num-em-iterations 50
done

# ------------------------------------------------------------
# Stage 3: FilterMutectCalls + PASS filtering
# ------------------------------------------------------------
echo "=== Stage 3: FilterMutectCalls (PASS only) ==="
for vcf in "${raw_vcfs[@]}"; do
  sample="$(basename "$vcf" .raw.vcf.gz)"
  orient="$ORIENT_DIR/${sample}_orientation.tar.gz"
  tmp="$FILTERED_DIR/${sample}.tmp.vcf.gz"
  out="$FILTERED_DIR/${sample}.filtered.vcf.gz"
  stats="${vcf}.stats"

  [[ -s "$out" ]] && { echo "[Stage3] SKIP $sample"; continue; }
  [[ -s "$orient" ]] || die "Stage 3 requires orientation file for $sample: $orient"
  [[ -s "$stats" ]] || die "Stage 3 requires stats file for $sample: $stats"
  run gatk --java-options "-Xmx${JAVA_MEM_GB}g" FilterMutectCalls \
    -R "$REF_FASTA" \
    -V "$vcf" \
    --orientation-bias-artifact-priors "$orient" \
    --stats "$stats" \
    -O "$tmp"

  run bcftools view -f PASS "$tmp" -Oz -o "$out"
  run tabix -f -p vcf "$out"
  if [[ "$DRY_RUN" == "0" ]]; then
    # Keep only final PASS-filtered file to avoid stale intermediates.
    rm -f "$tmp"
  fi
done

# ------------------------------------------------------------
# Stage 4: VariantAnnotator (QUAL, MQ, QD, RankSum)
# ------------------------------------------------------------
echo "=== Stage 4: VariantAnnotator ==="
filtered_vcfs=( "$FILTERED_DIR"/*.filtered.vcf.gz )
[[ "${#filtered_vcfs[@]}" -gt 0 ]] || die "Stage 4 has no filtered VCFs in: $FILTERED_DIR"
for vcf in "${filtered_vcfs[@]}"; do
  sample="$(basename "$vcf" .filtered.vcf.gz)"
  out="$SCORES_DIR/${sample}.scores.vcf.gz"
  [[ -s "$out" ]] && { echo "[Stage4] SKIP $sample"; continue; }
  run gatk --java-options "-Xmx${JAVA_MEM_GB}g" VariantAnnotator \
    -R "$REF_FASTA" \
    -V "$vcf" \
    -O "$out" \
    --annotation QualByDepth \
    --annotation MappingQuality \
    --annotation BaseQualityRankSumTest \
    --annotation MappingQualityRankSumTest
  run tabix -f -p vcf "$out"
done

# ------------------------------------------------------------
# Stage 5: bcftools normalise (split + left-align)
# ------------------------------------------------------------
echo "=== Stage 5: bcftools normalise ==="
scores_vcfs=( "$SCORES_DIR"/*.scores.vcf.gz )
[[ "${#scores_vcfs[@]}" -gt 0 ]] || die "Stage 5 has no scored VCFs in: $SCORES_DIR"
for vcf in "${scores_vcfs[@]}"; do
  sample="$(basename "$vcf" .scores.vcf.gz)"
  out="$NORM_DIR/${sample}.norm.vcf.gz"
  [[ -s "$out" ]] && { echo "[Stage5] SKIP $sample"; continue; }
  run bcftools norm -m-any -f "$REF_FASTA" "$vcf" -Oz -o "$out"
  run tabix -f -p vcf "$out"
done

# ------------------------------------------------------------
# Stage 6: Funcotator annotation
# ------------------------------------------------------------
echo "=== Stage 6: Funcotator annotation ==="
norm_vcfs=( "$NORM_DIR"/*.norm.vcf.gz )
[[ "${#norm_vcfs[@]}" -gt 0 ]] || die "Stage 6 has no normalised VCFs in: $NORM_DIR"
for vcf in "${norm_vcfs[@]}"; do
  sample="$(basename "$vcf" .norm.vcf.gz)"
  out="$ANNOT_DIR/${sample}.func.vcf.gz"
  [[ -s "$out" ]] && { echo "[Stage6] SKIP $sample"; continue; }
  run gatk --java-options "-Xmx${JAVA_MEM_GB}g" Funcotator \
    --variant "$vcf" \
    --reference "$REF_FASTA" \
    --data-sources-path "$FUNCOTATOR_DS" \
    --ref-version hg38 \
    --output "$out" \
    --output-file-format VCF
  run tabix -f -p vcf "$out"
done

# ------------------------------------------------------------
# Stage 7: gnomAD AF annotation with bcftools annotate
# ------------------------------------------------------------
echo "=== Stage 7: gnomAD AF annotation ==="
annot_vcfs=( "$ANNOT_DIR"/*.func.vcf.gz )
[[ "${#annot_vcfs[@]}" -gt 0 ]] || die "Stage 7 has no Funcotator VCFs in: $ANNOT_DIR"
for vcf in "${annot_vcfs[@]}"; do
  sample="$(basename "$vcf" .func.vcf.gz)"
  out="$ANNOT_GNOMAD_DIR/${sample}.func.af.vcf.gz"
  # Keep AF in a dedicated INFO/GNOMAD_AF field so downstream filters can read one stable column.
  [[ -s "$out" ]] && { echo "[Stage7] SKIP $sample"; continue; }
  run bcftools annotate \
    -a "$GNOMAD_AF_VCF" \
    # Copy AF from the gnomAD VCF into a separate INFO tag.
    -c CHROM,POS,REF,ALT,INFO/GNOMAD_AF:=INFO/AF \
    -h <(echo '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="gnomAD AF from somatic-hg38_af-only-gnomad.hg38.vcf.gz">') \
    -Oz -o "$out" \
    "$vcf"
  run tabix -f -p vcf "$out"
done

missing_samples=()
for sample in "${raw_samples[@]}"; do
  [[ -s "$ANNOT_GNOMAD_DIR/${sample}.func.af.vcf.gz" ]] || missing_samples+=( "$sample" )
done
[[ "${#missing_samples[@]}" -eq 0 ]] || die "Final completeness check failed; missing gnomAD-annotated VCFs for: ${missing_samples[*]}"

echo "=== Pipeline finished ==="
