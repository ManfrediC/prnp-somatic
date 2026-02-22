#!/usr/bin/env bash
# ============================================================
# Controls-only read/base quality metrics (no PoN):
#   Stage 1: Build BED variant sites from annotated VCFs
#   Stage 2: Link recalibrated BAMs
#   Stage 3: Remove duplicate reads (samtools -F 0x400)
#   Stage 4: Collect per-site metrics with bam-readcount
#   Stage 5: Parse readcount files to TSV metrics
# ============================================================
set -euo pipefail
shopt -s nullglob

# Keep caller-provided values so ENV_FILE does not silently override them.
CLI_DRY_RUN="${DRY_RUN-}"
CLI_THREADS="${THREADS-}"

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
THREADS="${CLI_THREADS:-${THREADS:-4}}"

# -----------------------
# Path defaults (repo-relative)
# -----------------------
MUTECT2_CONTROLS_OUT_ROOT="${MUTECT2_CONTROLS_OUT_ROOT:-runs/mutect2_controls_no_pon}"
if [[ "$MUTECT2_CONTROLS_OUT_ROOT" == run/* ]]; then
  MUTECT2_CONTROLS_OUT_ROOT="runs/${MUTECT2_CONTROLS_OUT_ROOT#run/}"
fi

OUT_ROOT="${OUT_ROOT:-$MUTECT2_CONTROLS_OUT_ROOT}"

READCOUNT_QC_ROOT="${READCOUNT_QC_ROOT:-$OUT_ROOT/readcount_qc}"
READCOUNT_VCF_DIR="${READCOUNT_VCF_DIR:-$OUT_ROOT/annot_with_gnomad}"
READCOUNT_BAM_DIR="${READCOUNT_BAM_DIR:-results/final_bam}"
READCOUNT_REF_FASTA="${READCOUNT_REF_FASTA:-${REF_FASTA:-resources/chr2_chr4_chr20.fasta}}"
READCOUNT_TO_TSV_PY="${READCOUNT_TO_TSV_PY:-src/pipelines/4_readcount_to_tsv.py}"

to_abs() {
  local p="$1"
  # Allow both repo-relative and absolute paths in ENV overrides.
  if [[ "$p" = /* ]]; then
    echo "$p"
  else
    echo "$REPO_ROOT/$p"
  fi
}

READCOUNT_QC_ROOT="$(to_abs "$READCOUNT_QC_ROOT")"
READCOUNT_VCF_DIR="$(to_abs "$READCOUNT_VCF_DIR")"
READCOUNT_BAM_DIR="$(to_abs "$READCOUNT_BAM_DIR")"
READCOUNT_REF_FASTA="$(to_abs "$READCOUNT_REF_FASTA")"
READCOUNT_TO_TSV_PY="$(to_abs "$READCOUNT_TO_TSV_PY")"

BEDS_DIR="$READCOUNT_QC_ROOT/beds"
BAM_WORK_DIR="$READCOUNT_QC_ROOT/bam_work"
BAM_NODUP_DIR="$READCOUNT_QC_ROOT/bam_nodup"
READCOUNTS_DIR="$READCOUNT_QC_ROOT/readcounts"
METRICS_DIR="$READCOUNT_QC_ROOT/metrics"

mkdir -p "$BEDS_DIR" "$BAM_WORK_DIR" "$BAM_NODUP_DIR" "$READCOUNTS_DIR" "$METRICS_DIR"

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

run_to_file() {
  local out="$1"
  shift
  echo "+ $* > $out"
  if [[ "$DRY_RUN" == "1" ]]; then
    return 0
  fi
  "$@" > "$out"
}

sample_from_vcf() {
  local path base sample
  path="$1"
  base="$(basename "$path")"
  # Preferred input naming from Stage 7.
  sample="${base%.func.af.vcf.gz}"
  if [[ "$sample" == "$base" ]]; then
    # Fallback if using generic .vcf.gz annotations.
    sample="${base%.vcf.gz}"
  fi
  # Return only the sample stem; all downstream filenames are rebuilt from this value.
  echo "$sample"
}

# -----------------------
# Tool + file checks
# -----------------------
have bcftools      || die "bcftools not in PATH (did you activate conda?)"
have samtools      || die "samtools not in PATH (did you activate conda?)"
have bam-readcount || die "bam-readcount not in PATH (did you activate conda?)"
have python3       || die "python3 not in PATH"

[[ -d "$READCOUNT_VCF_DIR" ]] || die "Missing READCOUNT_VCF_DIR: $READCOUNT_VCF_DIR"
[[ -d "$READCOUNT_BAM_DIR" ]] || die "Missing READCOUNT_BAM_DIR: $READCOUNT_BAM_DIR"
[[ -s "$READCOUNT_REF_FASTA" ]] || die "Missing READCOUNT_REF_FASTA: $READCOUNT_REF_FASTA"
[[ -s "$READCOUNT_REF_FASTA.fai" ]] || die "Missing FASTA index: $READCOUNT_REF_FASTA.fai"
[[ -s "$READCOUNT_TO_TSV_PY" ]] || die "Missing parser script: $READCOUNT_TO_TSV_PY"

vcfs=( "$READCOUNT_VCF_DIR"/*.func.af.vcf.gz )
if [[ "${#vcfs[@]}" -eq 0 ]]; then
  # Backward-compatible fallback for older output naming.
  vcfs=( "$READCOUNT_VCF_DIR"/*.vcf.gz )
fi
[[ "${#vcfs[@]}" -gt 0 ]] || die "No annotated VCFs found in: $READCOUNT_VCF_DIR"

samples=()
for vcf in "${vcfs[@]}"; do
  [[ -f "$vcf" ]] || continue
  samples+=( "$(sample_from_vcf "$vcf")" )
done

echo "== Controls readcount metrics (no PoN) =="
echo "Repo root:          $REPO_ROOT"
echo "READCOUNT_QC_ROOT:  $READCOUNT_QC_ROOT"
echo "READCOUNT_VCF_DIR:  $READCOUNT_VCF_DIR"
echo "READCOUNT_BAM_DIR:  $READCOUNT_BAM_DIR"
echo "READCOUNT_REF_FASTA:$READCOUNT_REF_FASTA"
echo "READCOUNT_TO_TSV_PY:$READCOUNT_TO_TSV_PY"
echo "THREADS:            $THREADS"
echo

# ------------------------------------------------------------
# Stage 1: Build BED variant sites from annotated VCFs
# ------------------------------------------------------------
echo "=== Stage 1: BED from VCF ==="
for vcf in "${vcfs[@]}"; do
  [[ -f "$vcf" ]] || continue
  sample="$(sample_from_vcf "$vcf")"
  bed="$BEDS_DIR/${sample}.bed"
  [[ -s "$bed" ]] && { echo "[Stage1] SKIP $sample"; continue; }
  # BED columns here are: CHROM, POS0 (0-based start), POS (1-based end), ALT.
  # This is the format expected by bam-readcount -l for per-site extraction.
  run_to_file "$bed" bcftools query -f '%CHROM\t%POS0\t%POS\t%ALT\n' "$vcf"
done

# ------------------------------------------------------------
# Stage 2: Link recalibrated BAMs
# ------------------------------------------------------------
echo "=== Stage 2: Link BAMs ==="
for sample in "${samples[@]}"; do
  src_bam="$READCOUNT_BAM_DIR/${sample}.bam"
  dst_bam="$BAM_WORK_DIR/${sample}.bam"
  dst_bai="${dst_bam}.bai"

  [[ -s "$src_bam" ]] || die "Missing BAM for $sample: $src_bam"
  if [[ -s "${src_bam}.bai" ]]; then
    src_bai="${src_bam}.bai"
  elif [[ -s "$READCOUNT_BAM_DIR/${sample}.bam.bai" ]]; then
    # Accept either BAM index naming convention.
    src_bai="$READCOUNT_BAM_DIR/${sample}.bam.bai"
  else
    die "Missing BAM index for $sample: ${src_bam}.bai"
  fi

  [[ -e "$dst_bam" ]] || run ln -s "$src_bam" "$dst_bam"
  [[ -e "$dst_bai" ]] || run ln -s "$src_bai" "$dst_bai"
done

# ------------------------------------------------------------
# Stage 3: Remove duplicate reads (samtools -F 0x400)
# ------------------------------------------------------------
echo "=== Stage 3: Remove duplicates ==="
for sample in "${samples[@]}"; do
  in_bam="$BAM_WORK_DIR/${sample}.bam"
  out_bam="$BAM_NODUP_DIR/${sample}.nodup.bam"
  # 0x400 is the SAM duplicate flag; dropping these matches the methods text.
  [[ -s "$out_bam" ]] && { echo "[Stage3] SKIP $sample"; continue; }
  run_to_file "$out_bam" samtools view -@ "$THREADS" -b -F 0x400 "$in_bam"
  run samtools index -@ "$THREADS" "$out_bam"
done

# ------------------------------------------------------------
# Stage 4: Collect metrics with bam-readcount
# ------------------------------------------------------------
echo "=== Stage 4: bam-readcount ==="
for sample in "${samples[@]}"; do
  bam="$BAM_NODUP_DIR/${sample}.nodup.bam"
  bed="$BEDS_DIR/${sample}.bed"
  out="$READCOUNTS_DIR/${sample}.txt"

  [[ -f "$out" ]] && { echo "[Stage4] SKIP $sample"; continue; }
  [[ -s "$bam" ]] || die "Missing deduplicated BAM for $sample: $bam"
  [[ -f "$bed" ]] || die "Missing BED for $sample: $bed"

  if [[ ! -s "$bed" ]]; then
    echo "[Stage4] EMPTY BED $sample -> writing empty readcount file"
    if [[ "$DRY_RUN" == "0" ]]; then
      : > "$out"
    fi
    continue
  fi

  run_to_file "$out" bam-readcount -f "$READCOUNT_REF_FASTA" -l "$bed" "$bam"
done

# ------------------------------------------------------------
# Stage 5: Parse readcount files to TSV metrics
# ------------------------------------------------------------
echo "=== Stage 5: Parse readcount files ==="
run python3 "$READCOUNT_TO_TSV_PY" --input-dir "$READCOUNTS_DIR" --output-dir "$METRICS_DIR"

if [[ "$DRY_RUN" == "0" ]]; then
  missing_samples=()
  for sample in "${samples[@]}"; do
    # Parser writes one metrics TSV per sample.
    [[ -s "$METRICS_DIR/${sample}_metrics.tsv" ]] || missing_samples+=( "$sample" )
  done
  [[ "${#missing_samples[@]}" -eq 0 ]] || die "Final completeness check failed; missing metrics TSVs for: ${missing_samples[*]}"
fi

echo "=== Readcount metrics pipeline finished ==="
