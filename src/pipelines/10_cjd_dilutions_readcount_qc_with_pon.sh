#!/usr/bin/env bash
# ============================================================
# CJD+dilutions read/base quality metrics (with PoN):
#   Stage 1: Build BED variant sites from annotated VCFs
#   Stage 2: Link recalibrated BAMs
#   Stage 3: Remove duplicate reads (samtools -F 0x400)
#   Stage 4: Collect per-site metrics with bam-readcount
# ============================================================
set -euo pipefail
shopt -s nullglob

# Keep caller-provided values so ENV_FILE does not silently override them.
# Precedence: explicit CLI env vars > values from ENV_FILE > script defaults.
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
  # Optional: this stage can run entirely from defaults if ENV_FILE is absent.
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
MUTECT2_WITH_PON_OUT_ROOT="${MUTECT2_WITH_PON_OUT_ROOT:-runs/mutect2_cjd_dilutions_with_pon}"
if [[ "$MUTECT2_WITH_PON_OUT_ROOT" == run/* ]]; then
  # Backward compatibility for historical singular "run/" prefixes.
  MUTECT2_WITH_PON_OUT_ROOT="runs/${MUTECT2_WITH_PON_OUT_ROOT#run/}"
fi

MUTECT2_WITH_PON_READCOUNT_ROOT="${MUTECT2_WITH_PON_READCOUNT_ROOT:-$MUTECT2_WITH_PON_OUT_ROOT}"
if [[ "$MUTECT2_WITH_PON_READCOUNT_ROOT" == run/* ]]; then
  MUTECT2_WITH_PON_READCOUNT_ROOT="runs/${MUTECT2_WITH_PON_READCOUNT_ROOT#run/}"
fi

WITH_PON_READCOUNT_GROUPS="${WITH_PON_READCOUNT_GROUPS:-${WITH_PON_GROUPS:-cjd dilutions}}"
WITH_PON_READCOUNT_BAM_DIR="${WITH_PON_READCOUNT_BAM_DIR:-${FINAL_BAM_DIR:-results/final_bam}}"
WITH_PON_READCOUNT_REF_FASTA="${WITH_PON_READCOUNT_REF_FASTA:-${REF_FASTA:-resources/chr2_chr4_chr20.fasta}}"

# VCF input is expected from Stage 9 output.
WITH_PON_READCOUNT_VCF_SUBDIR="${WITH_PON_READCOUNT_VCF_SUBDIR:-annot_with_gnomad}"
# Keep this default aligned with Stage 9 to avoid mismatched hand-off paths.


to_abs() {
  local p="$1"
  # Allow both repo-relative and absolute paths in ENV overrides.
  if [[ "$p" = /* ]]; then
    echo "$p"
  else
    echo "$REPO_ROOT/$p"
  fi
}

READCOUNT_ROOT="$(to_abs "$MUTECT2_WITH_PON_READCOUNT_ROOT")"
WITH_PON_READCOUNT_BAM_DIR="$(to_abs "$WITH_PON_READCOUNT_BAM_DIR")"
WITH_PON_READCOUNT_REF_FASTA="$(to_abs "$WITH_PON_READCOUNT_REF_FASTA")"

# -----------------------
# Tiny helpers
# -----------------------
have() { command -v "$1" >/dev/null 2>&1; }
die() { echo "ERROR: $*" >&2; exit 1; }

run() {
  echo "+ $*"
  if [[ "$DRY_RUN" == "1" ]]; then
    # DRY_RUN prints commands only, allowing safe plan inspection.
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
  sample="${base%.func.af.vcf.gz}"
  if [[ "$sample" == "$base" ]]; then
    sample="${base%.vcf.gz}"
  fi
  # Sample stem is reused across BED/BAM/readcount filenames.
  echo "$sample"
}

# -----------------------
# Tool + file checks
# -----------------------
have bcftools      || die "bcftools not in PATH (did you activate conda?)"
have samtools      || die "samtools not in PATH (did you activate conda?)"
have bam-readcount || die "bam-readcount not in PATH (did you activate conda?)"

[[ -d "$READCOUNT_ROOT" ]] || die "Missing root directory: $READCOUNT_ROOT"
[[ -d "$WITH_PON_READCOUNT_BAM_DIR" ]] || die "Missing BAM directory: $WITH_PON_READCOUNT_BAM_DIR"
[[ -s "$WITH_PON_READCOUNT_REF_FASTA" ]] || die "Missing FASTA: $WITH_PON_READCOUNT_REF_FASTA"
[[ -s "$WITH_PON_READCOUNT_REF_FASTA.fai" ]] || die "Missing FASTA index: $WITH_PON_READCOUNT_REF_FASTA.fai"

read -r -a groups <<< "$WITH_PON_READCOUNT_GROUPS"
[[ "${#groups[@]}" -gt 0 ]] || die "WITH_PON_READCOUNT_GROUPS is empty"
# Group list controls independent cjd/dilutions processing branches.

echo "== CJD+dilutions readcount metrics (with PoN) =="
echo "Repo root:              $REPO_ROOT"
echo "READCOUNT_ROOT:         $READCOUNT_ROOT"
echo "WITH_PON_GROUPS:        $WITH_PON_READCOUNT_GROUPS"
echo "WITH_PON_READCOUNT_BAM_DIR:  $WITH_PON_READCOUNT_BAM_DIR"
echo "WITH_PON_READCOUNT_REF_FASTA: $WITH_PON_READCOUNT_REF_FASTA"
echo "THREADS:                $THREADS"
echo

for group in "${groups[@]}"; do
  echo "=== Group: $group ==="

  VCF_DIR="$READCOUNT_ROOT/$group/$WITH_PON_READCOUNT_VCF_SUBDIR"
  READCOUNT_QC_ROOT="$READCOUNT_ROOT/$group/readcount_qc"

  BEDS_DIR="$READCOUNT_QC_ROOT/beds"
  BAM_WORK_DIR="$READCOUNT_QC_ROOT/bam_work"
  BAM_NODUP_DIR="$READCOUNT_QC_ROOT/bam_nodup"
  READCOUNTS_DIR="$READCOUNT_QC_ROOT/readcounts"

  mkdir -p "$BEDS_DIR" "$BAM_WORK_DIR" "$BAM_NODUP_DIR" "$READCOUNTS_DIR"

  [[ -d "$VCF_DIR" ]] || die "Missing VCF directory for group '$group': $VCF_DIR"

  vcfs=( "$VCF_DIR"/*.func.af.vcf.gz )
  if [[ "${#vcfs[@]}" -eq 0 ]]; then
    vcfs=( "$VCF_DIR"/*.vcf.gz )
  fi
  [[ "${#vcfs[@]}" -gt 0 ]] || die "No annotated VCFs found for group '$group' in: $VCF_DIR"
  # Accept both current and legacy VCF suffix conventions.

  samples=()
  for vcf in "${vcfs[@]}"; do
    [[ -f "$vcf" ]] || continue
    samples+=( "$(sample_from_vcf "$vcf")" )
  done

  # ------------------------------------------------------------
  # Stage 1: Build BED variant sites from annotated VCFs
  # ------------------------------------------------------------
  echo "=== Stage 1 ($group): BED from VCF ==="
  for vcf in "${vcfs[@]}"; do
    [[ -f "$vcf" ]] || continue
    sample="$(sample_from_vcf "$vcf")"
    bed="$BEDS_DIR/${sample}.bed"
    [[ -s "$bed" ]] && { echo "[Stage1][$group] SKIP $sample"; continue; }
    # BED columns: CHROM, POS0 (0-based start), POS (1-based end), ALT.
    run_to_file "$bed" bcftools query -f '%CHROM\t%POS0\t%POS\t%ALT\n' "$vcf"
  done

  # ------------------------------------------------------------
  # Stage 2: Link recalibrated BAMs
  # ------------------------------------------------------------
  echo "=== Stage 2 ($group): Link BAMs ==="
  for sample in "${samples[@]}"; do
    src_bam="$WITH_PON_READCOUNT_BAM_DIR/${sample}.bam"
    dst_bam="$BAM_WORK_DIR/${sample}.bam"
    dst_bai="${dst_bam}.bai"

    [[ -s "$src_bam" ]] || die "Missing BAM for $sample: $src_bam"
    if [[ -s "${src_bam}.bai" ]]; then
      src_bai="${src_bam}.bai"
    elif [[ -s "$WITH_PON_READCOUNT_BAM_DIR/${sample}.bam.bai" ]]; then
      src_bai="$WITH_PON_READCOUNT_BAM_DIR/${sample}.bam.bai"
    else
      die "Missing BAM index for $sample: ${src_bam}.bai"
    fi

    # Use symlinks to avoid duplicating large BAM/BAM index files in run folders.
    [[ -e "$dst_bam" ]] || run ln -s "$src_bam" "$dst_bam"
    [[ -e "$dst_bai" ]] || run ln -s "$src_bai" "$dst_bai"
  done

  # ------------------------------------------------------------
  # Stage 3: Remove duplicate reads (samtools -F 0x400)
  # ------------------------------------------------------------
  echo "=== Stage 3 ($group): Remove duplicates ==="
  for sample in "${samples[@]}"; do
    in_bam="$BAM_WORK_DIR/${sample}.bam"
    out_bam="$BAM_NODUP_DIR/${sample}.nodup.bam"
    [[ -s "$out_bam" ]] && { echo "[Stage3][$group] SKIP $sample"; continue; }
    run_to_file "$out_bam" samtools view -@ "$THREADS" -b -F 0x400 "$in_bam"
    run samtools index -@ "$THREADS" "$out_bam"
  done

  # ------------------------------------------------------------
  # Stage 4: Collect metrics with bam-readcount
  # ------------------------------------------------------------
  echo "=== Stage 4 ($group): bam-readcount ==="
  for sample in "${samples[@]}"; do
    bam="$BAM_NODUP_DIR/${sample}.nodup.bam"
    bed="$BEDS_DIR/${sample}.bed"
    out="$READCOUNTS_DIR/${sample}.txt"

    # Presence of output file marks completion (empty files are valid for empty BEDs).
    [[ -f "$out" ]] && { echo "[Stage4][$group] SKIP $sample"; continue; }
    if [[ "$DRY_RUN" == "0" ]]; then
      [[ -s "$bam" ]] || die "Missing deduplicated BAM for $sample: $bam"
      [[ -f "$bed" ]] || die "Missing BED for $sample: $bed"

      if [[ ! -s "$bed" ]]; then
        echo "[Stage4][$group] EMPTY BED $sample -> writing empty readcount file"
        # Keep placeholder output so downstream completeness checks remain deterministic.
        : > "$out"
        continue
      fi
    fi

    run_to_file "$out" bam-readcount -f "$WITH_PON_READCOUNT_REF_FASTA" -l "$bed" "$bam"
  done

  if [[ "$DRY_RUN" == "0" ]]; then
    # Final per-group contract: one readcount file per discovered sample.
    missing_samples=()
    for sample in "${samples[@]}"; do
      [[ -f "$READCOUNTS_DIR/${sample}.txt" ]] || missing_samples+=( "$sample" )
    done
    [[ "${#missing_samples[@]}" -eq 0 ]] || die "Final completeness check failed for group '$group'; missing readcount files for: ${missing_samples[*]}"
  fi

  echo "=== Group finished: $group ==="
  echo
done

echo "=== CJD+dilutions readcount collection (with PoN) finished ==="
