#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ------------------------------------------------------------
# Mutect2 tumour-only (no PoN), CONTROLS ONLY.
# Inputs: results/final_bam/Ctrl*.bam (+ .bai)
# Outputs: run/mutect2_controls_no_pon/{vcf,f1r2}/
# ------------------------------------------------------------

# -----------------------
# Batch selection (optional; controls are detected by BAM name)
# -----------------------
# Not used yet; kept only if you later want to restrict controls by batch.
# BATCHES=(
#   CJD_16_samples
#   CJD_8_samples
#   first_CJD_seq
#   sequencing_of_dilutions
# )

# -----------------------
# Defaults (override via env or ENV_FILE)
# -----------------------
DRY_RUN="${DRY_RUN:-0}"
THREADS="${THREADS:-8}"
JAVA_MEM_GB="${JAVA_MEM_GB:-8}"

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
# Inputs/outputs (repo-relative defaults)
# -----------------------
FINAL_BAM_DIR="${FINAL_BAM_DIR:-results/final_bam}"
MUTECT2_CONTROLS_OUT_ROOT="${MUTECT2_CONTROLS_OUT_ROOT:-run/mutect2_controls_no_pon}"
OUT_ROOT="${OUT_ROOT:-$MUTECT2_CONTROLS_OUT_ROOT}"

REF_FASTA="${REF_FASTA:-resources/chr2_chr4_chr20.fasta}"
INTERVALS="${INTERVALS:-resources/capture_targets.interval_list}"

# Paper text says “resource bundle”; Mutect2 tumour-only typically uses gnomAD AF-only.
# If you haven't added it to resources yet, set GNOMAD_VCF in your env.
GNOMAD_VCF="${GNOMAD_VCF:-resources/somatic-hg38_af-only-gnomad.hg38.vcf.gz}"

mkdir -p "$OUT_ROOT/vcf" "$OUT_ROOT/f1r2"

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
# Tool + file checks (keep minimal)
# -----------------------
have gatk     || die "gatk not in PATH (did you activate conda?)"
have samtools || die "samtools not in PATH (did you activate conda?)"

[[ -s "$REF_FASTA" ]]   || die "Missing REF_FASTA: $REF_FASTA"
[[ -s "$INTERVALS" ]]   || die "Missing INTERVALS: $INTERVALS"
[[ -s "$GNOMAD_VCF" ]]  || die "Missing GNOMAD_VCF: $GNOMAD_VCF (set GNOMAD_VCF in ENV_FILE if needed)"
[[ -s "${GNOMAD_VCF}.tbi" || -s "${GNOMAD_VCF}.csi" ]] || die "Missing index for GNOMAD_VCF: ${GNOMAD_VCF}.tbi/.csi"

# -----------------------
# Main loop: controls only
# -----------------------
bams=( "$FINAL_BAM_DIR"/Ctrl*.bam )
[[ "${#bams[@]}" -gt 0 ]] || die "No control BAMs found under: $FINAL_BAM_DIR (expected Ctrl*.bam)"

echo "== Mutect2 tumour-only (controls only) =="
echo "Repo root:      $REPO_ROOT"
echo "BAM dir:        $FINAL_BAM_DIR"
echo "Output root:    $OUT_ROOT"
echo "REF_FASTA:      $REF_FASTA"
echo "INTERVALS:      $INTERVALS"
echo "GNOMAD_VCF:     $GNOMAD_VCF"
echo "THREADS:        $THREADS"
echo "JAVA_MEM_GB:    $JAVA_MEM_GB"
echo

for bam in "${bams[@]}"; do
  [[ -f "$bam" ]] || continue
  sample="$(basename "$bam" .bam)"

  echo "== $sample =="

  # Ensure BAM index exists
  if [[ ! -f "${bam}.bai" ]]; then
    run samtools index -@ "$THREADS" "$bam"
  fi

  out_vcf="$OUT_ROOT/vcf/${sample}.raw.vcf.gz"
  out_f1r2="$OUT_ROOT/f1r2/${sample}.f1r2.tar.gz"

  if [[ -s "$out_vcf" ]]; then
    echo "SKIP: $out_vcf exists"
    continue
  fi

  run gatk --java-options "-Xmx${JAVA_MEM_GB}g" Mutect2 \
    -R "$REF_FASTA" \
    -I "$bam" \
    -tumor-sample "$sample" \
    --germline-resource "$GNOMAD_VCF" \
    --af-of-alleles-not-in-resource 0.0000025 \
    --f1r2-tar-gz "$out_f1r2" \
    --intervals "$INTERVALS" \
    -O "$out_vcf"

  echo
done

echo "DONE"
