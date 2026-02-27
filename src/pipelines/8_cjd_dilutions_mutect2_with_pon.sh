#!/usr/bin/env bash
# ============================================================
# Mutect2 tumour-only with PoN for CJD samples and dilutions.
# Inputs: results/final_bam/{CJD*.bam, dilution BAMs}
# Outputs: runs/mutect2_cjd_dilutions_with_pon/{cjd,dilutions}/{vcf,f1r2}/
# ============================================================
set -euo pipefail
shopt -s nullglob

# Keep caller-provided values so ENV_FILE does not silently override them.
# Precedence: explicit CLI env vars > values from ENV_FILE > script defaults.
CLI_DRY_RUN="${DRY_RUN-}"
CLI_THREADS="${THREADS-}"
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
  # Optional: script remains runnable with defaults when no env file is present.
  # shellcheck disable=SC1090
  source "$ENV_FILE"
fi

# -----------------------
# Defaults (override via env or ENV_FILE)
# -----------------------
DRY_RUN="${CLI_DRY_RUN:-${DRY_RUN:-0}}"
THREADS="${CLI_THREADS:-${THREADS:-8}}"
JAVA_MEM_GB="${CLI_JAVA_MEM_GB:-${JAVA_MEM_GB:-8}}"

# -----------------------
# Path defaults (repo-relative)
# -----------------------
MUTECT2_WITH_PON_OUT_ROOT="${MUTECT2_WITH_PON_OUT_ROOT:-runs/mutect2_cjd_dilutions_with_pon}"
if [[ "$MUTECT2_WITH_PON_OUT_ROOT" == run/* ]]; then
  # Backward compatibility for historical singular "run/" paths.
  MUTECT2_WITH_PON_OUT_ROOT="runs/${MUTECT2_WITH_PON_OUT_ROOT#run/}"
fi

FINAL_BAM_DIR="${FINAL_BAM_DIR:-results/final_bam}"
REF_FASTA="${REF_FASTA:-resources/chr2_chr4_chr20.fasta}"
INTERVALS="${INTERVALS:-resources/capture_targets.interval_list}"
GNOMAD_VCF="${GNOMAD_VCF:-resources/somatic-hg38_af-only-gnomad.hg38.vcf.gz}"
PON_VCF="${PON_VCF:-results/mutect2_controls_pon/panel_of_normals/CJD_controls_PoN.vcf.gz}"

# Dilution order kept explicit for reproducibility and stable output ordering.
DILUTION_SAMPLES="${DILUTION_SAMPLES:-NA100_undil NA100_1to10 NA99A1_undil A100_1to2 NA99A1_1to5 NA995A05_undil NA100_1to2}"

to_abs() {
  local p="$1"
  # Accept both absolute and repo-relative overrides from ENV/CLI.
  if [[ "$p" = /* ]]; then
    echo "$p"
  else
    echo "$REPO_ROOT/$p"
  fi
}

OUT_ROOT="$(to_abs "$MUTECT2_WITH_PON_OUT_ROOT")"
FINAL_BAM_DIR="$(to_abs "$FINAL_BAM_DIR")"
REF_FASTA="$(to_abs "$REF_FASTA")"
INTERVALS="$(to_abs "$INTERVALS")"
GNOMAD_VCF="$(to_abs "$GNOMAD_VCF")"
PON_VCF="$(to_abs "$PON_VCF")"

CJD_VCF_DIR="$OUT_ROOT/cjd/vcf"
CJD_F1R2_DIR="$OUT_ROOT/cjd/f1r2"
DIL_VCF_DIR="$OUT_ROOT/dilutions/vcf"
DIL_F1R2_DIR="$OUT_ROOT/dilutions/f1r2"

mkdir -p "$CJD_VCF_DIR" "$CJD_F1R2_DIR" "$DIL_VCF_DIR" "$DIL_F1R2_DIR"

# -----------------------
# Tiny helpers
# -----------------------
have() { command -v "$1" >/dev/null 2>&1; }
die() { echo "ERROR: $*" >&2; exit 1; }

run() {
  echo "+ $*"
  if [[ "$DRY_RUN" == "1" ]]; then
    # DRY_RUN prints commands only, so reviewers can inspect execution plans.
    return 0
  fi
  "$@"
}

has_vcf_index() {
  local vcf="$1"
  # Accept common index sidecars produced by GATK/htslib tooling.
  [[ -s "${vcf}.tbi" || -s "${vcf}.csi" || -s "${vcf}.idx" ]]
}

# -----------------------
# Tool + file checks
# -----------------------
have gatk     || die "gatk not in PATH (did you activate conda?)"
have samtools || die "samtools not in PATH (did you activate conda?)"

[[ -d "$FINAL_BAM_DIR" ]] || die "Missing FINAL_BAM_DIR: $FINAL_BAM_DIR"
[[ -s "$REF_FASTA" ]] || die "Missing REF_FASTA: $REF_FASTA"
[[ -s "$INTERVALS" ]] || die "Missing INTERVALS: $INTERVALS"
[[ -s "$GNOMAD_VCF" ]] || die "Missing GNOMAD_VCF: $GNOMAD_VCF"
[[ -s "$PON_VCF" ]] || die "Missing PON_VCF: $PON_VCF"

[[ -s "${GNOMAD_VCF}.tbi" || -s "${GNOMAD_VCF}.csi" ]] || die "Missing index for GNOMAD_VCF: ${GNOMAD_VCF}.tbi/.csi"
if ! has_vcf_index "$PON_VCF"; then
  # Build PoN index on demand to avoid manual pre-indexing requirements.
  run gatk IndexFeatureFile -I "$PON_VCF"
fi
has_vcf_index "$PON_VCF" || die "Missing index for PON_VCF after indexing attempt: $PON_VCF(.tbi/.csi/.idx)"

cjd_bams=( "$FINAL_BAM_DIR"/CJD*.bam )
[[ "${#cjd_bams[@]}" -gt 0 ]] || die "No CJD BAMs found under: $FINAL_BAM_DIR (expected CJD*.bam)"

read -r -a dilutions <<< "$DILUTION_SAMPLES"
[[ "${#dilutions[@]}" -gt 0 ]] || die "DILUTION_SAMPLES is empty"
# Explicit sample order keeps logs and output processing deterministic.

dil_bams=()
for sample in "${dilutions[@]}"; do
  bam="$FINAL_BAM_DIR/${sample}.bam"
  [[ -s "$bam" ]] || die "Missing dilution BAM for sample '$sample': $bam"
  dil_bams+=( "$bam" )
done

echo "== Mutect2 with PoN (CJD + dilutions) =="
echo "Repo root:        $REPO_ROOT"
echo "FINAL_BAM_DIR:    $FINAL_BAM_DIR"
echo "OUT_ROOT:         $OUT_ROOT"
echo "REF_FASTA:        $REF_FASTA"
echo "INTERVALS:        $INTERVALS"
echo "GNOMAD_VCF:       $GNOMAD_VCF"
echo "PON_VCF:          $PON_VCF"
echo "CJD count:        ${#cjd_bams[@]}"
echo "Dilution count:   ${#dil_bams[@]}"
echo "THREADS:          $THREADS"
echo "JAVA_MEM_GB:      $JAVA_MEM_GB"
echo

run_mutect() {
  local bam="$1"
  local vcf_dir="$2"
  local f1r2_dir="$3"
  local sample out_vcf out_f1r2

  sample="$(basename "$bam" .bam)"
  # Sample name is derived from BAM stem and reused across all output artefacts.

  if [[ ! -s "${bam}.bai" ]]; then
    # Build BAM index lazily so reruns do not spend time re-indexing existing files.
    run samtools index -@ "$THREADS" "$bam"
  fi

  out_vcf="$vcf_dir/${sample}.raw.vcf.gz"
  out_f1r2="$f1r2_dir/${sample}.f1r2.tar.gz"

  # Skip if VCF already exists; this keeps the stage resumable.
  if [[ -s "$out_vcf" ]]; then
    echo "SKIP: $out_vcf exists"
    return 0
  fi

  run gatk --java-options "-Xmx${JAVA_MEM_GB}g" Mutect2 \
    -R "$REF_FASTA" \
    -I "$bam" \
    -tumor-sample "$sample" \
    --panel-of-normals "$PON_VCF" \
    --germline-resource "$GNOMAD_VCF" \
    --af-of-alleles-not-in-resource 0.0000025 \
    --f1r2-tar-gz "$out_f1r2" \
    --intervals "$INTERVALS" \
    -O "$out_vcf"
  # Mutect2 writes compressed VCF directly; downstream stages consume *.raw.vcf.gz.
}

echo "=== CJD samples ==="
for bam in "${cjd_bams[@]}"; do
  run_mutect "$bam" "$CJD_VCF_DIR" "$CJD_F1R2_DIR"
done

echo
echo "=== Dilution samples ==="
for bam in "${dil_bams[@]}"; do
  run_mutect "$bam" "$DIL_VCF_DIR" "$DIL_F1R2_DIR"
done

if [[ "$DRY_RUN" == "0" ]]; then
  # Final cross-group completeness check: every declared BAM must yield one raw VCF.
  missing=()
  for bam in "${cjd_bams[@]}" "${dil_bams[@]}"; do
    sample="$(basename "$bam" .bam)"
    if [[ -s "$CJD_VCF_DIR/${sample}.raw.vcf.gz" || -s "$DIL_VCF_DIR/${sample}.raw.vcf.gz" ]]; then
      continue
    fi
    missing+=( "$sample" )
  done
  [[ "${#missing[@]}" -eq 0 ]] || die "Final completeness check failed; missing raw VCFs for: ${missing[*]}"
fi

echo "=== Mutect2 with PoN finished ==="
