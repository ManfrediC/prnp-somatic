#!/usr/bin/env bash
# ============================================================
# Controls-only Panel of Normals (PoN) creation:
#   Stage 1: Verify/index control filtered VCF inputs
#   Stage 2: Merge control VCFs into one multi-sample VCF
#   Stage 3: Build PoN with GATK CreateSomaticPanelOfNormals
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

PON_INPUT_DIR="${PON_INPUT_DIR:-$MUTECT2_CONTROLS_OUT_ROOT/filtered}"
PON_OUTPUT_ROOT="${PON_OUTPUT_ROOT:-results/mutect2_controls_pon/panel_of_normals}"
PON_MERGED_VCF="${PON_MERGED_VCF:-$PON_OUTPUT_ROOT/controls_multisample.filtered.vcf.gz}"
PON_VCF="${PON_VCF:-$PON_OUTPUT_ROOT/CJD_controls_PoN.vcf.gz}"
PON_CONTROLS="${PON_CONTROLS:-Ctrl1 Ctrl2 Ctrl3 Ctrl4 Ctrl5 Ctrl7}"

REF_FASTA="${REF_FASTA:-resources/chr2_chr4_chr20.fasta}"

to_abs() {
  local p="$1"
  if [[ "$p" = /* ]]; then
    echo "$p"
  else
    echo "$REPO_ROOT/$p"
  fi
}

PON_INPUT_DIR="$(to_abs "$PON_INPUT_DIR")"
PON_OUTPUT_ROOT="$(to_abs "$PON_OUTPUT_ROOT")"
PON_MERGED_VCF="$(to_abs "$PON_MERGED_VCF")"
PON_VCF="$(to_abs "$PON_VCF")"
REF_FASTA="$(to_abs "$REF_FASTA")"

mkdir -p "$PON_OUTPUT_ROOT"

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

has_vcf_index() {
  local vcf="$1"
  [[ -s "${vcf}.tbi" || -s "${vcf}.csi" ]]
}

# -----------------------
# Tool + file checks
# -----------------------
have bcftools || die "bcftools not in PATH (did you activate conda?)"
have gatk     || die "gatk not in PATH (did you activate conda?)"

[[ -d "$PON_INPUT_DIR" ]] || die "Missing PON_INPUT_DIR: $PON_INPUT_DIR"
[[ -s "$REF_FASTA" ]] || die "Missing REF_FASTA: $REF_FASTA"
[[ -s "$REF_FASTA.fai" ]] || die "Missing FASTA index: $REF_FASTA.fai"

read -r -a controls <<< "$PON_CONTROLS"
[[ "${#controls[@]}" -gt 0 ]] || die "PON_CONTROLS is empty"

input_vcfs=()
for sample in "${controls[@]}"; do
  vcf="$PON_INPUT_DIR/${sample}.filtered.vcf.gz"
  [[ -s "$vcf" ]] || die "Missing filtered VCF for control sample $sample: $vcf"
  input_vcfs+=( "$vcf" )
done

echo "== Controls PoN creation =="
echo "Repo root:        $REPO_ROOT"
echo "PON_INPUT_DIR:    $PON_INPUT_DIR"
echo "PON_OUTPUT_ROOT:  $PON_OUTPUT_ROOT"
echo "PON_MERGED_VCF:   $PON_MERGED_VCF"
echo "PON_VCF:          $PON_VCF"
echo "REF_FASTA:        $REF_FASTA"
echo "PON_CONTROLS:     ${controls[*]}"
echo "JAVA_MEM_GB:      $JAVA_MEM_GB"
echo

# ------------------------------------------------------------
# Stage 1: Ensure indexes for all control filtered VCFs
# ------------------------------------------------------------
echo "=== Stage 1: Validate/index control VCFs ==="
for vcf in "${input_vcfs[@]}"; do
  if has_vcf_index "$vcf"; then
    echo "[Stage1] OK $(basename "$vcf")"
    continue
  fi
  run bcftools index -f "$vcf"
done

# ------------------------------------------------------------
# Stage 2: Merge controls into one multi-sample VCF
# ------------------------------------------------------------
echo "=== Stage 2: Merge control VCFs ==="
if [[ -s "$PON_MERGED_VCF" ]] && has_vcf_index "$PON_MERGED_VCF"; then
  echo "[Stage2] SKIP merged VCF exists: $PON_MERGED_VCF"
else
  run bcftools merge -m all -O z -o "$PON_MERGED_VCF" "${input_vcfs[@]}"
  run bcftools index -f "$PON_MERGED_VCF"
fi

# ------------------------------------------------------------
# Stage 3: Create PoN VCF
# ------------------------------------------------------------
echo "=== Stage 3: CreateSomaticPanelOfNormals ==="
if [[ -s "$PON_VCF" ]] && has_vcf_index "$PON_VCF"; then
  echo "[Stage3] SKIP PoN VCF exists: $PON_VCF"
else
  run gatk --java-options "-Xmx${JAVA_MEM_GB}g" CreateSomaticPanelOfNormals \
    -V "$PON_MERGED_VCF" \
    -R "$REF_FASTA" \
    -O "$PON_VCF"
  run bcftools index -f "$PON_VCF"
fi

if [[ "$DRY_RUN" == "0" ]]; then
  [[ -s "$PON_MERGED_VCF" ]] || die "Missing merged output: $PON_MERGED_VCF"
  has_vcf_index "$PON_MERGED_VCF" || die "Missing merged VCF index: $PON_MERGED_VCF(.tbi/.csi)"
  [[ -s "$PON_VCF" ]] || die "Missing PoN output: $PON_VCF"
  has_vcf_index "$PON_VCF" || die "Missing PoN index: $PON_VCF(.tbi/.csi)"
  echo -n "PoN site count: "
  bcftools view -H "$PON_VCF" | wc -l
fi

echo "=== Controls PoN creation finished ==="
