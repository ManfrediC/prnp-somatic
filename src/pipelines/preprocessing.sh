#!/usr/bin/env bash
set -euo pipefail

# -----------------------
# User toggles (batches)
# -----------------------
BATCHES=(
  # CJD_16_samples
   CJD_8_samples
  # first_CJD_seq
  # sequencing_of_dilutions
)

# -----------------------
# Defaults (override via env)
# -----------------------
THREADS="${THREADS:-8}"
JAVA_MEM_GB="${JAVA_MEM_GB:-8}"
FORCE="${FORCE:-0}"          # 1 to overwrite existing outputs
DRY_RUN="${DRY_RUN:-0}"      # 1 to print commands only

# -----------------------
# Repo root resolution
# -----------------------
find_repo_root() {
  local start="$1"
  local d="$start"
  while [[ "$d" != "/" ]]; do
    [[ -f "$d/Makefile" ]] && echo "$d" && return 0 # walk up until Makefile is found
    d="$(dirname "$d")"
  done
  echo "ERROR: could not find repo root (Makefile not found)" >&2 # sends output to stderr
  return 1
}

REPO_ROOT="$(find_repo_root "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)")"

# -----------------------
# Paths (override via env)
# -----------------------
SAMPLES_TSV="${SAMPLES_TSV:-$REPO_ROOT/config/preprocessing_samples.tsv}"
FASTQ_DIR="${FASTQ_DIR:-$REPO_ROOT/fastq}"
RUNS_DIR="${RUNS_DIR:-$REPO_ROOT/runs/preprocessing}"
FINAL_BAM_DIR="${FINAL_BAM_DIR:-$REPO_ROOT/results/final_bam}"

REF_FASTA="${REF_FASTA:-$REPO_ROOT/resources/chr2_chr4_chr20.fasta}"
DBSNP_VCF="${DBSNP_VCF:-$REPO_ROOT/resources/dbsnp_146.hg38.vcf.gz}"
MILLS_VCF="${MILLS_VCF:-$REPO_ROOT/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz}"
ADAPTERS_FA="${ADAPTERS_FA:-$REPO_ROOT/resources/adapters/TruSeq3-PE.fa}"

# -----------------------
# Helpers
# -----------------------

# timestamped logging
log() { echo "[$(date -Is)] $*"; }

# hard failure
die() { echo "ERROR: $*" >&2; exit 1; }

# check if a command exists
have() { command -v "$1" >/dev/null 2>&1; }

run_cmd() {
  if [[ "$DRY_RUN" == "1" ]]; then
    echo "+ $*"
  else
    eval "$@"
  fi
}

in_list() {
  local x="$1"; shift
  for item in "$@"; do [[ "$x" == "$item" ]] && return 0; done
  return 1
}

ensure_file() { [[ -s "$1" ]] || die "Missing/empty file: $1"; }

maybe_skip() {
  local out="$1"
  [[ "$FORCE" == "1" ]] && return 1
  [[ -s "$out" ]]
}

# -----------------------
# Preflight: validate that key inputs exist
# -----------------------
preflight() {
  ensure_file "$SAMPLES_TSV"
  ensure_file "$REF_FASTA"
  ensure_file "${REF_FASTA}.fai"

  ensure_file "$DBSNP_VCF"
  ensure_file "${DBSNP_VCF}.tbi"
  ensure_file "$MILLS_VCF"
  ensure_file "${MILLS_VCF}.tbi"

  have bwa || die "bwa not found in PATH"
  have samtools || die "samtools not found in PATH"
  have gatk || die "gatk not found in PATH"
}

# -----------------------
# Per-sample stages (stubs)
# -----------------------
process_sample() {
  local batch="$1" sample="$2" r1="$3" r2="$4"

  local sample_dir="$RUNS_DIR/$batch/$sample"
  local log_dir="$sample_dir/logs"
  mkdir -p "$log_dir" "$FINAL_BAM_DIR"

  # Resolve FASTQs (paths in TSV are repo-relative)
  local r1_path="$REPO_ROOT/$r1"
  local r2_path="$REPO_ROOT/$r2"
  ensure_file "$r1_path"
  ensure_file "$r2_path"

  # Output naming conventions
  local trimmed_r1="$sample_dir/${sample}.trim.R1.fastq.gz"
  local trimmed_r2="$sample_dir/${sample}.trim.R2.fastq.gz"
  local final_bam="$sample_dir/${sample}.bwa.picard.markedDup.recal.bam"
  local final_bai="${final_bam}.bai"

  # Stage: trim (stub)
  if ! maybe_skip "$trimmed_r1"; then
    log "[$batch/$sample] trim"
    # run_cmd "trimmomatic PE -threads $THREADS ... > '$log_dir/trim.log' 2>&1"
  fi

  # Stage: align+sort (stub)
  # Stage: add RG (stub)
  # Stage: mark dup (stub)
  # Stage: BQSR (stub)

  # Publish
  if [[ -s "$final_bam" && -s "$final_bai" ]]; then # Only publish if the run produced a BAM + BAI in the per-sample run directory
    if ! maybe_skip "$FINAL_BAM_DIR/${sample}.bam"; then
      log "[$batch/$sample] publish final BAM"
      run_cmd "cp -f '$final_bam' '$FINAL_BAM_DIR/${sample}.bam'"
      run_cmd "cp -f '$final_bai' '$FINAL_BAM_DIR/${sample}.bam.bai'"
    fi
  fi
}

# -----------------------
# TSV loop
# -----------------------
main() {
  preflight

  log "Repo root: $REPO_ROOT"
  log "Selected batches: ${BATCHES[*]}"

  # Parse header indices
  local header
  header="$(grep -v '^#' "$SAMPLES_TSV" | head -n 1)"
  [[ -n "$header" ]] || die "No header found in $SAMPLES_TSV"

  # Expected columns (we keep it strict)
  [[ "$header" == $'batch_id\tsample_id\tr1\tr2' ]] || die "Unexpected header: $header"

  # Iterate rows
  tail -n +2 "$SAMPLES_TSV" | grep -v '^#' | while IFS=$'\t' read -r batch sample r1 r2; do
    [[ -z "${batch// }" ]] && continue
    in_list "$batch" "${BATCHES[@]}" || continue
    process_sample "$batch" "$sample" "$r1" "$r2"
  done
}

main "$@" # passes all command-line arguments to main
