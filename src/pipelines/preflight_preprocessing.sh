#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Preflight checks for preprocessing pipeline.
# Keeps checks out of preprocessing.sh (by design).
# ------------------------------------------------------------

# -----------------------
# Batch selection (toggle by commenting/uncommenting)
# -----------------------
BATCHES=(
   CJD_16_samples
   CJD_8_samples
   first_CJD_seq
   sequencing_of_dilutions
)

# -----------------------
# Repo root discovery
# -----------------------
find_repo_root() {
  local start="$1"
  local d="$start"
  while [[ "$d" != "/" ]]; do
    [[ -f "$d/Makefile" ]] && echo "$d" && return 0 # look for the Makefile
    d="$(dirname "$d")"
  done
  echo "ERROR: could not find repo root (Makefile not found)" >&2
  return 1
}

# folder this script is located in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# moves to the dir containing the Makefile
REPO_ROOT="$(find_repo_root "$(cd "$SCRIPT_DIR/../.." && pwd)")"

# -----------------------
# Config loading
# -----------------------

# If you already set an environment variable ENV_FILE, use it.
ENV_FILE_DEFAULT="$REPO_ROOT/config/preprocessing.env"

# Otherwise default to config/preprocessing.env
ENV_FILE="${ENV_FILE:-$ENV_FILE_DEFAULT}"

if [[ -f "$ENV_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$ENV_FILE"
else
  echo "ERROR: missing config file: $ENV_FILE" >&2
  echo "Hint: copy config/preprocessing.env.example -> config/preprocessing.env" >&2
  exit 1
fi

# Defaults: if they are not explicitly set in the env, use the following
ADAPTERS_FA="${ADAPTERS_FA:-resources/TruSeq3-PE.fa}"
DBSNP_VCF="${DBSNP_VCF:-resources/dbsnp_146.hg38.vcf.gz}"
FASTQ_DIR="${FASTQ_DIR:-fastq}"
FINAL_BAM_DIR="${FINAL_BAM_DIR:-results/final_bam}"
JAVA_MEM_GB="${JAVA_MEM_GB:-8}"
MILLS_VCF="${MILLS_VCF:-resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz}"
REF_DICT="${REF_DICT:-resources/chr2_chr4_chr20.dict}"
REF_FASTA="${REF_FASTA:-resources/chr2_chr4_chr20.fasta}"
RUNS_DIR="${RUNS_DIR:-runs/preprocessing}"
SAMPLES_TSV="${SAMPLES_TSV:-config/preprocessing_samples.tsv}"
THREADS="${THREADS:-8}"

# Convert relative paths to absolute
abs_path() {
  local p="$1"
  if [[ "$p" = /* ]]; then echo "$p"; else echo "$REPO_ROOT/$p"; fi
}

# Check whether a command exists
have() { command -v "$1" >/dev/null 2>&1; }

# Error out consistently
die() { echo "ERROR: $*" >&2; exit 1; }

in_list() {
  local x="$1"; shift
  for item in "$@"; do [[ "$x" == "$item" ]] && return 0; done
  return 1
}

# Require a file to exist and be non-empty
require_file() {
  local f="$1"
  [[ -s "$f" ]] || die "Missing/empty file: $f"
}

# Require a tool to exist
require_tool() {
  local t="$1"
  have "$t" || die "Tool not found in PATH: $t"
}

# print human-readable configuration summary
print_summary() {
  echo "== preprocessing preflight =="
  echo "Repo root:        $REPO_ROOT"
  echo "Env file:         $ENV_FILE"
  echo "Samples TSV:      $SAMPLES_TSV"
  echo "FASTQ_DIR:        $FASTQ_DIR"
  echo "RUNS_DIR:         $RUNS_DIR"
  echo "FINAL_BAM_DIR:    $FINAL_BAM_DIR"
  echo "REF_FASTA:        $REF_FASTA"
  echo "DBSNP_VCF:        $DBSNP_VCF"
  echo "MILLS_VCF:        $MILLS_VCF"
  echo "ADAPTERS_FA:      $ADAPTERS_FA"
  echo "THREADS:          $THREADS"
  echo "JAVA_MEM_GB:      $JAVA_MEM_GB"
  echo "Selected batches: ${BATCHES[*]:-(none)}"
}

# Tool checks (catches missing conda activation or wrong PATH)
check_tools() {
  # Keep aligned with what preprocessing.sh will use.
  require_tool bwa
  require_tool samtools
  require_tool gatk
  # If you call trimmomatic explicitly, uncomment:
  # require_tool trimmomatic
}

# Checks whether resources are present
check_resources() {
  local ref dbsnp mills adapters
  ref="$(abs_path "$REF_FASTA")"
  dict="$(abs_path "$REF_DICT")"
  dbsnp="$(abs_path "$DBSNP_VCF")"
  mills="$(abs_path "$MILLS_VCF")"
  adapters="$(abs_path "$ADAPTERS_FA")"

  require_file "$ref"
  require_file "${ref}.fai"
  require_file "$dict"

  require_file "$dbsnp"
  require_file "${dbsnp}.tbi"
  require_file "$mills"
  require_file "${mills}.tbi"
  require_file "$adapters"
}

# validating the samples TSV and FASTQs
check_samples_tsv() {
  local tsv
  tsv="$(abs_path "$SAMPLES_TSV")"
  require_file "$tsv"

  local header
  header="$(grep -v '^#' "$tsv" | head -n 1)"
  [[ "$header" == $'batch_id\tsample_id\tr1\tr2' ]] || die "Unexpected TSV header: $header"

  [[ "${#BATCHES[@]}" -gt 0 ]] || die "No batches selected (BATCHES array is empty)."

  local tmp
  tmp="$(mktemp)"
  trap 'rm -f "'"$tmp"'"' RETURN

  # Normalise to 4 columns
  tail -n +2 "$tsv" | grep -v '^#' | awk -F'\t' 'NF>0{print $1"\t"$2"\t"$3"\t"$4}' > "$tmp"

  # Confirm each selected batch exists
  for b in "${BATCHES[@]}"; do
    awk -F'\t' -v want="$b" '$1==want {found=1} END{exit !found}' "$tmp" \
    || die "Selected batch not found in TSV: $b"
  done

  # FASTQ existence checks
  local missing=0
  while IFS=$'\t' read -r batch sample r1 r2; do
    [[ -z "${batch// }" ]] && continue
    in_list "$batch" "${BATCHES[@]}" || continue

    [[ -n "$sample" ]] || die "Empty sample_id in TSV row: batch=$batch"
    [[ -n "$r1" && -n "$r2" ]] || die "Missing r1/r2 in TSV for $batch/$sample"

    local r1p r2p
    r1p="$(abs_path "$r1")"
    r2p="$(abs_path "$r2")"
    [[ -f "$r1p" ]] || { echo "MISSING: $r1p" >&2; missing=1; }
    [[ -f "$r2p" ]] || { echo "MISSING: $r2p" >&2; missing=1; }
  done < "$tmp"
  [[ "$missing" -eq 0 ]] || die "One or more FASTQ files are missing (see MISSING lines above)."

  # Enforce unique sample_id within selected batches (prevents overwriting results/final_bam/<sample>.bam)
  local dups
  dups="$(while IFS=$'\t' read -r batch sample r1 r2; do
            in_list "$batch" "${BATCHES[@]}" || continue
            echo "$sample"
          done < "$tmp" | sort | uniq -d)"
  [[ -z "$dups" ]] || die "Duplicate sample_id(s) within selected batches: $dups"

  echo "== Selected sample counts =="
  for b in "${BATCHES[@]}"; do
    local n
    n="$(awk -F'\t' -v want="$b" '$1==want {c++} END{print c+0}' "$tmp")"
    echo "  $b: $n"
  done
}

# runs the checklist in order
main() {
  print_summary
  check_tools
  check_resources
  check_samples_tsv
  echo "OK: preflight passed."
}

main "$@"
