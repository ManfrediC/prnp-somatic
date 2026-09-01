#!/usr/bin/env bash
# Stop on command failures, unset variables and failed pipeline steps.
set -euo pipefail

# ------------------------------------------------------------
# Collect fixed-site read counts from the high and low mixtures.
# ------------------------------------------------------------

# Resolve project paths from this script, regardless of the working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd -P)"
cd "$REPO_ROOT"

# Use the frozen marker set; this stage makes no recovery decisions.
SAMPLES_TSV="src2/spikein/samples.tsv"
MARKERS="results2/spikein/markers/informative_markers.tsv"
MARKER_SETTINGS="results2/spikein/markers/run_settings.tsv"
SOURCE_SETTINGS="results2/spikein/discovery/run_settings.tsv"
REF_FASTA="resources/references/snv/chr2_chr4_chr20.fasta"
OUT_ROOT="${OUT_ROOT:-results2/spikein/readcount_qc/mixtures}"
DRY_RUN="${DRY_RUN:-1}"
MAX_DEPTH=10000000
export LC_ALL=C

# Report failed checks and print file-producing commands before execution.
die() { echo "ERROR: $*" >&2; exit 1; }
run_to_file() {
  local out="$1"
  shift
  printf '+ %s > %s\n' "$*" "$out"
  [[ "$DRY_RUN" == "1" ]] || "$@" > "$out"
}

# Keep each validation block small and named for the condition it checks.
check_inputs() {
  # Check the tool and inputs before creating any output.
  [[ "$#" -eq 0 ]] || die "Use DRY_RUN and OUT_ROOT, not positional arguments"
  [[ "$DRY_RUN" == "0" || "$DRY_RUN" == "1" ]] || die "DRY_RUN must be 0 or 1"
  command -v bam-readcount >/dev/null || die "bam-readcount not in PATH (activate prnp-spikein)"
  for input in Makefile "$SAMPLES_TSV" "$MARKERS" "$MARKER_SETTINGS" "$SOURCE_SETTINGS" \
    "$REF_FASTA" "$REF_FASTA.fai"; do
    [[ -s "$input" ]] || die "Missing input: $input"
  done
}

check_frozen_inputs() {
  # Require the exact marker table frozen by stage 5.
  mapfile -t marker_hashes < <(awk -F '\t' '$1 == "informative_markers_sha256" {print $2}' \
    "$MARKER_SETTINGS")
  [[ "${#marker_hashes[@]}" -eq 1 ]] ||
    die "Expected one informative_markers_sha256 in $MARKER_SETTINGS"
  [[ "$(sha256sum "$MARKERS" | cut -d ' ' -f1)" == "${marker_hashes[0]}" ]] ||
    die "Marker table differs from the set frozen by stage 5"

  # Require the manifest used for source genotyping.
  mapfile -t manifest_hashes < <(awk -F '\t' '$1 ~ /\/src2\/spikein\/samples.tsv$/ {print $2}' \
    "$SOURCE_SETTINGS")
  [[ "${#manifest_hashes[@]}" -eq 1 ]] ||
    die "Expected one samples.tsv hash in $SOURCE_SETTINGS"
  [[ "$(sha256sum "$SAMPLES_TSV" | cut -d ' ' -f1)" == "${manifest_hashes[0]}" ]] ||
    die "Sample manifest differs from the source-genotyping manifest"

  # Check the marker columns used to build the bam-readcount site list.
  IFS=$'\t' read -r _ chromosome position ref _ < "$MARKERS"
  [[ "$chromosome" == "chromosome" && "$position" == "position" && "$ref" == "ref" ]] ||
    die "Unexpected marker-table columns"
}

set_output_paths() {
  # Keep all new files below results2/spikein and refuse existing outputs.
  OUT_ROOT="$(realpath -m "$OUT_ROOT")"
  case "$OUT_ROOT" in
    "$REPO_ROOT/results2/spikein/"*) ;;
    *) die "OUT_ROOT must be a directory below results2/spikein" ;;
  esac
  [[ ! -e "$OUT_ROOT" ]] || die "Output directory already exists: $OUT_ROOT"
  READCOUNTS_DIR="$OUT_ROOT/readcounts"
  SITES="$OUT_ROOT/sites.tsv"
}

read_mixture_samples() {
  # Require exactly one high and one low mixture, using roles rather than filenames.
  mixture_samples="$(awk -F '\t' '$1 == "high" || $1 == "low" {
    print $1 "\t" $2 "\t" $3
  }' "$SAMPLES_TSV")"
  [[ "$(cut -f1 <<< "$mixture_samples" | sort)" == $'high\nlow' ]] ||
    die "Expected one high and one low mixture in $SAMPLES_TSV"
  [[ "$(cut -f2 <<< "$mixture_samples" | sort -u | wc -l)" -eq 2 ]] ||
    die "High and low mixtures must have distinct sample IDs"
  while IFS=$'\t' read -r role sample bam; do
    [[ -s "$bam" && -s "$bam.bai" ]] || die "Missing BAM or canonical BAI for $role: $bam"
  done <<< "$mixture_samples"
  high_bam="$(awk -F '\t' '$1 == "high" {print $3}' "$SAMPLES_TSV")"
  low_bam="$(awk -F '\t' '$1 == "low" {print $3}' "$SAMPLES_TSV")"
  [[ "$(realpath "$high_bam")" != "$(realpath "$low_bam")" ]] ||
    die "High and low mixtures must use distinct BAMs"
}

build_sites() {
  # bam-readcount takes 1-based inclusive CHROM/POS/POS, not BED coordinates.
  sites="$(awk -F '\t' 'BEGIN {OFS = "\t"} NR > 1 {print $2, $3, $3}' \
    "$MARKERS" | sort -u)"
  marker_rows="$(awk 'END {print NR - 1}' "$MARKERS")"
  [[ "$marker_rows" -gt 0 && "$(wc -l <<< "$sites")" -eq "$marker_rows" ]] ||
    die "Marker sites must be non-empty and unique"
}

validate_counts() {
  local role="$1"
  local counts="$2"
  # Require every requested site once, with the frozen reference base.
  diff -u <(cut -f1,2 "$SITES") <(cut -f1,2 "$counts" | sort) ||
    die "Read-count sites differ from requested sites for $role"
  awk -F '\t' 'NR == FNR {ref[$2 FS $3] = $4; next}
    toupper($3) != ref[$1 FS $2] {exit 1}' "$MARKERS" "$counts" ||
    die "Read-count reference differs from the frozen marker set for $role"
  awk -v cap="$MAX_DEPTH" '$4 >= cap {exit 1}' "$counts" ||
    die "Reported depth reached the bam-readcount cap for $role"
}

main() {
  # Complete all preflight checks before creating output files.
  check_inputs "$@"
  check_frozen_inputs
  set_output_paths
  read_mixture_samples
  build_sites

  # A dry run creates nothing. Actual runs retain commands and raw tool messages.
  if [[ "$DRY_RUN" == "0" ]]; then
    mkdir -p "$(dirname "$OUT_ROOT")"
    mkdir "$OUT_ROOT"
    mkdir "$READCOUNTS_DIR"
    exec > >(tee "$OUT_ROOT/run.log") 2>&1
    git rev-parse HEAD
    sha256sum "$0" "$MARKERS" "$MARKER_SETTINGS" "$SAMPLES_TSV" \
      "$SOURCE_SETTINGS" "$REF_FASTA"
    printf '%s\n' "$mixture_samples"
  fi
  run_to_file "$SITES" printf '%s\n' "$sites"

  # Count directly from each original BAM. bam-readcount excludes flagged duplicates.
  # Explicit quality and depth settings preserve the established convention.
  while IFS=$'\t' read -r role sample bam; do
    echo "=== $role: $sample ==="
    counts="$READCOUNTS_DIR/$sample.txt"
    run_to_file "$counts" bam-readcount -q 0 -b 0 -d "$MAX_DEPTH" \
      -f "$REF_FASTA" -l "$SITES" "$bam"
    [[ "$DRY_RUN" == "1" ]] || validate_counts "$role" "$counts"
  done <<< "$mixture_samples"

  # Leave raw counts for the separate parser and recovery stage.
  echo "=== Mixture read-count collection finished (DRY_RUN=$DRY_RUN) ==="
}

# Run the stage after all functions have been defined.
main "$@"
