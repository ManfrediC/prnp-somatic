#!/usr/bin/env bash
# Stop on command failures, unset variables and failed pipeline steps.
set -euo pipefail

# ------------------------------------------------------------
# Collect fixed-site read counts from the two pure sources.
# ------------------------------------------------------------

# Resolve project paths from this script, regardless of the working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd -P)"
cd "$REPO_ROOT"

# Use the revised candidate table; mixture evidence cannot enter this stage.
SAMPLES_TSV="src2/spikein/samples.tsv"
MARKERS="results2/spikein/discovery/candidates_revised/candidate_markers.tsv"
REF_FASTA="resources/references/snv/chr2_chr4_chr20.fasta"
OUT_ROOT="${OUT_ROOT:-results2/spikein/readcount_qc/pure}"
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

# Check the tool and inputs before creating any output.
[[ "$#" -eq 0 ]] || die "Use DRY_RUN and OUT_ROOT, not positional arguments"
[[ "$DRY_RUN" == "0" || "$DRY_RUN" == "1" ]] || die "DRY_RUN must be 0 or 1"
command -v bam-readcount >/dev/null || die "bam-readcount not in PATH (activate prnp-spikein)"
for input in Makefile "$SAMPLES_TSV" "$MARKERS" "$REF_FASTA" "$REF_FASTA.fai"; do
  [[ -s "$input" ]] || die "Missing input: $input"
done

# Keep all new files below results2/spikein and refuse existing outputs.
OUT_ROOT="$(realpath -m "$OUT_ROOT")"
case "$OUT_ROOT" in
  "$REPO_ROOT/results2/spikein/"*) ;;
  *) die "OUT_ROOT must be a directory below results2/spikein" ;;
esac
[[ ! -e "$OUT_ROOT" ]] || die "Output directory already exists: $OUT_ROOT"
READCOUNTS_DIR="$OUT_ROOT/readcounts"
SITES="$OUT_ROOT/sites.tsv"

# Require exactly one donor and one WT row, using roles rather than filenames.
pure_samples="$(awk -F '\t' '$1 == "pure_donor" || $1 == "pure_wt" {
  print $1 "\t" $2 "\t" $3
}' "$SAMPLES_TSV")"
[[ "$(cut -f1 <<< "$pure_samples" | sort)" == $'pure_donor\npure_wt' ]] ||
  die "Expected one pure_donor and one pure_wt in $SAMPLES_TSV"
while IFS=$'\t' read -r role sample bam; do
  [[ -s "$bam" && -s "$bam.bai" ]] || die "Missing BAM or canonical BAI for $role: $bam"
done <<< "$pure_samples"

# bam-readcount takes 1-based inclusive CHROM/POS/POS, not BED coordinates.
sites="$(awk -F '\t' 'BEGIN {OFS = "\t"}
  NR == 1 {if ($2 != "chromosome" || $3 != "position") exit 1; next}
  {print $2, $3, $3}
' "$MARKERS" | sort -u)"
[[ -n "$sites" ]] || die "No candidate sites in $MARKERS (A117V must be included)"

# A dry run creates nothing. Actual runs retain commands and raw tool messages.
if [[ "$DRY_RUN" == "0" ]]; then
  mkdir -p "$(dirname "$OUT_ROOT")"
  mkdir "$OUT_ROOT"
  mkdir "$READCOUNTS_DIR"
  exec > >(tee "$OUT_ROOT/run.log") 2>&1
  git rev-parse HEAD
  sha256sum "$MARKERS" "$SAMPLES_TSV"
fi
run_to_file "$SITES" printf '%s\n' "$sites"

# Count directly from each original BAM. bam-readcount excludes flagged duplicates.
# Explicit quality and depth settings preserve the established convention.
while IFS=$'\t' read -r role sample bam; do
  echo "=== $role: $sample ==="
  counts="$READCOUNTS_DIR/$sample.txt"
  run_to_file "$counts" bam-readcount -q 0 -b 0 -d "$MAX_DEPTH" \
    -f "$REF_FASTA" -l "$SITES" "$bam"

  # Require every requested site exactly once.
  if [[ "$DRY_RUN" == "0" ]]; then
    diff -u <(cut -f1,2 "$SITES") <(cut -f1,2 "$counts" | sort) ||
      die "Read-count sites differ from requested sites for $role"
    awk -v cap="$MAX_DEPTH" '$4 >= cap {exit 1}' "$counts" ||
      die "Reported depth reached the bam-readcount cap for $role"
  fi
done <<< "$pure_samples"

# Leave raw counts for the separate parser and pure-source validation step.
echo "=== Pure-source read-count collection finished (DRY_RUN=$DRY_RUN) ==="
