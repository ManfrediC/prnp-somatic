#!/usr/bin/env bash
# Stop on command failures, unset variables and failed pipeline steps.
set -euo pipefail

# ------------------------------------------------------------
# Apply the selected orientation and call filters to both mixtures.
# Publish results only after both filtered VCFs pass validation.
# ------------------------------------------------------------

# Resolve project paths from this script, regardless of the working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd -P)"
cd "$REPO_ROOT"

# Use the final stage-8 callset and the marker set frozen before mixture analysis.
SAMPLES_TSV="src2/spikein/samples.tsv"
MARKERS="results2/spikein/markers/informative_markers.tsv"
REF_FASTA="resources/references/snv/chr2_chr4_chr20.fasta"
INPUT_ROOT="${INPUT_ROOT:-results2/spikein/mutect2}"
OUT_ROOT="${OUT_ROOT:-results2/spikein/filtermutectcalls}"
DRY_RUN="${DRY_RUN:-1}"
JAVA_MEM_GB="${JAVA_MEM_GB:-8}"
LOG_REDIRECTED=0
export LC_ALL=C

# Report failed checks and print commands before execution.
die() { echo "ERROR: $*" >&2; exit 1; }
run() {
  printf '+ '
  printf '%q ' "$@"
  printf '\n'
  [[ "$DRY_RUN" == "1" ]] || "$@"
}

setting_value() {
  awk -F '\t' -v key="$1" '$1 == key {print $2}' "$INPUT_ROOT/run_settings.tsv"
}

set_paths() {
  # Keep all inputs and outputs below results2/spikein.
  INPUT_ROOT="$(realpath "$INPUT_ROOT")"
  OUT_ROOT="$(realpath -m "$OUT_ROOT")"
  case "$INPUT_ROOT" in
    "$REPO_ROOT/results2/spikein/"*) ;;
    *) die "INPUT_ROOT must be below results2/spikein" ;;
  esac
  case "$OUT_ROOT" in
    "$REPO_ROOT/results2/spikein/"*) ;;
    *) die "OUT_ROOT must be below results2/spikein" ;;
  esac
  [[ "$INPUT_ROOT" != "$OUT_ROOT" ]] || die "Input and output directories must differ"
  [[ ! -e "$OUT_ROOT" ]] || die "Output directory already exists: $OUT_ROOT"

  # Stage beside the final directory so publication is one move.
  STAGING="$OUT_ROOT.tmp.$$"
  case "$STAGING" in
    "$REPO_ROOT/results2/spikein/"*.tmp.*) ;;
    *) die "Unsafe temporary output path: $STAGING" ;;
  esac
  FILTERED_DIR="$STAGING/filtered"
  ORIENTATION_DIR="$STAGING/orientation"
  WORK_DIR="$STAGING/work"
  GATK_JAVA_OPTIONS="-Xmx${JAVA_MEM_GB}g -Djava.io.tmpdir=$WORK_DIR/tmp"
}

read_mixture_samples() {
  # Require exactly one high and one low mixture.
  mixture_samples="$(awk -F '\t' '$1 == "high" || $1 == "low" {
    print $1 "\t" $2
  }' "$SAMPLES_TSV")"
  [[ "$(cut -f1 <<< "$mixture_samples" | sort)" == $'high\nlow' ]] ||
    die "Expected one high and one low mixture in $SAMPLES_TSV"
  [[ "$(cut -f2 <<< "$mixture_samples" | sort -u | wc -l)" -eq 2 ]] ||
    die "High and low mixtures must have distinct sample IDs"
}

check_inputs() {
  # Check tools, fixed inputs and the selected caller settings.
  [[ "$#" -eq 0 ]] || die "Use DRY_RUN, JAVA_MEM_GB, INPUT_ROOT and OUT_ROOT, not positional arguments"
  [[ "$DRY_RUN" == "0" || "$DRY_RUN" == "1" ]] || die "DRY_RUN must be 0 or 1"
  [[ "$JAVA_MEM_GB" =~ ^[1-9][0-9]*$ ]] || die "JAVA_MEM_GB must be a positive integer"
  for tool in gatk bcftools sha256sum; do
    command -v "$tool" >/dev/null || die "$tool not in PATH (activate prnp-spikein)"
  done
  for input in "$SAMPLES_TSV" "$MARKERS" "$REF_FASTA" "$REF_FASTA.fai" \
    "$INPUT_ROOT/run_settings.tsv"; do
    [[ -s "$input" ]] || die "Missing input: $input"
  done

  [[ "$(setting_value initial_tumor_lod)" == "0" ]] || die "Stage 8 did not use initial tumour LOD 0"
  [[ "$(setting_value max_population_af)" == "1.0" ]] || die "Stage 8 did not use population AF 1.0"
  [[ "$(setting_value max_reads_per_alignment_start)" == "0" ]] || die "Stage 8 did not disable the read cap"
  emission_lod="$(setting_value tumor_lod_to_emit)"
  [[ -z "$emission_lod" || "$emission_lod" == "3" ]] || die "Stage 8 did not use emission LOD 3"
  activity_threshold="$(setting_value active_probability_threshold)"
  [[ -z "$activity_threshold" || "$activity_threshold" == "0.002" ]] ||
    die "Stage 8 did not use the default activity-probability threshold"
}

check_callset() {
  # Require the frozen marker hash and complete stage-8 outputs for both mixtures.
  marker_hash="$(sha256sum "$MARKERS" | cut -d ' ' -f1)"
  [[ "$(setting_value informative_markers_sha256)" == "$marker_hash" ]] ||
    die "Stage 8 marker hash differs from the frozen marker table"

  while IFS=$'\t' read -r role sample; do
    raw_vcf="$INPUT_ROOT/vcf/$sample.raw.vcf.gz"
    for input in "$raw_vcf" "$raw_vcf.tbi" "$raw_vcf.stats" \
      "$INPUT_ROOT/f1r2/$sample.f1r2.tar.gz"; do
      [[ -s "$input" ]] || die "Missing stage-8 input for $role: $input"
    done
    bcftools view -Ov "$raw_vcf" >/dev/null || die "Unreadable VCF: $raw_vcf"
    [[ "$(bcftools query -l "$raw_vcf")" == "$sample" ]] ||
      die "Unexpected sample in $raw_vcf"
  done <<< "$mixture_samples"
}

cleanup_partial_output() {
  # Remove only this run's temporary directory after a failure.
  local status=$?
  trap - EXIT
  if [[ "$LOG_REDIRECTED" == "1" ]]; then
    exec 1>&3 2>&4
    exec 3>&- 4>&-
    LOG_REDIRECTED=0
  fi
  if [[ "$status" -ne 0 && -d "${STAGING:-}" ]]; then
    rm -rf -- "$STAGING"
  fi
  exit "$status"
}

filter_sample() {
  local role="$1"
  local sample="$2"
  local raw_vcf="$INPUT_ROOT/vcf/$sample.raw.vcf.gz"
  local orientation="$ORIENTATION_DIR/$sample.orientation.tar.gz"
  local filtered="$FILTERED_DIR/$sample.filtered.vcf.gz"

  # Learn the orientation model, then apply the selected clustered-events limit.
  echo "=== $role: $sample ==="
  run gatk --java-options "$GATK_JAVA_OPTIONS" LearnReadOrientationModel \
    -I "$INPUT_ROOT/f1r2/$sample.f1r2.tar.gz" \
    -O "$orientation" \
    --num-em-iterations 50
  run gatk --java-options "$GATK_JAVA_OPTIONS" FilterMutectCalls \
    -R "$REF_FASTA" \
    -V "$raw_vcf" \
    --orientation-bias-artifact-priors "$orientation" \
    --stats "$raw_vcf.stats" \
    --max-events-in-region 5 \
    -O "$filtered"

  # Require a readable, indexed VCF with the expected sample and record count.
  if [[ "$DRY_RUN" == "0" ]]; then
    [[ -s "$orientation" && -s "$filtered" && -s "$filtered.tbi" ]] ||
      die "Missing filtered output for $sample"
    bcftools view -Ov "$filtered" >/dev/null || die "Unreadable VCF: $filtered"
    [[ "$(bcftools query -l "$filtered")" == "$sample" ]] ||
      die "Unexpected sample in $filtered"
    [[ "$(bcftools view -H "$raw_vcf" | wc -l)" -eq "$(bcftools view -H "$filtered" | wc -l)" ]] ||
      die "FilterMutectCalls changed the record count for $sample"
  fi
}

write_settings() {
  # Record the selected filter setting and the matching stage-8 input.
  {
    printf 'key\tvalue\n'
    printf 'input_root\t%s\nreference\t%s\n' "$INPUT_ROOT" "$REF_FASTA"
    printf 'orientation_em_iterations\t50\nmax_events_in_region\t5\n'
    printf 'java_mem_gb\t%s\ninformative_markers_sha256\t%s\n' "$JAVA_MEM_GB" "$marker_hash"
    printf 'git_commit\t%s\n' "$(git rev-parse HEAD)"
  } > "$STAGING/run_settings.tsv"
}

main() {
  local role sample
  set_paths
  read_mixture_samples
  check_inputs "$@"
  check_callset

  # An actual run writes only to the temporary output tree.
  if [[ "$DRY_RUN" == "0" ]]; then
    mkdir -p "$FILTERED_DIR" "$ORIENTATION_DIR" "$WORK_DIR/tmp"
    trap cleanup_partial_output EXIT
    exec 3>&1 4>&2
    exec > "$STAGING/run.log" 2>&1
    LOG_REDIRECTED=1
    write_settings
    printf '%s\n' "$mixture_samples"
    gatk --version
    bcftools --version
  fi

  # Filter the high and low mixtures independently.
  while IFS=$'\t' read -r role sample; do
    filter_sample "$role" "$sample"
  done <<< "$mixture_samples"

  # Remove temporary work, hash the validated files and publish the result.
  if [[ "$DRY_RUN" == "0" ]]; then
    rm -rf -- "$WORK_DIR"
    exec 1>&3 2>&4
    exec 3>&- 4>&-
    LOG_REDIRECTED=0
    (
      cd "$STAGING"
      find . -type f ! -name sha256.txt -printf '%P\0' | sort -z | xargs -0 sha256sum
    ) > "$STAGING/sha256.txt"
    mv "$STAGING" "$OUT_ROOT"
  fi

  echo "=== Mixture filtering finished (DRY_RUN=$DRY_RUN) ==="
}

# Run the stage after all functions have been defined.
main "$@"
