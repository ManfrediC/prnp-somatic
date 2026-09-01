#!/usr/bin/env bash
# Stop on command failures, unset variables and failed pipeline steps.
set -euo pipefail

# ------------------------------------------------------------
# Run the established tumour-only Mutect2 call on both mixtures.
# Marker recovery is assessed separately against the frozen set.
# ------------------------------------------------------------

# Resolve project paths from this script, regardless of the working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd -P)"
cd "$REPO_ROOT"

# Use the authoritative sample roles and the marker set frozen before mixture analysis.
SAMPLES_TSV="src2/spikein/samples.tsv"
MARKERS="results2/spikein/markers/informative_markers.tsv"
MARKER_SETTINGS="results2/spikein/markers/run_settings.tsv"
REF_FASTA="resources/references/snv/chr2_chr4_chr20.fasta"
INTERVALS="resources/intervals/capture_targets.interval_list"
GNOMAD_VCF="resources/population/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
PON_VCF="results/mutect2_controls_pon/panel_of_normals/CJD_controls_PoN.vcf.gz"
OUT_ROOT="${OUT_ROOT:-results2/spikein/mutect2}"
DRY_RUN="${DRY_RUN:-1}"
JAVA_MEM_GB="${JAVA_MEM_GB:-8}"
export LC_ALL=C

# Report failed checks and print commands before execution.
die() { echo "ERROR: $*" >&2; exit 1; }
run() {
  printf '+ '
  printf '%q ' "$@"
  printf '\n'
  [[ "$DRY_RUN" == "1" ]] || "$@"
}

check_inputs() {
  # Check settings, tools and fixed project inputs before creating output.
  [[ "$#" -eq 0 ]] || die "Use DRY_RUN, JAVA_MEM_GB and OUT_ROOT, not positional arguments"
  [[ "$DRY_RUN" == "0" || "$DRY_RUN" == "1" ]] || die "DRY_RUN must be 0 or 1"
  [[ "$JAVA_MEM_GB" =~ ^[1-9][0-9]*$ ]] || die "JAVA_MEM_GB must be a positive integer"
  for tool in gatk bcftools; do
    command -v "$tool" >/dev/null || die "$tool not in PATH (activate prnp-spikein)"
  done
  for input in Makefile "$SAMPLES_TSV" "$MARKERS" "$MARKER_SETTINGS" \
    "$REF_FASTA" "$REF_FASTA.fai" "$INTERVALS" "$GNOMAD_VCF" \
    "$GNOMAD_VCF.tbi" "$PON_VCF" "$PON_VCF.tbi"; do
    [[ -s "$input" ]] || die "Missing input: $input"
  done
}

check_frozen_markers() {
  # Require the exact marker table fixed by stage 5 before mixtures were examined.
  mapfile -t marker_hashes < <(awk -F '\t' '$1 == "informative_markers_sha256" {print $2}' \
    "$MARKER_SETTINGS")
  [[ "${#marker_hashes[@]}" -eq 1 ]] ||
    die "Expected one informative_markers_sha256 in $MARKER_SETTINGS"
  marker_hash="$(sha256sum "$MARKERS" | cut -d ' ' -f1)"
  [[ "$marker_hash" == "${marker_hashes[0]}" ]] ||
    die "Marker table differs from the set frozen by stage 5"
}

set_output_paths() {
  # Keep all new files below results2/spikein and refuse existing outputs.
  OUT_ROOT="$(realpath -m "$OUT_ROOT")"
  [[ "$OUT_ROOT" != *[[:space:]]* ]] || die "OUT_ROOT must not contain whitespace"
  case "$OUT_ROOT" in
    "$REPO_ROOT/results2/spikein/"*) ;;
    *) die "OUT_ROOT must be a directory below results2/spikein" ;;
  esac
  [[ ! -e "$OUT_ROOT" ]] || die "Output directory already exists: $OUT_ROOT"
  VCF_DIR="$OUT_ROOT/vcf"
  F1R2_DIR="$OUT_ROOT/f1r2"
  TMP_DIR="$OUT_ROOT/work/tmp"
  GATK_JAVA_OPTIONS="-Xmx${JAVA_MEM_GB}g -Djava.io.tmpdir=$TMP_DIR"
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
  mapfile -t mixture_bams < <(cut -f3 <<< "$mixture_samples")
  [[ "$(realpath "${mixture_bams[0]}")" != "$(realpath "${mixture_bams[1]}")" ]] ||
    die "High and low mixtures must use distinct BAMs"
}

run_mutect2() {
  local role="$1"
  local sample="$2"
  local bam="$3"
  local out_vcf="$VCF_DIR/$sample.raw.vcf.gz"
  local out_f1r2="$F1R2_DIR/$sample.f1r2.tar.gz"

  # Match the established PoN-enabled call across the full capture panel.
  echo "=== $role: $sample ==="
  run gatk --java-options "$GATK_JAVA_OPTIONS" Mutect2 \
    -R "$REF_FASTA" \
    -I "$bam" \
    -tumor-sample "$sample" \
    --panel-of-normals "$PON_VCF" \
    --germline-resource "$GNOMAD_VCF" \
    --af-of-alleles-not-in-resource 0.0000025 \
    --initial-tumor-lod 0 \
    --max-population-af 1.0 \
    --max-reads-per-alignment-start 0 \
    --f1r2-tar-gz "$out_f1r2" \
    --intervals "$INTERVALS" \
    --tmp-dir "$TMP_DIR" \
    -O "$out_vcf"

  # Require a readable, indexed single-sample VCF and orientation data.
  if [[ "$DRY_RUN" == "0" ]]; then
    [[ -s "$out_vcf" && -s "$out_vcf.tbi" && -s "$out_f1r2" ]] ||
      die "Missing Mutect2 output for $sample"
    bcftools view -Ov "$out_vcf" >/dev/null || die "Unreadable VCF: $out_vcf"
    [[ "$(bcftools query -l "$out_vcf")" == "$sample" ]] ||
      die "Unexpected VCF sample for $sample"
  fi
}

main() {
  # Complete all preflight checks before creating output files.
  check_inputs "$@"
  check_frozen_markers
  set_output_paths
  read_mixture_samples

  # A dry run creates nothing. Actual runs retain commands and raw tool messages.
  if [[ "$DRY_RUN" == "0" ]]; then
    mkdir -p "$(dirname "$OUT_ROOT")"
    mkdir "$OUT_ROOT"
    mkdir "$VCF_DIR" "$F1R2_DIR"
    mkdir -p "$TMP_DIR"
    exec > >(tee "$OUT_ROOT/run.log") 2>&1
    {
      printf 'key\tvalue\n'
      printf 'reference\t%s\nintervals\t%s\n' "$REF_FASTA" "$INTERVALS"
      printf 'germline_resource\t%s\npanel_of_normals\t%s\n' "$GNOMAD_VCF" "$PON_VCF"
      printf 'af_of_alleles_not_in_resource\t0.0000025\n'
      printf 'initial_tumor_lod\t0\n'
      printf 'max_population_af\t1.0\n'
      printf 'max_reads_per_alignment_start\t0\n'
      printf 'java_mem_gb\t%s\ninformative_markers_sha256\t%s\n' "$JAVA_MEM_GB" "$marker_hash"
      printf 'git_commit\t%s\n' "$(git rev-parse HEAD)"
    } > "$OUT_ROOT/run_settings.tsv"
    printf '%s\n' "$mixture_samples"
    gatk --version
    bcftools --version
  fi

  # Call both mixtures independently; no marker or recovery decision is made here.
  while IFS=$'\t' read -r role sample bam; do
    run_mutect2 "$role" "$sample" "$bam"
  done <<< "$mixture_samples"

  echo "=== Mixture Mutect2 calling finished (DRY_RUN=$DRY_RUN) ==="
}

# Run the stage after all functions have been defined.
main "$@"
