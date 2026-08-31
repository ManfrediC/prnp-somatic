#!/usr/bin/env bash
# Stop on command failures, unset variables and failed pipeline steps.
set -euo pipefail

# ------------------------------------------------------------
# Discover germline variants in the two pure sources.
# Check all four samples; call only the donor and WT.
# Marker selection and mixture counting follow separately.
# ------------------------------------------------------------

# Default to checking inputs and printing commands without creating files.
DRY_RUN="${DRY_RUN:-1}"

# Set the Java heap for each sequential GATK process.
JAVA_MEM_GB="${JAVA_MEM_GB:-8}"

# Resolve project paths from this script, regardless of the working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd -P)"
# Require the repository Makefile at this location.
[[ -f "$REPO_ROOT/Makefile" ]] || { echo "ERROR: repository Makefile not found" >&2; exit 1; }
cd "$REPO_ROOT"

# Use explicit sample roles and a separate directory for new outputs.
SAMPLES_TSV="${SAMPLES_TSV:-src2/spikein/samples.tsv}"
OUT_ROOT="${OUT_ROOT:-results2/spikein/discovery}"

# Use the project reference and the complete capture panel.
REF_FASTA="resources/references/snv/chr2_chr4_chr20.fasta"
REF_DICT="resources/references/snv/chr2_chr4_chr20.dict"
INTERVALS="resources/intervals/capture_targets.interval_list"

# Check available tools and report failed checks.
have() { command -v "$1" >/dev/null 2>&1; }
die() { echo "ERROR: $*" >&2; exit 1; }

# Print quoted commands and execute them only when DRY_RUN is 0.
run() {
  printf '+ '
  printf '%q ' "$@"
  printf '\n'
  [[ "$DRY_RUN" == "1" ]] || "$@"
}

# Extract contig names and lengths in header order.
sequence_dictionary() {
  awk -F '\t' '$1 == "@SQ" {
    name = ""; size = ""
    for (i = 2; i <= NF; i++) {
      if ($i ~ /^SN:/) name = substr($i, 4)
      if ($i ~ /^LN:/) size = substr($i, 4)
    }
    sub(/\r$/, "", size)
    print name "\t" size
  }'
}

# Check execution settings before reading inputs.
[[ "$#" -eq 0 ]] || die "Use environment settings, not positional arguments"
[[ "$DRY_RUN" == "0" || "$DRY_RUN" == "1" ]] || die "DRY_RUN must be 0 or 1"
[[ "$JAVA_MEM_GB" =~ ^[1-9][0-9]*$ ]] || die "JAVA_MEM_GB must be a positive integer"

# Require the caller and the tools used to check inputs and outputs.
for tool in gatk samtools bcftools; do
  have "$tool" || die "$tool not in PATH (activate prnp-spikein)"
done

# Resolve the output path before checking its location.
OUT_ROOT="$(realpath -m "$OUT_ROOT")"
# The installed GATK launcher splits --java-options on whitespace.
[[ "$OUT_ROOT" != *[[:space:]]* ]] || die "OUT_ROOT must not contain whitespace (GATK Java options)"
# Keep outputs below results2/spikein and refuse existing directories.
case "$OUT_ROOT" in
  "$REPO_ROOT/results2/spikein/"*) ;;
  *) die "OUT_ROOT must be a directory below results2/spikein" ;;
esac
[[ ! -e "$OUT_ROOT" ]] || die "Output directory already exists: $OUT_ROOT; use a fresh directory"
# Require readable metadata, reference files and capture intervals.
for path in "$SAMPLES_TSV" "$REF_FASTA" "$REF_FASTA.fai" "$REF_DICT" "$INTERVALS"; do
  [[ -r "$path" && -s "$path" ]] || die "Missing, empty or unreadable input: $path"
done
SAMPLES_TSV="$(realpath -e "$SAMPLES_TSV")"

# Require matching contig names, lengths and order in reference metadata.
ref_contigs="$(sequence_dictionary < "$REF_DICT")"
[[ -n "$ref_contigs" ]] || die "Empty reference sequence dictionary"
[[ "$ref_contigs" == "$(awk -F '\t' '{print $1 "\t" $2}' "$REF_FASTA.fai")" ]] || die "Reference dictionary and FAI disagree"
[[ "$ref_contigs" == "$(sequence_dictionary < "$INTERVALS")" ]] || die "Capture and reference dictionaries disagree"
# Check that capture intervals lie within their reference contigs.
awk -F '\t' '
  NR == FNR { size[$1] = $2; next }
  /^@/ { next }
  {
    if (NF < 5 || !($1 in size) || $2 !~ /^[0-9]+$/ ||
        $3 !~ /^[0-9]+$/ || $2 < 1 || $3 < $2 || $3 > size[$1]) exit 1
    count++
  }
  END { if (!count) exit 1 }
' "$REF_FASTA.fai" "$INTERVALS" || die "Invalid capture interval coordinates"
# Use one position to check index readability later.
first_contig="$(head -n 1 "$REF_FASTA.fai" | cut -f 1)"

# Require the expected five-column manifest.
IFS= read -r header < "$SAMPLES_TSV"
[[ "${header%$'\r'}" == $'role\tsample_id\tbam\tdonor_fraction\tfraction_provenance' ]] || die "Unexpected sample manifest header"

# Store checked inputs by role, independently of filenames and row order.
declare -A bam_by_role=() sample_by_role=() index_by_role=()

# Check all four BAMs; mixture data do not enter marker discovery.
while IFS=$'\t' read -r role sample bam fraction provenance extra; do
  # Require complete rows, recognised roles and safe sample names.
  [[ -n "$role" && -n "$sample" && -n "$bam" && -n "$fraction" && -n "$provenance" && -z "$extra" ]] || die "Expected five populated manifest columns"
  case "$role" in
    pure_donor|pure_wt|high|low) ;;
    *) die "Unexpected sample role: $role" ;;
  esac
  [[ -z "${bam_by_role[$role]+set}" ]] || die "Duplicate sample role: $role"
  [[ "$sample" =~ ^[A-Za-z0-9_-]+$ ]] || die "Invalid sample ID: $sample"
  # Require distinct BAM files, including when paths are hard links.
  [[ -r "$bam" && -s "$bam" ]] || die "Missing, empty or unreadable BAM for $role: $bam"
  bam="$(realpath -e "$bam")"
  for previous in "${bam_by_role[@]}"; do
    [[ ! "$bam" -ef "$previous" ]] || die "BAM file reused across roles: $bam"
  done
  # Require distinct sample IDs.
  for previous in "${sample_by_role[@]}"; do
    [[ "$sample" != "$previous" ]] || die "Sample ID reused across roles: $sample"
  done

  # Check the header and end marker, without scanning every BAM record.
  samtools quickcheck -v "$bam" || die "BAM quickcheck failed for $role: $bam"

  # Require matching reference contigs and declared coordinate order.
  bam_header="$(samtools view -H "$bam")"
  [[ "$(sequence_dictionary <<< "$bam_header")" == "$ref_contigs" ]] || die "BAM/reference contigs disagree for $role"
  awk -F '\t' '$1 == "@HD" {for (i = 2; i <= NF; i++) if ($i == "SO:coordinate") ok = 1} END {exit !ok}' <<< "$bam_header" || die "BAM is not marked coordinate-sorted: $bam"

  # Require every read group to name the sample in the manifest.
  sample_names="$(awk -F '\t' '$1 == "@RG" {
    sample = ""
    for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) sample = substr($i, 4)
    if (sample == "") exit 1
    print sample
  }' <<< "$bam_header" | LC_ALL=C sort -u)" || die "Read group without SM for $role"
  [[ "$sample_names" == "$sample" ]] || die "BAM SM does not match manifest for $role: expected $sample, found $sample_names"

  # Read the canonical .bam.bai beside the BAM; leave archived indexes alone.
  index="$bam.bai"
  [[ -e "$index" ]] || die "Missing existing BAI for $role: $index; this stage does not create indexes"
  [[ -r "$index" && -s "$index" ]] || die "Unreadable or empty index: $index"

  # Check a regional query and discard its count.
  samtools view -c -X "$bam" "$index" "$first_contig:1-1" >/dev/null || die "Cannot read BAM index for $role: $index"

  # Save the checked inputs by sample role.
  bam_by_role[$role]="$bam"
  sample_by_role[$role]="$sample"
  index_by_role[$role]="$index"
done < <(tail -n +2 "$SAMPLES_TSV")

# Require all four roles before proceeding.
[[ "${#bam_by_role[@]}" -eq 4 ]] || die "Manifest must contain exactly pure_donor, pure_wt, high and low"

# Prepare output paths; dry runs print commands without creating files.
WORK_DIR="$OUT_ROOT/work"
TMP_DIR="$WORK_DIR/tmp"
# Share Java settings across the four GATK commands.
GATK_JAVA_OPTIONS="-Xmx${JAVA_MEM_GB}g -Djava.io.tmpdir=$TMP_DIR"
run mkdir -p "$(dirname "$OUT_ROOT")"
run mkdir "$OUT_ROOT"
run mkdir -p "$TMP_DIR"

# Log actual runs and keep temporary files below the output directory.
if [[ "$DRY_RUN" == "0" ]]; then
  exec > >(tee "$OUT_ROOT/run.log") 2>&1
  export TMPDIR="$TMP_DIR" TMP="$TMP_DIR" TEMP="$TMP_DIR"
fi

# Report checked inputs and settings.
echo "== Pure donor + WT germline discovery =="
echo "Preflight passed for all four manifest roles"
echo "Repo root:      $REPO_ROOT"
echo "Samples TSV:    $SAMPLES_TSV"
echo "Output root:    $OUT_ROOT"
echo "REF_FASTA:      $REF_FASTA"
echo "INTERVALS:      $INTERVALS"
echo "JAVA_MEM_GB:    $JAVA_MEM_GB"
echo "DRY_RUN:        $DRY_RUN"
echo "Only pure_donor and pure_wt enter the callers"

# Record provenance only during an actual run.
if [[ "$DRY_RUN" == "0" ]]; then
  # Record caller settings and source inputs.
  {
    printf 'key\tvalue\n'
    printf 'reference\t%s\nintervals\t%s\n' "$REF_FASTA" "$INTERVALS"
    printf 'java_mem_gb\t%s\n' "$JAVA_MEM_GB"
    printf 'ploidy\t2\nhc_max_reads_per_alignment_start\t0\nhc_native_pair_hmm_threads\t1\n'
    printf 'git_commit\t%s\n' "$(git rev-parse HEAD)"
    # Hash metadata, code, reference files and the two pure BAM/index pairs.
    for path in \
      "$SAMPLES_TSV" "$SCRIPT_DIR/1_sources_haplotypecaller.sh" \
      "$REPO_ROOT/src/pipelines/1_controls_mutect2_no_pon.sh" \
      "$REF_FASTA" "$REF_FASTA.fai" "$REF_DICT" "$INTERVALS" \
      "${bam_by_role[pure_donor]}" "${index_by_role[pure_donor]}" \
      "${bam_by_role[pure_wt]}" "${index_by_role[pure_wt]}"; do
      printf 'sha256:%s\t%s\n' "$path" "$(sha256sum "$path" | cut -d ' ' -f 1)"
    done
  } > "$OUT_ROOT/run_settings.tsv"
  # Log sample roles, Git state and tool versions.
  cat "$SAMPLES_TSV"
  GIT_OPTIONAL_LOCKS=0 git status --short
  gatk --version
  samtools --version
  bcftools --version
fi

# Call the two pure sources sequentially, using their explicit roles.
gvcfs=()
for role in pure_donor pure_wt; do
  # Select checked inputs and keep the gVCF path for joint genotyping.
  sample="${sample_by_role[$role]}"
  bam="${bam_by_role[$role]}"
  index="${index_by_role[$role]}"
  out_vcf="$OUT_ROOT/${sample}.g.vcf.gz"
  gvcfs+=( "$out_vcf" )
  # Emit diploid gVCFs across the full capture panel with one native PairHMM thread.
  # Disable downsampling of reads with the same alignment start.
  run gatk --java-options "$GATK_JAVA_OPTIONS" HaplotypeCaller \
    -R "$REF_FASTA" \
    -I "$bam" \
    --read-index "$index" \
    -ERC GVCF \
    --sample-ploidy 2 \
    --max-reads-per-alignment-start 0 \
    --native-pair-hmm-threads 1 \
    --intervals "$INTERVALS" \
    --tmp-dir "$TMP_DIR" \
    -O "$out_vcf"
done

# Combine source gVCFs, retaining separate sample records.
run gatk --java-options "$GATK_JAVA_OPTIONS" CombineGVCFs \
  -R "$REF_FASTA" \
  -V "${gvcfs[0]}" -V "${gvcfs[1]}" \
  --intervals "$INTERVALS" \
  --tmp-dir "$TMP_DIR" \
  -O "$OUT_ROOT/combined.g.vcf.gz"

# Genotype donor and WT at the same sites; marker selection follows later.
run gatk --java-options "$GATK_JAVA_OPTIONS" GenotypeGVCFs \
  -R "$REF_FASTA" \
  -V "$OUT_ROOT/combined.g.vcf.gz" \
  --sample-ploidy 2 \
  --intervals "$INTERVALS" \
  --tmp-dir "$TMP_DIR" \
  -O "$OUT_ROOT/joint.vcf.gz"

# Check generated files only after an actual run.
if [[ "$DRY_RUN" == "0" ]]; then
  # Set expected joint sample names, independently of column order.
  joint_samples="$(printf '%s\n' "${sample_by_role[pure_donor]}" "${sample_by_role[pure_wt]}" | LC_ALL=C sort)"
  # Check four VCFs and their indexes; a joint VCF may contain no variants.
  for vcf in "${gvcfs[@]}" "$OUT_ROOT/combined.g.vcf.gz" "$OUT_ROOT/joint.vcf.gz"; do
    [[ -s "$vcf" && -s "$vcf.tbi" ]] || die "Missing VCF or index: $vcf"
    bcftools view -Ov "$vcf" >/dev/null || die "Unreadable VCF: $vcf"
    bcftools view -H -r "$first_contig:1-1" "$vcf" >/dev/null || die "Unreadable VCF index: $vcf.tbi"
    # Require both sources in joint files and one in each source gVCF.
    expected="$joint_samples"
    if [[ "$vcf" == "${gvcfs[0]}" ]]; then expected="${sample_by_role[pure_donor]}"; fi
    if [[ "$vcf" == "${gvcfs[1]}" ]]; then expected="${sample_by_role[pure_wt]}"; fi
    # Compare sample names regardless of column order.
    actual="$(bcftools query -l "$vcf" | LC_ALL=C sort)"
    [[ "$actual" == "$expected" ]] || die "Unexpected VCF samples: $vcf"
  done
  # Report completed source genotyping.
  echo "DONE: source genotyping only; marker selection has not run"
else
  # Report a dry run without implying that calling ran.
  echo "DRY RUN COMPLETE: no outputs created and no variant calling performed"
fi
