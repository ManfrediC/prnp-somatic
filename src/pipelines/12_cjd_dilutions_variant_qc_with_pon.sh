#!/usr/bin/env bash
# ============================================================
# CJD+dilutions variant extraction and QC (with PoN):
#   Stage 1: Select PASS variants from annotated VCFs
#   Stage 2: Export per-sample tables with GATK VariantsToTable
#   Stage 3: Integrate readcount metrics and apply QC in R
# ============================================================
set -euo pipefail
shopt -s nullglob

# Keep caller-provided values so ENV_FILE does not silently override them.
CLI_DRY_RUN="${DRY_RUN-}"
CLI_JAVA_MEM_GB="${JAVA_MEM_GB-}"
CLI_ENABLE_AAF_FILTER="${ENABLE_AAF_FILTER-}"
CLI_AAF_THRESHOLD="${AAF_THRESHOLD-}"
CLI_WITH_PON_VARIANT_QC_GROUPS="${WITH_PON_VARIANT_QC_GROUPS-}"

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
ENABLE_AAF_FILTER="${CLI_ENABLE_AAF_FILTER:-${ENABLE_AAF_FILTER:-1}}"
AAF_THRESHOLD="${CLI_AAF_THRESHOLD:-${AAF_THRESHOLD:-0.0081}}"

# Threshold defaults mirror the methods section and can be overridden per run.
MIN_ALT_COUNT="${MIN_ALT_COUNT:-10}"
MIN_DP="${MIN_DP:-100}"
MIN_STRAND_ALT="${MIN_STRAND_ALT:-3}"
MIN_MEAN_BQ="${MIN_MEAN_BQ:-20}"
MIN_MEAN_MQ="${MIN_MEAN_MQ:-20}"
MAX_POP_FREQ="${MAX_POP_FREQ:-0.001}"
MAX_BINOM_P="${MAX_BINOM_P:-1e-6}"

# -----------------------
# Path defaults (repo-relative)
# -----------------------
MUTECT2_WITH_PON_OUT_ROOT="${MUTECT2_WITH_PON_OUT_ROOT:-runs/mutect2_cjd_dilutions_with_pon}"
if [[ "$MUTECT2_WITH_PON_OUT_ROOT" == run/* ]]; then
  MUTECT2_WITH_PON_OUT_ROOT="runs/${MUTECT2_WITH_PON_OUT_ROOT#run/}"
fi

MUTECT2_WITH_PON_VARIANT_QC_ROOT="${MUTECT2_WITH_PON_VARIANT_QC_ROOT:-$MUTECT2_WITH_PON_OUT_ROOT}"
if [[ "$MUTECT2_WITH_PON_VARIANT_QC_ROOT" == run/* ]]; then
  MUTECT2_WITH_PON_VARIANT_QC_ROOT="runs/${MUTECT2_WITH_PON_VARIANT_QC_ROOT#run/}"
fi

WITH_PON_VARIANT_QC_GROUPS="${CLI_WITH_PON_VARIANT_QC_GROUPS:-${WITH_PON_VARIANT_QC_GROUPS:-${WITH_PON_GROUPS:-cjd dilutions}}}"
WITH_PON_VARIANT_QC_VCF_SUBDIR="${WITH_PON_VARIANT_QC_VCF_SUBDIR:-annot_with_gnomad}"
WITH_PON_VARIANT_QC_METRICS_SUBDIR="${WITH_PON_VARIANT_QC_METRICS_SUBDIR:-readcount_qc/metrics}"

WITH_PON_VARIANT_QC_RESULTS_ROOT="${WITH_PON_VARIANT_QC_RESULTS_ROOT:-results/mutect2_cjd_dilutions_with_pon/variant_qc}"
WITH_PON_VARIANT_QC_R_SCRIPT="${WITH_PON_VARIANT_QC_R_SCRIPT:-src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R}"
MANUAL_POP_FREQ_TSV="${MANUAL_POP_FREQ_TSV:-resources/annotations/manual_population_freq.tsv}"


to_abs() {
  local p="$1"
  # Allow both repo-relative and absolute paths in ENV overrides.
  if [[ "$p" = /* ]]; then
    echo "$p"
  else
    echo "$REPO_ROOT/$p"
  fi
}

VARIANT_QC_ROOT="$(to_abs "$MUTECT2_WITH_PON_VARIANT_QC_ROOT")"
WITH_PON_VARIANT_QC_RESULTS_ROOT="$(to_abs "$WITH_PON_VARIANT_QC_RESULTS_ROOT")"
WITH_PON_VARIANT_QC_R_SCRIPT="$(to_abs "$WITH_PON_VARIANT_QC_R_SCRIPT")"
MANUAL_POP_FREQ_TSV="$(to_abs "$MANUAL_POP_FREQ_TSV")"

mkdir -p "$WITH_PON_VARIANT_QC_RESULTS_ROOT"

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

sample_from_vcf() {
  local path base sample
  path="$1"
  base="$(basename "$path")"
  sample="${base%.func.af.vcf.gz}"
  if [[ "$sample" == "$base" ]]; then
    sample="${base%.vcf.gz}"
  fi
  echo "$sample"
}

# -----------------------
# Tool + file checks
# -----------------------
have gatk    || die "gatk not in PATH (did you activate conda?)"
have tabix   || die "tabix not in PATH (did you activate conda?)"
have Rscript || die "Rscript not in PATH"

[[ -d "$VARIANT_QC_ROOT" ]] || die "Missing variant-QC root directory: $VARIANT_QC_ROOT"
[[ -s "$WITH_PON_VARIANT_QC_R_SCRIPT" ]] || die "Missing R QC script: $WITH_PON_VARIANT_QC_R_SCRIPT"
[[ -s "$MANUAL_POP_FREQ_TSV" ]] || die "Missing manual population frequency file: $MANUAL_POP_FREQ_TSV"

read -r -a groups <<< "$WITH_PON_VARIANT_QC_GROUPS"
[[ "${#groups[@]}" -gt 0 ]] || die "WITH_PON_VARIANT_QC_GROUPS is empty"

echo "== CJD+dilutions variant extraction + QC (with PoN) =="
echo "Repo root:                         $REPO_ROOT"
echo "VARIANT_QC_ROOT:                   $VARIANT_QC_ROOT"
echo "WITH_PON_GROUPS:                   $WITH_PON_VARIANT_QC_GROUPS"
echo "WITH_PON_VARIANT_QC_RESULTS_ROOT:  $WITH_PON_VARIANT_QC_RESULTS_ROOT"
echo "WITH_PON_VARIANT_QC_R_SCRIPT:      $WITH_PON_VARIANT_QC_R_SCRIPT"
echo "MANUAL_POP_FREQ_TSV:               $MANUAL_POP_FREQ_TSV"
echo "JAVA_MEM_GB:                       $JAVA_MEM_GB"
echo "ENABLE_AAF_FILTER:                 $ENABLE_AAF_FILTER"
echo "AAF_THRESHOLD:                     $AAF_THRESHOLD"
echo "MIN_ALT_COUNT:                     $MIN_ALT_COUNT"
echo "MIN_DP:                            $MIN_DP"
echo "MIN_STRAND_ALT:                    $MIN_STRAND_ALT"
echo "MIN_MEAN_BQ:                       $MIN_MEAN_BQ"
echo "MIN_MEAN_MQ:                       $MIN_MEAN_MQ"
echo "MAX_POP_FREQ:                      $MAX_POP_FREQ"
echo "MAX_BINOM_P:                       $MAX_BINOM_P"
echo

# Common fields
INFO_FIELDS=(
  CHROM POS REF ALT FILTER QUAL MQ QD AF GNOMAD_AF FS
  BaseQualityRankSumTest MappingQualityRankSumTest AS_SB_TABLE STRANDQ FUNCOTATION
)
FORMAT_FIELDS=(GT DP AD F1R2 F2R1 SB)

info_args=()
for f in "${INFO_FIELDS[@]}"; do
  info_args+=( -F "$f" )
done

format_args=()
for f in "${FORMAT_FIELDS[@]}"; do
  format_args+=( -GF "$f" )
done

for group in "${groups[@]}"; do
  echo "=== Group: $group ==="

  GROUP_ROOT="$VARIANT_QC_ROOT/$group"
  VCF_DIR="$GROUP_ROOT/$WITH_PON_VARIANT_QC_VCF_SUBDIR"
  METRICS_DIR="$GROUP_ROOT/$WITH_PON_VARIANT_QC_METRICS_SUBDIR"
  SELECT_DIR="$GROUP_ROOT/variant_qc/select_variants"
  TABLE_DIR="$GROUP_ROOT/variant_qc/variant_tables"
  OUTPUT_DIR="$WITH_PON_VARIANT_QC_RESULTS_ROOT/$group"

  mkdir -p "$SELECT_DIR" "$TABLE_DIR" "$OUTPUT_DIR"

  [[ -d "$VCF_DIR" ]] || die "Missing input VCF dir for group '$group': $VCF_DIR"
  [[ -d "$METRICS_DIR" ]] || die "Missing metrics dir for group '$group': $METRICS_DIR"

  vcfs=( "$VCF_DIR"/*.func.af.vcf.gz )
  if [[ "${#vcfs[@]}" -eq 0 ]]; then
    vcfs=( "$VCF_DIR"/*.vcf.gz )
  fi
  [[ "${#vcfs[@]}" -gt 0 ]] || die "No input VCFs found for group '$group' in: $VCF_DIR"

  samples=()
  for vcf in "${vcfs[@]}"; do
    [[ -f "$vcf" ]] || continue
    samples+=( "$(sample_from_vcf "$vcf")" )
  done
  [[ "${#samples[@]}" -gt 0 ]] || die "No valid sample names could be derived for group '$group'"

  # ------------------------------------------------------------
  # Stage 1: Select PASS variants
  # ------------------------------------------------------------
  echo "=== Stage 1 ($group): Select PASS variants ==="
  for vcf in "${vcfs[@]}"; do
    [[ -f "$vcf" ]] || continue
    sample="$(sample_from_vcf "$vcf")"
    out="$SELECT_DIR/${sample}.PASS.vcf.gz"
    [[ -s "$out" ]] && { echo "[Stage1][$group] SKIP $sample"; continue; }
    run gatk --java-options "-Xmx${JAVA_MEM_GB}g" SelectVariants \
      -V "$vcf" \
      --exclude-filtered \
      -O "$out"
    run tabix -f -p vcf "$out"
  done

  # ------------------------------------------------------------
  # Stage 2: Export tables with VariantsToTable
  # ------------------------------------------------------------
  echo "=== Stage 2 ($group): VariantsToTable export ==="
  for sample in "${samples[@]}"; do
    pass_vcf="$SELECT_DIR/${sample}.PASS.vcf.gz"
    out_tsv="$TABLE_DIR/${sample}.withPoN.tsv"

    if [[ ! -s "$pass_vcf" && "$DRY_RUN" != "1" ]]; then
      die "Missing PASS VCF for $sample: $pass_vcf"
    fi
    [[ -s "$out_tsv" ]] && { echo "[Stage2][$group] SKIP $sample"; continue; }

    run gatk --java-options "-Xmx${JAVA_MEM_GB}g" VariantsToTable \
      -V "$pass_vcf" \
      -O "$out_tsv" \
      "${info_args[@]}" \
      "${format_args[@]}"
  done

  if [[ "$DRY_RUN" == "0" ]]; then
    missing_tables=()
    for sample in "${samples[@]}"; do
      [[ -s "$TABLE_DIR/${sample}.withPoN.tsv" ]] || missing_tables+=( "$sample" )
    done
    [[ "${#missing_tables[@]}" -eq 0 ]] || die "Missing VariantsToTable output for group '$group': ${missing_tables[*]}"
  fi

  # ------------------------------------------------------------
  # Stage 3: R-based QC integration and filtering
  # ------------------------------------------------------------
  echo "=== Stage 3 ($group): R QC integration ==="
  # Dilution runs are used to derive the AAF threshold, so do not apply it here.
  group_enable_aaf_filter="$ENABLE_AAF_FILTER"
  if [[ "$group" == "dilutions" ]]; then
    group_enable_aaf_filter="0"
  fi
  echo "[Stage3][$group] enable_aaf_filter=$group_enable_aaf_filter (global default: $ENABLE_AAF_FILTER)"

  run Rscript "$WITH_PON_VARIANT_QC_R_SCRIPT" \
    --variant-dir "$TABLE_DIR" \
    --metrics-dir "$METRICS_DIR" \
    --manual-freq "$MANUAL_POP_FREQ_TSV" \
    --output-dir "$OUTPUT_DIR" \
    --enable-aaf-filter "$group_enable_aaf_filter" \
    --aaf-threshold "$AAF_THRESHOLD" \
    --min-alt-count "$MIN_ALT_COUNT" \
    --min-dp "$MIN_DP" \
    --min-strand-alt "$MIN_STRAND_ALT" \
    --min-mean-bq "$MIN_MEAN_BQ" \
    --min-mean-mq "$MIN_MEAN_MQ" \
    --max-pop-freq "$MAX_POP_FREQ" \
    --max-binom-p "$MAX_BINOM_P"

  if [[ "$DRY_RUN" == "0" ]]; then
    [[ -s "$OUTPUT_DIR/summary_combined_variants.tsv" ]] || die "Missing summary output: $OUTPUT_DIR/summary_combined_variants.tsv"
    [[ -s "$OUTPUT_DIR/filtered_variants.tsv" ]] || die "Missing filtered output: $OUTPUT_DIR/filtered_variants.tsv"
    [[ -s "$OUTPUT_DIR/filter_counts.tsv" ]] || die "Missing filter-count output: $OUTPUT_DIR/filter_counts.tsv"
  fi

  echo "=== Group finished: $group ==="
  echo
done

echo "=== CJD+dilutions variant extraction + QC (with PoN) finished ==="
