#!/usr/bin/env bash
# ============================================================
# Controls-only variant extraction and QC (no PoN):
#   Stage 1: Select PASS variants from annotated VCFs
#   Stage 2: Export per-sample tables with GATK VariantsToTable
#   Stage 3: Integrate readcount metrics and apply QC in R
# ============================================================
set -euo pipefail
shopt -s nullglob

# Keep caller-provided values so ENV_FILE does not silently override them.
# Precedence: explicit CLI env vars > values from ENV_FILE > script defaults.
CLI_DRY_RUN="${DRY_RUN-}"
CLI_JAVA_MEM_GB="${JAVA_MEM_GB-}"
CLI_ENABLE_AAF_FILTER="${ENABLE_AAF_FILTER-}"
CLI_AAF_THRESHOLD="${AAF_THRESHOLD-}"

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
  # ENV_FILE is optional: this script also runs with defaults only.
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
MUTECT2_CONTROLS_OUT_ROOT="${MUTECT2_CONTROLS_OUT_ROOT:-runs/mutect2_controls_no_pon}"
if [[ "$MUTECT2_CONTROLS_OUT_ROOT" == run/* ]]; then
  MUTECT2_CONTROLS_OUT_ROOT="runs/${MUTECT2_CONTROLS_OUT_ROOT#run/}"
fi

OUT_ROOT="${OUT_ROOT:-$MUTECT2_CONTROLS_OUT_ROOT}"

VARIANT_QC_ROOT="${VARIANT_QC_ROOT:-$OUT_ROOT/variant_qc}"
VARIANT_QC_VCF_DIR="${VARIANT_QC_VCF_DIR:-$OUT_ROOT/annot_with_gnomad}"
VARIANT_QC_READCOUNT_METRICS_DIR="${VARIANT_QC_READCOUNT_METRICS_DIR:-$OUT_ROOT/readcount_qc/metrics}"

VARIANT_QC_SELECT_DIR="${VARIANT_QC_SELECT_DIR:-$VARIANT_QC_ROOT/select_variants}"
VARIANT_QC_TABLE_DIR="${VARIANT_QC_TABLE_DIR:-$VARIANT_QC_ROOT/variant_tables}"
VARIANT_QC_OUTPUT_DIR="${VARIANT_QC_OUTPUT_DIR:-results/mutect2_controls_no_pon/variant_qc}"

VARIANT_QC_R_SCRIPT="${VARIANT_QC_R_SCRIPT:-src/pipelines/6_controls_variant_table_qc_no_pon.R}"
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

VARIANT_QC_VCF_DIR="$(to_abs "$VARIANT_QC_VCF_DIR")"
VARIANT_QC_READCOUNT_METRICS_DIR="$(to_abs "$VARIANT_QC_READCOUNT_METRICS_DIR")"
VARIANT_QC_SELECT_DIR="$(to_abs "$VARIANT_QC_SELECT_DIR")"
VARIANT_QC_TABLE_DIR="$(to_abs "$VARIANT_QC_TABLE_DIR")"
VARIANT_QC_OUTPUT_DIR="$(to_abs "$VARIANT_QC_OUTPUT_DIR")"
VARIANT_QC_R_SCRIPT="$(to_abs "$VARIANT_QC_R_SCRIPT")"
MANUAL_POP_FREQ_TSV="$(to_abs "$MANUAL_POP_FREQ_TSV")"

mkdir -p "$VARIANT_QC_SELECT_DIR" "$VARIANT_QC_TABLE_DIR" "$VARIANT_QC_OUTPUT_DIR"

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
  # Return just the sample stem so all stage outputs share one naming key.
  echo "$sample"
}

# -----------------------
# Tool + file checks
# -----------------------
have gatk    || die "gatk not in PATH (did you activate conda?)"
have tabix   || die "tabix not in PATH (did you activate conda?)"
have Rscript || die "Rscript not in PATH"

[[ -d "$VARIANT_QC_VCF_DIR" ]] || die "Missing VARIANT_QC_VCF_DIR: $VARIANT_QC_VCF_DIR"
[[ -d "$VARIANT_QC_READCOUNT_METRICS_DIR" ]] || die "Missing VARIANT_QC_READCOUNT_METRICS_DIR: $VARIANT_QC_READCOUNT_METRICS_DIR"
[[ -s "$VARIANT_QC_R_SCRIPT" ]] || die "Missing R QC script: $VARIANT_QC_R_SCRIPT"
[[ -s "$MANUAL_POP_FREQ_TSV" ]] || die "Missing manual population frequency file: $MANUAL_POP_FREQ_TSV"

vcfs=( "$VARIANT_QC_VCF_DIR"/*.func.af.vcf.gz )
if [[ "${#vcfs[@]}" -eq 0 ]]; then
  vcfs=( "$VARIANT_QC_VCF_DIR"/*.vcf.gz )
fi
[[ "${#vcfs[@]}" -gt 0 ]] || die "No input VCFs found in: $VARIANT_QC_VCF_DIR"

samples=()
for vcf in "${vcfs[@]}"; do
  [[ -f "$vcf" ]] || continue
  # Stage 2/3 iterate over this list; deriving once keeps sample ordering stable.
  samples+=( "$(sample_from_vcf "$vcf")" )
done
[[ "${#samples[@]}" -gt 0 ]] || die "No valid sample names could be derived from VCF filenames"

echo "== Controls variant extraction + QC (no PoN) =="
echo "Repo root:                        $REPO_ROOT"
echo "VARIANT_QC_VCF_DIR:               $VARIANT_QC_VCF_DIR"
echo "VARIANT_QC_READCOUNT_METRICS_DIR: $VARIANT_QC_READCOUNT_METRICS_DIR"
echo "VARIANT_QC_SELECT_DIR:            $VARIANT_QC_SELECT_DIR"
echo "VARIANT_QC_TABLE_DIR:             $VARIANT_QC_TABLE_DIR"
echo "VARIANT_QC_OUTPUT_DIR:            $VARIANT_QC_OUTPUT_DIR"
echo "VARIANT_QC_R_SCRIPT:              $VARIANT_QC_R_SCRIPT"
echo "MANUAL_POP_FREQ_TSV:              $MANUAL_POP_FREQ_TSV"
echo "JAVA_MEM_GB:                      $JAVA_MEM_GB"
echo "ENABLE_AAF_FILTER:                $ENABLE_AAF_FILTER"
echo "AAF_THRESHOLD:                    $AAF_THRESHOLD"
echo "MIN_ALT_COUNT:                    $MIN_ALT_COUNT"
echo "MIN_DP:                           $MIN_DP"
echo "MIN_STRAND_ALT:                   $MIN_STRAND_ALT"
echo "MIN_MEAN_BQ:                      $MIN_MEAN_BQ"
echo "MIN_MEAN_MQ:                      $MIN_MEAN_MQ"
echo "MAX_POP_FREQ:                     $MAX_POP_FREQ"
echo "MAX_BINOM_P:                      $MAX_BINOM_P"
echo
# Configuration echo above acts as a run manifest for troubleshooting and provenance.

# ------------------------------------------------------------
# Stage 1: Select PASS variants
# ------------------------------------------------------------
echo "=== Stage 1: Select PASS variants ==="
# Even if upstream files are expected to be PASS-only, this keeps the stage idempotent.
for vcf in "${vcfs[@]}"; do
  [[ -f "$vcf" ]] || continue
  sample="$(sample_from_vcf "$vcf")"
  out="$VARIANT_QC_SELECT_DIR/${sample}.PASS.vcf.gz"
  [[ -s "$out" ]] && { echo "[Stage1] SKIP $sample"; continue; }
  run gatk --java-options "-Xmx${JAVA_MEM_GB}g" SelectVariants \
    -V "$vcf" \
    --exclude-filtered \
    -O "$out"
  run tabix -f -p vcf "$out"
done

# ------------------------------------------------------------
# Stage 2: Export tables with VariantsToTable
# ------------------------------------------------------------
echo "=== Stage 2: VariantsToTable export ==="
INFO_FIELDS=(
  CHROM POS REF ALT FILTER QUAL MQ QD AF GNOMAD_AF FS
  BaseQualityRankSumTest MappingQualityRankSumTest AS_SB_TABLE STRANDQ FUNCOTATION
)
FORMAT_FIELDS=(GT DP AD F1R2 F2R1 SB)
# INFO fields: site-level annotations.
# FORMAT fields: per-sample depth/strand evidence needed for QC rules.

info_args=()
for f in "${INFO_FIELDS[@]}"; do
  info_args+=( -F "$f" )
done

format_args=()
for f in "${FORMAT_FIELDS[@]}"; do
  format_args+=( -GF "$f" )
done

for sample in "${samples[@]}"; do
  # Keep per-sample table filenames aligned with legacy *.noPoN.tsv outputs.
  pass_vcf="$VARIANT_QC_SELECT_DIR/${sample}.PASS.vcf.gz"
  out_tsv="$VARIANT_QC_TABLE_DIR/${sample}.noPoN.tsv"

  if [[ ! -s "$pass_vcf" && "$DRY_RUN" != "1" ]]; then
    die "Missing PASS VCF for $sample: $pass_vcf"
  fi
  [[ -s "$out_tsv" ]] && { echo "[Stage2] SKIP $sample"; continue; }

  run gatk --java-options "-Xmx${JAVA_MEM_GB}g" VariantsToTable \
    -V "$pass_vcf" \
    -O "$out_tsv" \
    "${info_args[@]}" \
    "${format_args[@]}"
done

if [[ "$DRY_RUN" == "0" ]]; then
  # Fail early if any sample table is missing before invoking R integration.
  missing_tables=()
  for sample in "${samples[@]}"; do
    [[ -s "$VARIANT_QC_TABLE_DIR/${sample}.noPoN.tsv" ]] || missing_tables+=( "$sample" )
  done
  [[ "${#missing_tables[@]}" -eq 0 ]] || die "Missing VariantsToTable output for: ${missing_tables[*]}"
fi

# ------------------------------------------------------------
# Stage 3: R-based QC integration and filtering
# ------------------------------------------------------------
echo "=== Stage 3: R QC integration ==="
# Keep QC logic centralised in the R script so criteria changes happen in one place.
# This stage merges table + readcount metrics and emits final filtered summaries.
run Rscript "$VARIANT_QC_R_SCRIPT" \
  --variant-dir "$VARIANT_QC_TABLE_DIR" \
  --metrics-dir "$VARIANT_QC_READCOUNT_METRICS_DIR" \
  --manual-freq "$MANUAL_POP_FREQ_TSV" \
  --output-dir "$VARIANT_QC_OUTPUT_DIR" \
  --enable-aaf-filter "$ENABLE_AAF_FILTER" \
  --aaf-threshold "$AAF_THRESHOLD" \
  --min-alt-count "$MIN_ALT_COUNT" \
  --min-dp "$MIN_DP" \
  --min-strand-alt "$MIN_STRAND_ALT" \
  --min-mean-bq "$MIN_MEAN_BQ" \
  --min-mean-mq "$MIN_MEAN_MQ" \
  --max-pop-freq "$MAX_POP_FREQ" \
  --max-binom-p "$MAX_BINOM_P"

if [[ "$DRY_RUN" == "0" ]]; then
  # Final contract: these three artefacts are required by downstream review steps.
  [[ -s "$VARIANT_QC_OUTPUT_DIR/summary_combined_variants.tsv" ]] || die "Missing summary output: $VARIANT_QC_OUTPUT_DIR/summary_combined_variants.tsv"
  [[ -s "$VARIANT_QC_OUTPUT_DIR/filtered_variants.tsv" ]] || die "Missing filtered output: $VARIANT_QC_OUTPUT_DIR/filtered_variants.tsv"
  [[ -s "$VARIANT_QC_OUTPUT_DIR/filter_counts.tsv" ]] || die "Missing filter-count output: $VARIANT_QC_OUTPUT_DIR/filter_counts.tsv"
fi

echo "=== Variant extraction + QC pipeline finished ==="
