#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Repository setup and local configuration
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$REPO_ROOT"

# Load local overrides when present so long reruns can be configured without
# editing the script itself.
CONFIG_FILE="${CONFIG_FILE:-${REPO_ROOT}/config/repeats.env}"
if [[ -f "$CONFIG_FILE" ]]; then
  set -a
  # shellcheck source=/dev/null
  source "$CONFIG_FILE"
  set +a
fi

# ---------------------------------------------------------------------------
# User-facing workflow settings
# ---------------------------------------------------------------------------

REPEAT_BAM_DIR="${REPEAT_BAM_DIR:-results/final_bam}"
REPEAT_RESULTS_ROOT="${REPEAT_RESULTS_ROOT:-results/repeats}"
REPEAT_REF_FASTA="${REPEAT_REF_FASTA:-resources/references/snv/chr2_chr4_chr20.fasta}"
PRNP_ORR_BED="${PRNP_ORR_BED:-resources/repeats/prnp_orr.hg38.bed}"
PRNP_EH_CATALOG="${PRNP_EH_CATALOG:-resources/repeats/prnp_orr.expansionhunter.json}"
RUN_GANGSTR="${RUN_GANGSTR:-0}"
GANGSTR_BIN="${GANGSTR_BIN:-}"
GANGSTR_REGIONS_BED="${GANGSTR_REGIONS_BED:-resources/repeats/prnp_orr.gangstr.bed}"
# These may be set explicitly, but empty defaults mean "estimate from each BAM"
# so GangSTR stays usable across cohorts with different library properties.
GANGSTR_READLENGTH="${GANGSTR_READLENGTH:-}"
GANGSTR_INSERTMEAN="${GANGSTR_INSERTMEAN:-}"
GANGSTR_INSERTSDEV="${GANGSTR_INSERTSDEV:-}"
GANGSTR_MAX_PROC_READ="${GANGSTR_MAX_PROC_READ:-1000000}"
PRNP_TOTAL_REFERENCE_REPEATS="${PRNP_TOTAL_REFERENCE_REPEATS:-5}"
PRNP_VARIABLE_REPEAT_OFFSET="${PRNP_VARIABLE_REPEAT_OFFSET:-3}"
REPEAT_THREADS="${REPEAT_THREADS:-4}"
ORR_MIN_READS="${ORR_MIN_READS:-10}"
FORCE="${FORCE:-0}"
ARCHIVE_EXISTING_RUN="${ARCHIVE_EXISTING_RUN:-0}"
ARCHIVE_RUN_LABEL="${ARCHIVE_RUN_LABEL:-}"
EXPECTED_SAMPLE_COUNT="${EXPECTED_SAMPLE_COUNT:-32}"

normalize_path_setting() {
  # Support Windows-style absolute paths when the workflow is launched from WSL.
  local value="$1"
  local normalized

  if [[ -z "$value" ]]; then
    echo ""
    return 0
  fi

  # Quick pass for common Windows drive-letter paths like C:\path\to\file.
  if [[ "$value" == [A-Za-z]:\\* ]]; then
    if command -v wslpath >/dev/null 2>&1; then
      normalized="$(wslpath -u "$value")"
      echo "$normalized"
      return 0
    fi
  fi

  # Replace backslashes for any remaining Windows-style path separators.
  echo "${value//\\//}"
}

# Normalise configured paths before constructing any derived output paths.
REPEAT_BAM_DIR="$(normalize_path_setting "$REPEAT_BAM_DIR")"
REPEAT_RESULTS_ROOT="$(normalize_path_setting "$REPEAT_RESULTS_ROOT")"
REPEAT_REF_FASTA="$(normalize_path_setting "$REPEAT_REF_FASTA")"
PRNP_ORR_BED="$(normalize_path_setting "$PRNP_ORR_BED")"
PRNP_EH_CATALOG="$(normalize_path_setting "$PRNP_EH_CATALOG")"
GANGSTR_BIN="$(normalize_path_setting "$GANGSTR_BIN")"
GANGSTR_REGIONS_BED="$(normalize_path_setting "$GANGSTR_REGIONS_BED")"

# ---------------------------------------------------------------------------
# Derived live-output paths
# ---------------------------------------------------------------------------

RAW_ROOT="${REPEAT_RESULTS_ROOT}/raw/expansionhunter"
GANGSTR_RAW_ROOT="${REPEAT_RESULTS_ROOT}/raw/gangstr"
LOG_ROOT="${REPEAT_RESULTS_ROOT}/logs"
MANIFEST_TSV="${REPEAT_RESULTS_ROOT}/sample_manifest.tsv"
RUN_SETTINGS_TSV="${REPEAT_RESULTS_ROOT}/run_settings.tsv"
OLD_RUNS_ROOT="${REPEAT_RESULTS_ROOT}/old_runs"
SAMPLE_CALLS_TSV="${REPEAT_RESULTS_ROOT}/sample_calls.tsv"
CANDIDATE_CALLS_TSV="${REPEAT_RESULTS_ROOT}/candidate_calls.tsv"
COHORT_SUMMARY_TSV="${REPEAT_RESULTS_ROOT}/cohort_summary.tsv"
SUBCLONAL_SUPPORT_TSV="${REPEAT_RESULTS_ROOT}/subclonal_read_support.tsv"
GANGSTR_CALLS_TSV="${REPEAT_RESULTS_ROOT}/gangstr_calls.tsv"
SOMATIC_SCREEN_TSV="${REPEAT_RESULTS_ROOT}/somatic_screen.tsv"
SUMMARY_PY="${SCRIPT_DIR}/03_summarize_prnp_orr.py"
SUBCLONAL_PY="${SCRIPT_DIR}/02_inspect_prnp_orr_subclonal.py"
GANGSTR_SUMMARY_PY="${SCRIPT_DIR}/04_summarize_gangstr.py"
SOMATIC_SCREEN_PY="${SCRIPT_DIR}/05_summarize_somatic_screen.py"
MANIFEST_TMP=""
RUN_TARGETS=(
  # These are the files/directories that define the current live repeat run.
  # If any exist, they are either archived first or the run stops.
  "$MANIFEST_TSV"
  "$RUN_SETTINGS_TSV"
  "$SAMPLE_CALLS_TSV"
  "$CANDIDATE_CALLS_TSV"
  "$COHORT_SUMMARY_TSV"
  "$SUBCLONAL_SUPPORT_TSV"
  "$GANGSTR_CALLS_TSV"
  "$SOMATIC_SCREEN_TSV"
  "$RAW_ROOT"
  "$GANGSTR_RAW_ROOT"
  "$LOG_ROOT"
)

if [[ "${1:-}" == "--preflight" ]]; then
  printf 'repeat_bam_dir\t%s\n' "$REPEAT_BAM_DIR"
  printf 'repeat_results_root\t%s\n' "$REPEAT_RESULTS_ROOT"
  printf 'repeat_ref_fasta\t%s\n' "$REPEAT_REF_FASTA"
  printf 'prnp_orr_bed\t%s\n' "$PRNP_ORR_BED"
  printf 'prnp_eh_catalog\t%s\n' "$PRNP_EH_CATALOG"
  printf 'gangstr_bin\t%s\n' "$GANGSTR_BIN"
  printf 'gangstr_regions_bed\t%s\n' "$GANGSTR_REGIONS_BED"
  exit 0
fi

# ---------------------------------------------------------------------------
# Basic helper functions
# ---------------------------------------------------------------------------

require_cmd() {
  # Stop immediately when a core dependency is missing so long runs do not fail
  # only after partial work has already been written.
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Required command not found: $1" >&2
    exit 1
  fi
}

require_cmd_with_hint() {
  local cmd_name="$1"
  local hint="$2"
  # Provide an actionable installation/activation hint for optional tools that
  # may live outside the main repeat environment.
  if ! command -v "$cmd_name" >/dev/null 2>&1; then
    echo "Required command not found: $cmd_name" >&2
    echo "$hint" >&2
    exit 1
  fi
}

resolve_python_cmd() {
  # Prefer python3 in WSL and modern Linux environments, while still
  # supporting environments where only `python` exists.
  if command -v python3 >/dev/null 2>&1; then
    echo python3
    return 0
  fi
  if command -v python >/dev/null 2>&1; then
    echo python
    return 0
  fi
  return 1
}

config_file_loaded() {
  # Keep the provenance table explicit about whether local overrides were used.
  [[ -f "$CONFIG_FILE" ]]
}

capture_expansionhunter_version() {
  local version_output
  # Normalize the tool-reported version string into one stable provenance value.
  version_output="$(ExpansionHunter --version 2>&1 || true)"
  printf '%s\n' "$version_output" \
    | grep -Eo 'ExpansionHunter v[0-9]+(\.[0-9]+)+' \
    | head -n 1 \
    || true
}

capture_gangstr_version() {
  local version_output
  # GangSTR reports a simple first-line version string, which is enough for
  # reproducibility tracking in run_settings.tsv.
  version_output="$("${GANGSTR_CMD[@]}" --version 2>&1 || true)"
  printf '%s\n' "$version_output" | head -n 1 || true
}

# GangSTR wants explicit library-size parameters. When the user does not supply
# them, estimate read length from a small prefix of alignments so reruns remain
# reproducible without hard-coding values for every sequencing batch.
estimate_bam_read_length() {
  local bam_path="$1"
  local status
  set +o pipefail
  samtools view "$bam_path" \
    | awk '
        NR <= 2000 && length($10) > 0 { sum += length($10); count += 1 }
        NR >= 2000 { exit }
        END {
          if (count == 0) exit 1
          printf "%.0f\n", sum / count
        }
      '
  status=$?
  set -o pipefail
  return "$status"
}

# Likewise estimate insert-size mean/SD from a bounded sample of proper pairs.
# The intent is not to build a perfect library model, only to provide stable
# sample-specific defaults that are much closer than a generic hard-coded value.
estimate_bam_insert_stats() {
  local bam_path="$1"
  local status
  set +o pipefail
  samtools view -f 67 -F 3852 "$bam_path" \
    | awk '
        {
          template_length = $9
          if (template_length < 0) template_length = -template_length
          if (template_length > 0 && template_length < 2000) {
            count += 1
            sum += template_length
            sumsq += template_length * template_length
            if (count >= 50000) exit
          }
        }
        END {
          if (count < 100) exit 1
          mean = sum / count
          variance = (sumsq / count) - (mean * mean)
          if (variance < 1) variance = 1
          printf "%.6f\t%.6f\n", mean, sqrt(variance)
        }
      '
  status=$?
  set -o pipefail
  return "$status"
}

require_file() {
  # Centralize file-existence checks so missing references fail with a uniform
  # error message before any expensive analysis starts.
  if [[ ! -f "$1" ]]; then
    echo "Required file not found: $1" >&2
    exit 1
  fi
}

# Remove an unfinished manifest if the script exits early. This keeps reruns
# from inheriting a partially written cohort manifest.
cleanup_tmp_manifest() {
  if [[ -n "$MANIFEST_TMP" && -f "$MANIFEST_TMP" ]]; then
    rm -f "$MANIFEST_TMP"
  fi
}

# Archive the previous live run before starting a new one. The workflow keeps
# only one live run under results/repeats; everything older moves to old_runs.
archive_existing_run() {
  local existing_targets=()
  local target archive_id archive_dir

  # First detect whether there is anything in the live output tree that would
  # otherwise be mixed into the next rerun.
  for target in "${RUN_TARGETS[@]}"; do
    if [[ -e "$target" ]]; then
      existing_targets+=("$target")
    fi
  done

  if [[ "${#existing_targets[@]}" -eq 0 ]]; then
    return 0
  fi

  if [[ "$ARCHIVE_EXISTING_RUN" != "1" ]]; then
    echo "Existing repeat outputs detected under $REPEAT_RESULTS_ROOT." >&2
    echo "Set ARCHIVE_EXISTING_RUN=1 to move them into $OLD_RUNS_ROOT before starting a clean rerun." >&2
    exit 1
  fi

  mkdir -p "$OLD_RUNS_ROOT"
  # Sanitize any user-provided archive label so it becomes a safe directory
  # name; otherwise fall back to a timestamped archive ID.
  archive_id="${ARCHIVE_RUN_LABEL//[^A-Za-z0-9_.-]/_}"
  if [[ -z "$archive_id" ]]; then
    archive_id="run_$(date -u +%Y%m%dT%H%M%SZ)"
  fi
  archive_dir="${OLD_RUNS_ROOT}/${archive_id}"
  if [[ -e "$archive_dir" ]]; then
    echo "Archive destination already exists: $archive_dir" >&2
    exit 1
  fi
  mkdir -p "$archive_dir"

  for target in "${existing_targets[@]}"; do
    mv "$target" "$archive_dir/"
  done

  echo "Archived previous repeat outputs to $archive_dir"
}

# Verify that the final live run is internally complete before the script
# reports success. This catches interrupted runs that produced only a subset of
# sample outputs or skipped one of the summary tables.
verify_latest_run() {
  local manifest_path="$1"
  local expected_count="$2"
  local manifest_unique_count sample_calls_count
  local subclonal_count gangstr_count somatic_screen_count
  local sample_id group bam bai orr_reads sample_raw_dir
  local eh_prefix eh_vcf eh_json
  local gangstr_prefix gangstr_vcf
  local failures=()

  # The manifest is the canonical cohort definition for this run, so we check
  # its unique sample count first.
  manifest_unique_count="$(
    tail -n +2 "$manifest_path" | cut -f1 | LC_ALL=C sort -u | wc -l | awk '{print $1}'
  )"
  if [[ -n "$expected_count" && "$manifest_unique_count" -ne "$expected_count" ]]; then
    failures+=(
      "Manifest unique sample count ${manifest_unique_count} does not match expected count ${expected_count}."
    )
  fi

  # For each manifest sample, confirm that the caller outputs exist.
  while IFS=$'\t' read -r sample_id group bam bai orr_reads; do
    if [[ "$sample_id" == "sample_id" ]]; then
      continue
    fi

    sample_raw_dir="${RAW_ROOT}/${sample_id}"
    eh_prefix="${sample_raw_dir}/${sample_id}"
    eh_vcf="${eh_prefix}.vcf"
    eh_json="${eh_prefix}.json"
    gangstr_prefix="${GANGSTR_RAW_ROOT}/${sample_id}/${sample_id}"
    gangstr_vcf="${gangstr_prefix}.vcf"

    if [[ ! -f "$eh_vcf" ]]; then
      failures+=("Missing ExpansionHunter VCF for ${sample_id}: ${eh_vcf}")
    fi
    if [[ ! -f "$eh_json" ]]; then
      failures+=("Missing ExpansionHunter JSON for ${sample_id}: ${eh_json}")
    fi
    if [[ "$RUN_GANGSTR" == "1" && ! -f "$gangstr_vcf" ]]; then
      failures+=("Missing GangSTR VCF for ${sample_id}: ${gangstr_vcf}")
    fi
  done <"$manifest_path"

  # The top-level summary files are what downstream interpretation will consume,
  # so they must all exist in the live run directory.
  for target in \
    "$RUN_SETTINGS_TSV" \
    "$SAMPLE_CALLS_TSV" \
    "$CANDIDATE_CALLS_TSV" \
    "$COHORT_SUMMARY_TSV" \
    "$SUBCLONAL_SUPPORT_TSV" \
    "$GANGSTR_CALLS_TSV" \
    "$SOMATIC_SCREEN_TSV"; do
    if [[ ! -f "$target" ]]; then
      failures+=("Missing summary output: ${target}")
    fi
  done

  # Cross-check the two cohort-level per-sample tables against the manifest so
  # we do not silently accept truncated summarization.
  if [[ -f "$SAMPLE_CALLS_TSV" ]]; then
    sample_calls_count="$(
      tail -n +2 "$SAMPLE_CALLS_TSV" | cut -f1 | LC_ALL=C sort -u | wc -l | awk '{print $1}'
    )"
    if [[ "$sample_calls_count" -ne "$manifest_unique_count" ]]; then
      failures+=(
        "sample_calls.tsv contains ${sample_calls_count} unique samples but manifest expects ${manifest_unique_count}."
      )
    fi
  fi

  if [[ -f "$SUBCLONAL_SUPPORT_TSV" ]]; then
    subclonal_count="$(
      tail -n +2 "$SUBCLONAL_SUPPORT_TSV" | cut -f1 | LC_ALL=C sort -u | wc -l | awk '{print $1}'
    )"
    if [[ "$subclonal_count" -ne "$manifest_unique_count" ]]; then
      failures+=(
        "subclonal_read_support.tsv contains ${subclonal_count} unique samples but manifest expects ${manifest_unique_count}."
      )
    fi
  fi

  if [[ -f "$GANGSTR_CALLS_TSV" ]]; then
    gangstr_count="$(
      tail -n +2 "$GANGSTR_CALLS_TSV" | cut -f1 | LC_ALL=C sort -u | wc -l | awk '{print $1}'
    )"
    if [[ "$gangstr_count" -ne "$manifest_unique_count" ]]; then
      failures+=(
        "gangstr_calls.tsv contains ${gangstr_count} unique samples but manifest expects ${manifest_unique_count}."
      )
    fi
  fi

  if [[ -f "$SOMATIC_SCREEN_TSV" ]]; then
    somatic_screen_count="$(
      tail -n +2 "$SOMATIC_SCREEN_TSV" | cut -f1 | LC_ALL=C sort -u | wc -l | awk '{print $1}'
    )"
    if [[ "$somatic_screen_count" -ne "$manifest_unique_count" ]]; then
      failures+=(
        "somatic_screen.tsv contains ${somatic_screen_count} unique samples but manifest expects ${manifest_unique_count}."
      )
    fi
  fi

  if [[ "${#failures[@]}" -gt 0 ]]; then
    printf '%s\n' "${failures[@]}" >&2
    exit 1
  fi

  echo "Verified repeat outputs for ${manifest_unique_count} samples."
}

trap cleanup_tmp_manifest EXIT

# ---------------------------------------------------------------------------
# Dependency and input validation
# ---------------------------------------------------------------------------

require_cmd samtools
require_cmd_with_hint \
  ExpansionHunter \
  "Activate the repeat caller environment first: conda activate prnp-repeats"
PYTHON_CMD="$(resolve_python_cmd)"
if [[ -z "${PYTHON_CMD:-}" ]]; then
  echo "Required command not found: python or python3" >&2
  exit 1
fi

require_file "$REPEAT_REF_FASTA"
require_file "${REPEAT_REF_FASTA}.fai"
require_file "$PRNP_ORR_BED"
require_file "$PRNP_EH_CATALOG"
require_file "$SUMMARY_PY"
require_file "$SUBCLONAL_PY"
require_file "$GANGSTR_SUMMARY_PY"
require_file "$SOMATIC_SCREEN_PY"

if [[ ! -d "$REPEAT_BAM_DIR" ]]; then
  echo "Repeat BAM directory not found: $REPEAT_BAM_DIR" >&2
  exit 1
fi

GANGSTR_CMD=(GangSTR)
if [[ "$RUN_GANGSTR" == "1" ]]; then
  if [[ -n "$GANGSTR_BIN" ]]; then
    require_file "$GANGSTR_BIN"
    GANGSTR_CMD=("$GANGSTR_BIN")
  else
    require_cmd_with_hint \
      GangSTR \
      "Update the repeat environment first: conda env update -f env/repeats.environment.yml"
  fi
  require_file "$GANGSTR_REGIONS_BED"
fi

# ---------------------------------------------------------------------------
# Prepare a clean live output tree
# ---------------------------------------------------------------------------

mkdir -p "$REPEAT_RESULTS_ROOT"
archive_existing_run
mkdir -p "$RAW_ROOT" "$LOG_ROOT"
if [[ "$RUN_GANGSTR" == "1" ]]; then
  mkdir -p "$GANGSTR_RAW_ROOT"
fi

# ---------------------------------------------------------------------------
# Discover the cohort BAMs that belong to this repeat run
# ---------------------------------------------------------------------------

mapfile -t BAMS < <(
  find "$REPEAT_BAM_DIR" -maxdepth 1 -type f \
    \( -name 'CJD*.bam' -o -name 'Ctrl*.bam' \) \
    | LC_ALL=C sort -u
)

if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "No CJD/Control BAMs found in: $REPEAT_BAM_DIR" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Build a fresh sample manifest with per-sample coverage sanity checks
# ---------------------------------------------------------------------------

MANIFEST_TMP="$(mktemp "${REPEAT_RESULTS_ROOT}/sample_manifest.tsv.tmp.XXXXXX")"
echo -e "sample_id\tgroup\tbam_path\tbai_path\torr_overlapping_reads" >"$MANIFEST_TMP"
declare -A SEEN_SAMPLE_IDS=()
SAMPLE_COUNT=0
for BAM in "${BAMS[@]}"; do
  SAMPLE_ID="$(basename "$BAM" .bam)"
  if [[ "$SAMPLE_ID" == CJD* ]]; then
    GROUP="cjd"
  elif [[ "$SAMPLE_ID" == Ctrl* ]]; then
    GROUP="control"
  else
    continue
  fi

  # Fail fast on duplicate sample IDs rather than letting repeated rows cascade
  # into duplicated outputs and incorrect cohort summaries.
  if [[ -n "${SEEN_SAMPLE_IDS[$SAMPLE_ID]:-}" ]]; then
    echo "Duplicate sample ID encountered while building manifest: $SAMPLE_ID" >&2
    exit 1
  fi
  SEEN_SAMPLE_IDS["$SAMPLE_ID"]=1

  BAI="${BAM}.bai"
  require_file "$BAI"
  # Count reads overlapping the ORR interval before launching the expensive
  # caller so we can stop early on obviously unusable samples.
  ORR_READS="$(samtools view -c -L "$PRNP_ORR_BED" "$BAM")"
  if [[ "$ORR_READS" -lt "$ORR_MIN_READS" ]]; then
    echo "Sample $SAMPLE_ID has only $ORR_READS reads over PRNP ORR (< $ORR_MIN_READS)." >&2
    exit 1
  fi
  echo -e "${SAMPLE_ID}\t${GROUP}\t${BAM}\t${BAI}\t${ORR_READS}" >>"$MANIFEST_TMP"
  SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
done

if [[ -n "$EXPECTED_SAMPLE_COUNT" && "$SAMPLE_COUNT" -ne "$EXPECTED_SAMPLE_COUNT" ]]; then
  echo "Expected ${EXPECTED_SAMPLE_COUNT} repeat samples but found ${SAMPLE_COUNT} in ${REPEAT_BAM_DIR}." >&2
  exit 1
fi

mv "$MANIFEST_TMP" "$MANIFEST_TSV"
MANIFEST_TMP=""

# ---------------------------------------------------------------------------
# Record run provenance for the live results directory
# ---------------------------------------------------------------------------

{
  echo -e "key\tvalue"
  echo -e "created_utc\t$(date -u +%FT%TZ)"
  echo -e "config_file\t$CONFIG_FILE"
  echo -e "config_file_loaded\t$(if config_file_loaded; then echo yes; else echo no; fi)"
  echo -e "bam_dir\t$REPEAT_BAM_DIR"
  echo -e "results_root\t$REPEAT_RESULTS_ROOT"
  echo -e "reference_fasta\t$REPEAT_REF_FASTA"
  echo -e "prnp_orr_bed\t$PRNP_ORR_BED"
  echo -e "prnp_orr_bed_sha256\t$(sha256sum "$PRNP_ORR_BED" | awk '{print $1}')"
  echo -e "eh_catalog\t$PRNP_EH_CATALOG"
  echo -e "eh_catalog_sha256\t$(sha256sum "$PRNP_EH_CATALOG" | awk '{print $1}')"
  echo -e "repeat_threads\t$REPEAT_THREADS"
  echo -e "orr_min_reads\t$ORR_MIN_READS"
  echo -e "run_gangstr\t$RUN_GANGSTR"
  echo -e "gangstr_bin\t${GANGSTR_BIN:-GangSTR}"
  echo -e "gangstr_regions_bed\t$GANGSTR_REGIONS_BED"
  echo -e "gangstr_readlength\t${GANGSTR_READLENGTH:-auto_from_bam}"
  echo -e "gangstr_insertmean\t${GANGSTR_INSERTMEAN:-auto_from_bam}"
  echo -e "gangstr_insertsdev\t${GANGSTR_INSERTSDEV:-auto_from_bam}"
  echo -e "gangstr_max_proc_read\t$GANGSTR_MAX_PROC_READ"
  echo -e "force\t$FORCE"
  echo -e "archive_existing_run\t$ARCHIVE_EXISTING_RUN"
  echo -e "archive_run_label\t${ARCHIVE_RUN_LABEL:-not_set}"
  echo -e "old_runs_root\t$OLD_RUNS_ROOT"
  echo -e "expected_sample_count\t$EXPECTED_SAMPLE_COUNT"
  echo -e "prnp_total_reference_repeats\t$PRNP_TOTAL_REFERENCE_REPEATS"
  echo -e "prnp_variable_repeat_offset\t$PRNP_VARIABLE_REPEAT_OFFSET"
  echo -e "sample_count\t$SAMPLE_COUNT"
  echo -e "expansionhunter_version\t$(capture_expansionhunter_version || echo unknown)"
  if [[ "$RUN_GANGSTR" == "1" ]]; then
    echo -e "gangstr_version\t$(capture_gangstr_version)"
  fi
  echo -e "samtools_version\t$(samtools --version | head -n 1)"
  echo -e "python_version\t$(${PYTHON_CMD} --version 2>&1)"
} >"$RUN_SETTINGS_TSV"

# ---------------------------------------------------------------------------
# Run ExpansionHunter and optional GangSTR for each manifest sample
# ---------------------------------------------------------------------------

while IFS=$'\t' read -r SAMPLE_ID GROUP BAM BAI ORR_READS <&3; do
  if [[ "$SAMPLE_ID" == "sample_id" ]]; then
    continue
  fi

  SAMPLE_RAW_DIR="${RAW_ROOT}/${SAMPLE_ID}"
  mkdir -p "$SAMPLE_RAW_DIR"

  EH_PREFIX="${SAMPLE_RAW_DIR}/${SAMPLE_ID}"
  EH_VCF="${EH_PREFIX}.vcf"
  EH_JSON="${EH_PREFIX}.json"

  # Reuse existing sample-level caller outputs unless FORCE=1. This makes
  # intentional resumptions cheap while still protecting fresh clean reruns.
  if [[
    "$FORCE" != "1" &&
    -f "$EH_VCF" &&
    -f "$EH_JSON"
  ]]; then
    echo "Skipping ExpansionHunter for $SAMPLE_ID (outputs exist)."
  else
    echo "Running ExpansionHunter for $SAMPLE_ID"
    ExpansionHunter \
      --reads "$BAM" \
      --reference "$REPEAT_REF_FASTA" \
      --variant-catalog "$PRNP_EH_CATALOG" \
      --output-prefix "$EH_PREFIX" \
      --analysis-mode streaming \
      --threads "$REPEAT_THREADS" \
      </dev/null \
      >"${LOG_ROOT}/${SAMPLE_ID}.expansionhunter.stdout.log" \
      2>"${LOG_ROOT}/${SAMPLE_ID}.expansionhunter.stderr.log"
  fi

  if [[ "$RUN_GANGSTR" == "1" ]]; then
    SAMPLE_GANGSTR_DIR="${GANGSTR_RAW_ROOT}/${SAMPLE_ID}"
    mkdir -p "$SAMPLE_GANGSTR_DIR"
    GANGSTR_PREFIX="${SAMPLE_GANGSTR_DIR}/${SAMPLE_ID}"
    GANGSTR_VCF="${GANGSTR_PREFIX}.vcf"

    # Resolve the GangSTR library parameters once per sample. This lets the
    # shell workflow remain cohort-agnostic while still recording what was used
    # into run_settings.tsv for provenance.
    gangstr_readlength="${GANGSTR_READLENGTH:-$(estimate_bam_read_length "$BAM")}"
    if [[ -n "$GANGSTR_INSERTMEAN" && -n "$GANGSTR_INSERTSDEV" ]]; then
      gangstr_insertmean="$GANGSTR_INSERTMEAN"
      gangstr_insertsdev="$GANGSTR_INSERTSDEV"
    else
      read -r gangstr_insertmean gangstr_insertsdev < <(estimate_bam_insert_stats "$BAM")
    fi

    if [[ "$FORCE" != "1" && -f "$GANGSTR_VCF" ]]; then
      echo "Skipping GangSTR for $SAMPLE_ID (outputs exist)."
    else
      echo "Running GangSTR for $SAMPLE_ID"
      "${GANGSTR_CMD[@]}" \
        --bam "$BAM" \
        --ref "$REPEAT_REF_FASTA" \
        --regions "$GANGSTR_REGIONS_BED" \
        --out "$GANGSTR_PREFIX" \
        --targeted \
        --nonuniform \
        --readlength "$gangstr_readlength" \
        --insertmean "$gangstr_insertmean" \
        --insertsdev "$gangstr_insertsdev" \
        --max-proc-read "$GANGSTR_MAX_PROC_READ" \
        </dev/null \
        >"${LOG_ROOT}/${SAMPLE_ID}.gangstr.stdout.log" \
        2>"${LOG_ROOT}/${SAMPLE_ID}.gangstr.stderr.log"
    fi
  fi
done 3<"$MANIFEST_TSV"

# ---------------------------------------------------------------------------
# Build cohort summary tables and verify the final live run
# ---------------------------------------------------------------------------

"$PYTHON_CMD" "$SUMMARY_PY" \
  --manifest "$MANIFEST_TSV" \
  --results-root "$REPEAT_RESULTS_ROOT" \
  --reference-total-repeats "$PRNP_TOTAL_REFERENCE_REPEATS" \
  --variable-repeat-offset "$PRNP_VARIABLE_REPEAT_OFFSET"

"$PYTHON_CMD" "$SUBCLONAL_PY" \
  --sample-calls "$SAMPLE_CALLS_TSV" \
  --output "$SUBCLONAL_SUPPORT_TSV"

"$PYTHON_CMD" "$GANGSTR_SUMMARY_PY" \
  --manifest "$MANIFEST_TSV" \
  --results-root "$REPEAT_RESULTS_ROOT" \
  --reference-total-repeats "$PRNP_TOTAL_REFERENCE_REPEATS" \
  --variable-repeat-offset "$PRNP_VARIABLE_REPEAT_OFFSET" \
  --gangstr-enabled "$RUN_GANGSTR"

"$PYTHON_CMD" "$SOMATIC_SCREEN_PY" \
  --results-root "$REPEAT_RESULTS_ROOT"

verify_latest_run "$MANIFEST_TSV" "$SAMPLE_COUNT"

echo "Repeat workflow complete."
