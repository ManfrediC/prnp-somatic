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
REPEAT_REF_FASTA="${REPEAT_REF_FASTA:-resources/chr2_chr4_chr20.fasta}"
PRNP_ORR_BED="${PRNP_ORR_BED:-resources/prnp_orr.hg38.bed}"
PRNP_EH_CATALOG="${PRNP_EH_CATALOG:-resources/repeats/prnp_orr.expansionhunter.json}"
REVIEWER_CONDA_ENV="${REVIEWER_CONDA_ENV:-}"
REVIEWER_BIN="${REVIEWER_BIN:-}"
PRNP_TOTAL_REFERENCE_REPEATS="${PRNP_TOTAL_REFERENCE_REPEATS:-5}"
PRNP_VARIABLE_REPEAT_OFFSET="${PRNP_VARIABLE_REPEAT_OFFSET:-3}"
REPEAT_THREADS="${REPEAT_THREADS:-4}"
ORR_MIN_READS="${ORR_MIN_READS:-10}"
RUN_REVIEWER="${RUN_REVIEWER:-1}"
FORCE="${FORCE:-0}"
ARCHIVE_EXISTING_RUN="${ARCHIVE_EXISTING_RUN:-0}"
ARCHIVE_RUN_LABEL="${ARCHIVE_RUN_LABEL:-}"
EXPECTED_SAMPLE_COUNT="${EXPECTED_SAMPLE_COUNT:-32}"

# ---------------------------------------------------------------------------
# Derived live-output paths
# ---------------------------------------------------------------------------

RAW_ROOT="${REPEAT_RESULTS_ROOT}/raw/expansionhunter"
REVIEW_ROOT="${REPEAT_RESULTS_ROOT}/review/reviewer"
LOG_ROOT="${REPEAT_RESULTS_ROOT}/logs"
MANIFEST_TSV="${REPEAT_RESULTS_ROOT}/sample_manifest.tsv"
RUN_SETTINGS_TSV="${REPEAT_RESULTS_ROOT}/run_settings.tsv"
OLD_RUNS_ROOT="${REPEAT_RESULTS_ROOT}/old_runs"
SAMPLE_CALLS_TSV="${REPEAT_RESULTS_ROOT}/sample_calls.tsv"
SAMPLE_REVIEW_TSV="${REPEAT_RESULTS_ROOT}/sample_review.tsv"
CANDIDATE_CALLS_TSV="${REPEAT_RESULTS_ROOT}/candidate_calls.tsv"
COHORT_SUMMARY_TSV="${REPEAT_RESULTS_ROOT}/cohort_summary.tsv"
SUMMARY_PY="${SCRIPT_DIR}/summarize_prnp_orr.py"
MANIFEST_TMP=""
RUN_TARGETS=(
  # These are the files/directories that define the current live repeat run.
  # If any exist, they are either archived first or the run stops.
  "$MANIFEST_TSV"
  "$RUN_SETTINGS_TSV"
  "$SAMPLE_CALLS_TSV"
  "$SAMPLE_REVIEW_TSV"
  "$CANDIDATE_CALLS_TSV"
  "$COHORT_SUMMARY_TSV"
  "$RAW_ROOT"
  "$REVIEW_ROOT"
  "$LOG_ROOT"
)

# ---------------------------------------------------------------------------
# Basic helper functions
# ---------------------------------------------------------------------------

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Required command not found: $1" >&2
    exit 1
  fi
}

require_cmd_with_hint() {
  local cmd_name="$1"
  local hint="$2"
  if ! command -v "$cmd_name" >/dev/null 2>&1; then
    echo "Required command not found: $cmd_name" >&2
    echo "$hint" >&2
    exit 1
  fi
}

find_conda_env_prefix() {
  local env_name="$1"
  conda env list | awk -v env_name="$env_name" '
    $1 == env_name { print $NF; found=1; exit }
    END { if (!found) exit 1 }
  '
}

config_file_loaded() {
  [[ -f "$CONFIG_FILE" ]]
}

capture_expansionhunter_version() {
  local version_output
  version_output="$(ExpansionHunter --version 2>&1 || true)"
  printf '%s\n' "$version_output" \
    | grep -Eo 'ExpansionHunter v[0-9]+(\.[0-9]+)+' \
    | head -n 1 \
    || true
}

capture_reviewer_version() {
  local version_output package_version

  version_output="$("${REVIEWER_CMD[@]}" --version 2>&1 || true)"
  package_version="$(
    printf '%s\n' "$version_output" \
      | grep -Eo 'REViewer v[0-9]+(\.[0-9]+)+' \
      | head -n 1 \
      || true
  )"
  if [[ -n "$package_version" ]]; then
    printf '%s\n' "$package_version"
    return 0
  fi

  if [[ -n "$REVIEWER_CONDA_ENV" ]] && command -v conda >/dev/null 2>&1; then
    package_version="$(
      conda list -n "$REVIEWER_CONDA_ENV" reviewer 2>/dev/null \
        | awk '$1 == "reviewer" { print "REViewer v" $2; exit }'
    )"
    if [[ -n "$package_version" ]]; then
      printf '%s\n' "$package_version"
      return 0
    fi
  fi

  printf '%s\n' "unknown"
}

require_file() {
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
  local manifest_unique_count sample_calls_count sample_review_count
  local sample_id group bam bai orr_reads sample_raw_dir sample_review_dir
  local eh_prefix eh_vcf eh_json eh_realigned_bam review_prefix review_svg
  local review_metrics review_phasing
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

  # For each manifest sample, confirm that the caller outputs exist and, when
  # enabled, that REViewer produced the expected review artifacts.
  while IFS=$'\t' read -r sample_id group bam bai orr_reads; do
    if [[ "$sample_id" == "sample_id" ]]; then
      continue
    fi

    sample_raw_dir="${RAW_ROOT}/${sample_id}"
    sample_review_dir="${REVIEW_ROOT}/${sample_id}"
    eh_prefix="${sample_raw_dir}/${sample_id}"
    eh_vcf="${eh_prefix}.vcf"
    eh_json="${eh_prefix}.json"
    eh_realigned_bam="${eh_prefix}_realigned.bam"
    review_prefix="${sample_review_dir}/${sample_id}.PRNP_ORR"
    review_svg="${review_prefix}.PRNP_ORR.svg"
    review_metrics="${review_prefix}.metrics.tsv"
    review_phasing="${review_prefix}.phasing.tsv"

    if [[ ! -f "$eh_vcf" ]]; then
      failures+=("Missing ExpansionHunter VCF for ${sample_id}: ${eh_vcf}")
    fi
    if [[ ! -f "$eh_json" ]]; then
      failures+=("Missing ExpansionHunter JSON for ${sample_id}: ${eh_json}")
    fi
    if [[ "$RUN_REVIEWER" == "1" ]]; then
      if [[ ! -f "$eh_realigned_bam" ]]; then
        failures+=("Missing ExpansionHunter realigned BAM for ${sample_id}: ${eh_realigned_bam}")
      fi
      if [[ ! -f "$review_svg" ]]; then
        failures+=("Missing REViewer SVG for ${sample_id}: ${review_svg}")
      fi
      if [[ ! -f "$review_metrics" ]]; then
        failures+=("Missing REViewer metrics for ${sample_id}: ${review_metrics}")
      fi
      if [[ ! -f "$review_phasing" ]]; then
        failures+=("Missing REViewer phasing table for ${sample_id}: ${review_phasing}")
      fi
    fi
  done <"$manifest_path"

  # The top-level summary files are what downstream interpretation will consume,
  # so they must all exist in the live run directory.
  for target in \
    "$RUN_SETTINGS_TSV" \
    "$SAMPLE_CALLS_TSV" \
    "$SAMPLE_REVIEW_TSV" \
    "$CANDIDATE_CALLS_TSV" \
    "$COHORT_SUMMARY_TSV"; do
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

  if [[ -f "$SAMPLE_REVIEW_TSV" ]]; then
    sample_review_count="$(
      tail -n +2 "$SAMPLE_REVIEW_TSV" | cut -f1 | LC_ALL=C sort -u | wc -l | awk '{print $1}'
    )"
    if [[ "$sample_review_count" -ne "$manifest_unique_count" ]]; then
      failures+=(
        "sample_review.tsv contains ${sample_review_count} unique samples but manifest expects ${manifest_unique_count}."
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
require_cmd python
REVIEWER_CMD=(REViewer)
if [[ "$RUN_REVIEWER" == "1" ]]; then
  # REViewer may live either in the active environment, a dedicated path, or a
  # separate conda environment. Resolve that once up front before the long run.
  if [[ -n "$REVIEWER_BIN" ]]; then
    require_file "$REVIEWER_BIN"
    REVIEWER_CMD=("$REVIEWER_BIN")
  elif [[ -n "$REVIEWER_CONDA_ENV" ]]; then
    require_cmd conda
    REVIEWER_PREFIX="$(find_conda_env_prefix "$REVIEWER_CONDA_ENV")"
    REVIEWER_BIN="${REVIEWER_PREFIX}/bin/REViewer"
    require_file "$REVIEWER_BIN"
    REVIEWER_CMD=("$REVIEWER_BIN")
  elif command -v conda >/dev/null 2>&1 && find_conda_env_prefix "prnp-reviewer" >/dev/null 2>&1; then
    # When no local repeat config is present, prefer the conventional dedicated
    # reviewer environment if it exists on this machine.
    REVIEWER_CONDA_ENV="prnp-reviewer"
    REVIEWER_PREFIX="$(find_conda_env_prefix "$REVIEWER_CONDA_ENV")"
    REVIEWER_BIN="${REVIEWER_PREFIX}/bin/REViewer"
    require_file "$REVIEWER_BIN"
    REVIEWER_CMD=("$REVIEWER_BIN")
  else
    require_cmd_with_hint \
      REViewer \
      "Install or point to REViewer, or set REVIEWER_CONDA_ENV=prnp-reviewer if it lives in a separate environment."
  fi
fi

require_file "$REPEAT_REF_FASTA"
require_file "${REPEAT_REF_FASTA}.fai"
require_file "$PRNP_ORR_BED"
require_file "$PRNP_EH_CATALOG"
require_file "$SUMMARY_PY"

if [[ ! -d "$REPEAT_BAM_DIR" ]]; then
  echo "Repeat BAM directory not found: $REPEAT_BAM_DIR" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Prepare a clean live output tree
# ---------------------------------------------------------------------------

mkdir -p "$REPEAT_RESULTS_ROOT"
archive_existing_run
mkdir -p "$RAW_ROOT" "$REVIEW_ROOT" "$LOG_ROOT"

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
  echo -e "run_reviewer\t$RUN_REVIEWER"
  echo -e "reviewer_conda_env\t$REVIEWER_CONDA_ENV"
  echo -e "reviewer_bin\t$REVIEWER_BIN"
  echo -e "force\t$FORCE"
  echo -e "archive_existing_run\t$ARCHIVE_EXISTING_RUN"
  echo -e "archive_run_label\t${ARCHIVE_RUN_LABEL:-not_set}"
  echo -e "old_runs_root\t$OLD_RUNS_ROOT"
  echo -e "expected_sample_count\t$EXPECTED_SAMPLE_COUNT"
  echo -e "prnp_total_reference_repeats\t$PRNP_TOTAL_REFERENCE_REPEATS"
  echo -e "prnp_variable_repeat_offset\t$PRNP_VARIABLE_REPEAT_OFFSET"
  echo -e "sample_count\t$SAMPLE_COUNT"
  echo -e "expansionhunter_version\t$(capture_expansionhunter_version || echo unknown)"
  if [[ "$RUN_REVIEWER" == "1" ]]; then
    echo -e "reviewer_version\t$(capture_reviewer_version)"
  fi
  echo -e "samtools_version\t$(samtools --version | head -n 1)"
  echo -e "python_version\t$(python --version 2>&1)"
} >"$RUN_SETTINGS_TSV"

# ---------------------------------------------------------------------------
# Run ExpansionHunter and optional REViewer for each manifest sample
# ---------------------------------------------------------------------------

while IFS=$'\t' read -r SAMPLE_ID GROUP BAM BAI ORR_READS <&3; do
  if [[ "$SAMPLE_ID" == "sample_id" ]]; then
    continue
  fi

  SAMPLE_RAW_DIR="${RAW_ROOT}/${SAMPLE_ID}"
  SAMPLE_REVIEW_DIR="${REVIEW_ROOT}/${SAMPLE_ID}"
  mkdir -p "$SAMPLE_RAW_DIR" "$SAMPLE_REVIEW_DIR"

  EH_PREFIX="${SAMPLE_RAW_DIR}/${SAMPLE_ID}"
  EH_VCF="${EH_PREFIX}.vcf"
  EH_JSON="${EH_PREFIX}.json"
  EH_REALIGNED_BAM="${EH_PREFIX}_realigned.bam"
  EH_REALIGNED_SORTED_BAM="${EH_PREFIX}_realigned.sorted.bam"
  REVIEW_PREFIX="${SAMPLE_REVIEW_DIR}/${SAMPLE_ID}.PRNP_ORR"
  REVIEW_SVG="${REVIEW_PREFIX}.PRNP_ORR.svg"
  REVIEW_METRICS="${REVIEW_PREFIX}.metrics.tsv"
  REVIEW_PHASING="${REVIEW_PREFIX}.phasing.tsv"

  # Reuse existing sample-level caller outputs unless FORCE=1. This makes
  # intentional resumptions cheap while still protecting fresh clean reruns.
  if [[
    "$FORCE" != "1" &&
    -f "$EH_VCF" &&
    -f "$EH_JSON" &&
    (
      "$RUN_REVIEWER" != "1" ||
      -f "$EH_REALIGNED_BAM"
    )
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

  if [[ "$RUN_REVIEWER" == "1" ]]; then
    # ExpansionHunter emits the realigned BAM that REViewer consumes. Treat its
    # absence as a hard failure because downstream review would be incomplete.
    if [[ ! -f "$EH_REALIGNED_BAM" ]]; then
      echo "Expected realigned BAM missing for $SAMPLE_ID: $EH_REALIGNED_BAM" >&2
      exit 1
    fi

    # Sort and index the realigned BAM once so REViewer gets a stable input.
    if [[ "$FORCE" == "1" || ! -f "$EH_REALIGNED_SORTED_BAM" ]]; then
      samtools sort -o "$EH_REALIGNED_SORTED_BAM" "$EH_REALIGNED_BAM"
      samtools index "$EH_REALIGNED_SORTED_BAM"
    fi

    # Reuse completed review artifacts when possible; otherwise regenerate them
    # from the current caller outputs.
    if [[
      "$FORCE" != "1" &&
      -f "$REVIEW_SVG" &&
      -f "$REVIEW_METRICS" &&
      -f "$REVIEW_PHASING"
    ]]; then
      echo "Skipping REViewer for $SAMPLE_ID (outputs exist)."
    else
      echo "Running REViewer for $SAMPLE_ID"
      "${REVIEWER_CMD[@]}" \
        --reads "$EH_REALIGNED_SORTED_BAM" \
        --vcf "$EH_VCF" \
        --reference "$REPEAT_REF_FASTA" \
        --catalog "$PRNP_EH_CATALOG" \
        --locus PRNP_ORR \
        --output-prefix "$REVIEW_PREFIX" \
        </dev/null \
        >"${LOG_ROOT}/${SAMPLE_ID}.reviewer.stdout.log" \
        2>"${LOG_ROOT}/${SAMPLE_ID}.reviewer.stderr.log"
    fi
  fi
done 3<"$MANIFEST_TSV"

# ---------------------------------------------------------------------------
# Build cohort summary tables and verify the final live run
# ---------------------------------------------------------------------------

python "$SUMMARY_PY" \
  --manifest "$MANIFEST_TSV" \
  --results-root "$REPEAT_RESULTS_ROOT" \
  --reference-total-repeats "$PRNP_TOTAL_REFERENCE_REPEATS" \
  --variable-repeat-offset "$PRNP_VARIABLE_REPEAT_OFFSET" \
  --reviewer-enabled "$RUN_REVIEWER"

verify_latest_run "$MANIFEST_TSV" "$SAMPLE_COUNT"

echo "Repeat workflow complete."
