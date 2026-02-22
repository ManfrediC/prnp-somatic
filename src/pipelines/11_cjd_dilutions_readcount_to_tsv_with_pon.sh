#!/usr/bin/env bash
# ============================================================
# CJD+dilutions readcount parsing (with PoN):
#   Stage 5: Convert bam-readcount text files into per-allele TSV metrics
# ============================================================
set -euo pipefail
shopt -s nullglob

# Keep caller-provided values so ENV_FILE does not silently override them.
CLI_DRY_RUN="${DRY_RUN-}"

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

# -----------------------
# Path defaults (repo-relative)
# -----------------------
MUTECT2_WITH_PON_OUT_ROOT="${MUTECT2_WITH_PON_OUT_ROOT:-runs/mutect2_cjd_dilutions_with_pon}"
if [[ "$MUTECT2_WITH_PON_OUT_ROOT" == run/* ]]; then
  MUTECT2_WITH_PON_OUT_ROOT="runs/${MUTECT2_WITH_PON_OUT_ROOT#run/}"
fi

MUTECT2_WITH_PON_READCOUNT_ROOT="${MUTECT2_WITH_PON_READCOUNT_ROOT:-$MUTECT2_WITH_PON_OUT_ROOT}"
if [[ "$MUTECT2_WITH_PON_READCOUNT_ROOT" == run/* ]]; then
  MUTECT2_WITH_PON_READCOUNT_ROOT="runs/${MUTECT2_WITH_PON_READCOUNT_ROOT#run/}"
fi

WITH_PON_READCOUNT_GROUPS="${WITH_PON_READCOUNT_GROUPS:-${WITH_PON_GROUPS:-cjd dilutions}}"
READCOUNT_TO_TSV_PY="${READCOUNT_TO_TSV_PY:-src/pipelines/4_readcount_to_tsv.py}"


to_abs() {
  local p="$1"
  # Allow both repo-relative and absolute paths in ENV overrides.
  if [[ "$p" = /* ]]; then
    echo "$p"
  else
    echo "$REPO_ROOT/$p"
  fi
}

READCOUNT_ROOT="$(to_abs "$MUTECT2_WITH_PON_READCOUNT_ROOT")"
READCOUNT_TO_TSV_PY="$(to_abs "$READCOUNT_TO_TSV_PY")"

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

# -----------------------
# Tool + file checks
# -----------------------
have python3 || die "python3 not in PATH"
[[ -d "$READCOUNT_ROOT" ]] || die "Missing root directory: $READCOUNT_ROOT"
[[ -s "$READCOUNT_TO_TSV_PY" ]] || die "Missing parser script: $READCOUNT_TO_TSV_PY"

read -r -a groups <<< "$WITH_PON_READCOUNT_GROUPS"
[[ "${#groups[@]}" -gt 0 ]] || die "WITH_PON_READCOUNT_GROUPS is empty"

echo "== CJD+dilutions readcount parsing (with PoN) =="
echo "Repo root:        $REPO_ROOT"
echo "READCOUNT_ROOT:   $READCOUNT_ROOT"
echo "WITH_PON_GROUPS:  $WITH_PON_READCOUNT_GROUPS"
echo "PARSER:           $READCOUNT_TO_TSV_PY"
echo

for group in "${groups[@]}"; do
  echo "=== Group: $group ==="

  READCOUNTS_DIR="$READCOUNT_ROOT/$group/readcount_qc/readcounts"
  METRICS_DIR="$READCOUNT_ROOT/$group/readcount_qc/metrics"

  [[ -d "$READCOUNTS_DIR" ]] || die "Missing readcount directory for group '$group': $READCOUNTS_DIR"
  mkdir -p "$METRICS_DIR"

  readcount_files=( "$READCOUNTS_DIR"/*.txt )
  [[ "${#readcount_files[@]}" -gt 0 ]] || die "No readcount files found for group '$group' in: $READCOUNTS_DIR"

  run python3 "$READCOUNT_TO_TSV_PY" --input-dir "$READCOUNTS_DIR" --output-dir "$METRICS_DIR"

  if [[ "$DRY_RUN" == "0" ]]; then
    missing=()
    for txt in "${readcount_files[@]}"; do
      sample="$(basename "$txt" .txt)"
      [[ -s "$METRICS_DIR/${sample}_metrics.tsv" ]] || missing+=( "$sample" )
    done
    [[ "${#missing[@]}" -eq 0 ]] || die "Final completeness check failed for group '$group'; missing metrics TSVs for: ${missing[*]}"
  fi

  echo "=== Group finished: $group ==="
  echo
done

echo "=== CJD+dilutions readcount parsing (with PoN) finished ==="
