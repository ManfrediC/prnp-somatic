#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${1:-manifest.tsv}"
OUT="${2:-manifest_qc.tsv}"

if [[ ! -f "$MANIFEST" ]]; then
  echo "ERROR: manifest not found: $MANIFEST" >&2
  exit 2
fi

expected_header=$'sample_id\tgroup\tbatch\tinput_dir'
read -r header < "$MANIFEST" || true
if [[ "$header" != "$expected_header" ]]; then
  echo "ERROR: manifest header mismatch." >&2
  echo "Expected: $expected_header" >&2
  echo "Found:    $header" >&2
  exit 2
fi

allowed_groups_regex='^(dilution|CJD|control)$'

# Candidate BAM naming styles
# We support both historical "short" names and full preprocessing names.
bam_candidate_suffixes=(
  ".bam"
  ".bwa.picard.markedDup.recal.bam"
)

picard_metrics_suffix=".bwa.picard.markedDup.metrics"

printf "sample_id\tgroup\tbatch\tinput_dir\tbam_path\tbam_style\tbam_exists\tbai_path\tbai_exists\tpicard_metrics_path\tpicard_metrics_exists\tpicard_metrics_required\terrors\n" > "$OUT"

while IFS=$'\t' read -r sample_id group batch input_dir || [[ -n "${sample_id:-}${group:-}${batch:-}${input_dir:-}" ]]; do
  if [[ -z "${sample_id}${group}${batch}${input_dir}" ]]; then
    continue
  fi

  errors=()

  if [[ -z "$sample_id" ]]; then errors+=("missing_sample_id"); fi
  if [[ -z "$group" ]]; then errors+=("missing_group"); fi
  if [[ -z "$batch" ]]; then errors+=("missing_batch"); fi
  if [[ -z "$input_dir" ]]; then errors+=("missing_input_dir"); fi

  if [[ -n "$group" && ! "$group" =~ $allowed_groups_regex ]]; then
    errors+=("invalid_group:${group}")
  fi

  if [[ -n "$input_dir" && ! -d "$input_dir" ]]; then
    errors+=("input_dir_not_found")
  fi

  input_dir="${input_dir%/}"

  # Resolve BAM path by trying known naming conventions in order.
  bam_path=""
  bam_style="none"
  for suf in "${bam_candidate_suffixes[@]}"; do
    candidate="${input_dir}/${sample_id}${suf}"
    if [[ -f "$candidate" ]]; then
      bam_path="$candidate"
      if [[ "$suf" == ".bam" ]]; then
        bam_style="short"
      else
        bam_style="long"
      fi
      break
    fi
  done

  bam_exists=0
  if [[ -n "$bam_path" ]]; then
    bam_exists=1
  else
    # Construct the first candidate path for reporting purposes
    bam_path="${input_dir}/${sample_id}${bam_candidate_suffixes[0]}"
    errors+=("missing_bam")
  fi

  # Index: accept either <bam>.bai OR <bam without .bam>.bai (covers *.recal.bai case).
  bai_exists=0
  bai_path="${bam_path}.bai"
  alt_bai_path="${bam_path%.bam}.bai"

  if [[ -f "$bai_path" ]]; then
    bai_exists=1
  elif [[ -f "$alt_bai_path" ]]; then
    bai_exists=1
    bai_path="$alt_bai_path"
  else
    errors+=("missing_bai")
  fi

  # Picard MarkDuplicates metrics are required for long-style BAMs (preprocessing outputs),
  # optional for short-style BAMs (legacy/manual BAMs).
  picard_metrics_path="${input_dir}/${sample_id}${picard_metrics_suffix}"
  picard_metrics_exists=0
  if [[ -f "$picard_metrics_path" ]]; then
    picard_metrics_exists=1
  fi

  picard_metrics_required=0
  if [[ "$bam_style" == "long" ]]; then
    picard_metrics_required=1
    if [[ "$picard_metrics_exists" -eq 0 ]]; then
      errors+=("missing_picard_metrics_required")
    fi
  fi

  err_str=""
  if [[ ${#errors[@]} -gt 0 ]]; then
    err_str="$(IFS=';'; echo "${errors[*]}")"
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%s\n" \
    "$sample_id" "$group" "$batch" "$input_dir" \
    "$bam_path" "$bam_style" "$bam_exists" \
    "$bai_path" "$bai_exists" \
    "$picard_metrics_path" "$picard_metrics_exists" "$picard_metrics_required" \
    "$err_str" >> "$OUT"
done < <(tail -n +2 "$MANIFEST")

# Failure definition:
# - missing BAM or BAI always fails
# - missing Picard metrics fails only if required (long-style)
# - any non-empty "errors" field fails
fail_count=$(awk -F'\t' '
  NR>1 {
    bam_ok=($7==1);
    bai_ok=($9==1);
    pic_req=($12==1);
    pic_ok=($11==1);
    if (!bam_ok || !bai_ok || (pic_req && !pic_ok) || $13!="") c++;
  }
  END { print c+0 }
' "$OUT")

total_count=$(awk -F'\t' 'NR>1 { c++ } END { print c+0 }' "$OUT")

echo "Validated $total_count manifest rows. Failures: $fail_count" >&2

if [[ "$fail_count" -gt 0 ]]; then
  exit 1
fi
exit 0
