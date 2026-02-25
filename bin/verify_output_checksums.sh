#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

MANIFEST="doc/reproducibility/final_outputs_manifest.tsv"
CHECKSUMS="doc/reproducibility/final_outputs.sha256"
MODE="write"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --checksums) CHECKSUMS="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

[[ -f "$MANIFEST" ]] || { echo "Manifest not found: $MANIFEST" >&2; exit 1; }
[[ "$MODE" =~ ^(write|check)$ ]] || { echo "--mode must be write|check" >&2; exit 2; }

tmp_file="$(mktemp)"
missing_required=0

while IFS=$'\t' read -r workflow path required notes; do
  [[ -z "${workflow}" ]] && continue
  [[ "${workflow}" == "workflow" ]] && continue
  [[ "${workflow:0:1}" == "#" ]] && continue

  if [[ -f "$path" ]]; then
    sha256sum "$path" >> "$tmp_file"
  else
    if [[ "$required" == "yes" ]]; then
      echo "MISSING REQUIRED: $path ($workflow)" >&2
      missing_required=1
    else
      echo "MISSING OPTIONAL: $path ($workflow)" >&2
    fi
  fi
done < "$MANIFEST"

if [[ "$missing_required" -ne 0 ]]; then
  rm -f "$tmp_file"
  exit 1
fi

if [[ "$MODE" == "write" ]]; then
  mv "$tmp_file" "$CHECKSUMS"
  echo "Wrote checksums: $CHECKSUMS"
  exit 0
fi

[[ -f "$CHECKSUMS" ]] || { echo "Checksums file not found: $CHECKSUMS" >&2; rm -f "$tmp_file"; exit 1; }
if diff -u "$CHECKSUMS" "$tmp_file"; then
  echo "Checksum verification passed"
  rm -f "$tmp_file"
else
  echo "Checksum verification failed" >&2
  rm -f "$tmp_file"
  exit 1
fi
