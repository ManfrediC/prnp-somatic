#!/usr/bin/env bash
# Generate or verify checksums for declared final outputs.
#
# Source of truth:
# - manifest: doc/reproducibility/final_outputs_manifest.tsv
# - checksums: doc/reproducibility/final_outputs.sha256
#
# Behavior:
# - in write mode, computes checksums for present outputs and writes checksum file
# - in check mode, recomputes checksums and diffs against stored checksum file
# - fails fast if any manifest entry marked required=yes is missing
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

# Build checksums in a temporary file first, then atomically move/diff.
tmp_file="$(mktemp)"
missing_required=0

while IFS=$'\t' read -r workflow path required notes; do
  # Skip blank/header/comment lines.
  [[ -z "${workflow}" ]] && continue
  [[ "${workflow}" == "workflow" ]] && continue
  [[ "${workflow:0:1}" == "#" ]] && continue

  # Compute checksum for each file that exists.
  if [[ -f "$path" ]]; then
    sha256sum "$path" >> "$tmp_file"
  else
    # Required missing files are fatal; optional missing files are reported only.
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
  # Overwrite stored checksums only after successful computation.
  mv "$tmp_file" "$CHECKSUMS"
  echo "Wrote checksums: $CHECKSUMS"
  exit 0
fi

[[ -f "$CHECKSUMS" ]] || { echo "Checksums file not found: $CHECKSUMS" >&2; rm -f "$tmp_file"; exit 1; }
# In check mode, diff recomputed checksums against committed checksum file.
if diff -u "$CHECKSUMS" "$tmp_file"; then
  echo "Checksum verification passed"
  rm -f "$tmp_file"
else
  echo "Checksum verification failed" >&2
  rm -f "$tmp_file"
  exit 1
fi
