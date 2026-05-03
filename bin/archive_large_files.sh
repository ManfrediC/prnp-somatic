#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash bin/archive_large_files.sh inventory [options]
  bash bin/archive_large_files.sh compress [options]
  bash bin/archive_large_files.sh restore [options]

Actions:
  inventory   List files that match the archive rules
  compress    Compress matching files in place
  restore     Restore .xz/.gz archives to their original filenames

Options:
  --root PATH           Repository root to scan (default: current directory)
  --dirs CSV            Comma-separated directories to scan
                        (default: runs,results,resources,raw/fastq)
  --exts CSV            Comma-separated extensions to target, without dots
                        (default: bam,bai,sam,fa,fasta,fna,gtf,vcf,tsv)
  --min-size-mb N       Minimum file size in MiB (default: 100)
  --format xz|gz        Archive format for compress (default: xz)
  --force               Recompress even if an archive already exists
  --dry-run             Show actions without modifying files
  --help                Show this message

Examples:
  bash bin/archive_large_files.sh inventory
  bash bin/archive_large_files.sh compress --dirs runs,results --exts bam,bai,sam
  bash bin/archive_large_files.sh restore --dirs runs,results
EOF
}

die() {
  printf 'Error: %s\n' "$*" >&2
  exit 1
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

ROOT="."
DIRS_CSV="runs,results,resources,raw/fastq"
EXTS_CSV="bam,bai,sam,fa,fasta,fna,gtf,vcf,tsv"
MIN_SIZE_MB=100
FORMAT="xz"
FORCE=0
DRY_RUN=0

ACTION="${1:-}"
if [[ -z "${ACTION}" ]]; then
  usage
  exit 1
fi
shift

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root)
      ROOT="$2"
      shift 2
      ;;
    --dirs)
      DIRS_CSV="$2"
      shift 2
      ;;
    --exts)
      EXTS_CSV="$2"
      shift 2
      ;;
    --min-size-mb)
      MIN_SIZE_MB="$2"
      shift 2
      ;;
    --format)
      FORMAT="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

case "${ACTION}" in
  inventory|compress|restore) ;;
  *)
    die "unknown action: ${ACTION}"
    ;;
esac

case "${FORMAT}" in
  xz|gz) ;;
  *)
    die "--format must be xz or gz"
    ;;
esac

[[ "${MIN_SIZE_MB}" =~ ^[0-9]+$ ]] || die "--min-size-mb must be an integer"
require_cmd find
require_cmd stat
require_cmd awk

if [[ "${ACTION}" == "compress" ]]; then
  if [[ "${FORMAT}" == "xz" ]]; then
    require_cmd xz
  else
    require_cmd gzip
  fi
fi

if [[ "${ACTION}" == "restore" ]]; then
  require_cmd xz
  require_cmd gzip
fi

IFS=',' read -r -a SCAN_DIRS <<< "${DIRS_CSV}"
IFS=',' read -r -a TARGET_EXTS <<< "${EXTS_CSV}"
MIN_SIZE_BYTES=$((MIN_SIZE_MB * 1024 * 1024))

trim() {
  local value="$1"
  value="${value#"${value%%[![:space:]]*}"}"
  value="${value%"${value##*[![:space:]]}"}"
  printf '%s' "${value}"
}

join_by() {
  local sep="$1"
  shift
  local out=""
  local first=1
  local item
  for item in "$@"; do
    if [[ ${first} -eq 1 ]]; then
      out="${item}"
      first=0
    else
      out="${out}${sep}${item}"
    fi
  done
  printf '%s' "${out}"
}

human_mb() {
  awk -v bytes="$1" 'BEGIN { printf "%.1f MiB", bytes / 1024 / 1024 }'
}

should_skip_compressed_name() {
  local path="$1"
  case "${path}" in
    *.gz|*.xz|*.zip|*.bz2|*.zst|*.bgz|*.cram|*.crai|*.tar|*.tar.gz)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

matches_target_extension() {
  local path="$1"
  local ext
  for ext in "${TARGET_EXTS[@]}"; do
    ext="$(trim "${ext}")"
    [[ -z "${ext}" ]] && continue
    if [[ "${path}" == *.${ext} ]]; then
      return 0
    fi
  done
  return 1
}

emit_candidates() {
  local dir
  local scan_dir
  for dir in "${SCAN_DIRS[@]}"; do
    scan_dir="$(trim "${dir}")"
    [[ -z "${scan_dir}" ]] && continue
    scan_dir="${ROOT%/}/${scan_dir}"
    [[ -d "${scan_dir}" ]] || continue
    while IFS= read -r path; do
      [[ -f "${path}" ]] || continue
      should_skip_compressed_name "${path}" && continue
      matches_target_extension "${path}" || continue
      local size
      size="$(stat -c %s "${path}")"
      (( size >= MIN_SIZE_BYTES )) || continue
      printf '%s\t%s\n' "${size}" "${path}"
    done < <(find "${scan_dir}" -type f | sort)
  done
}

emit_archives() {
  local dir
  local scan_dir
  for dir in "${SCAN_DIRS[@]}"; do
    scan_dir="$(trim "${dir}")"
    [[ -z "${scan_dir}" ]] && continue
    scan_dir="${ROOT%/}/${scan_dir}"
    [[ -d "${scan_dir}" ]] || continue
    while IFS= read -r path; do
      [[ -f "${path}" ]] || continue
      case "${path}" in
        *.xz|*.gz)
          local size
          size="$(stat -c %s "${path}")"
          printf '%s\t%s\n' "${size}" "${path}"
          ;;
      esac
    done < <(find "${scan_dir}" -type f | sort)
  done
}

inventory_candidates() {
  local count=0
  local total_bytes=0
  while IFS=$'\t' read -r size path; do
    [[ -n "${path:-}" ]] || continue
    printf '%10s  %s\n' "$(human_mb "${size}")" "${path#${ROOT%/}/}"
    count=$((count + 1))
    total_bytes=$((total_bytes + size))
  done < <(emit_candidates)
  printf 'Candidates: %d files, %s total\n' "${count}" "$(human_mb "${total_bytes}")"
}

compress_file() {
  local path="$1"
  local archive_path
  if [[ "${FORMAT}" == "xz" ]]; then
    archive_path="${path}.xz"
  else
    archive_path="${path}.gz"
  fi

  if [[ -e "${archive_path}" && ${FORCE} -ne 1 ]]; then
    printf 'skip    %s (archive exists)\n' "${path#${ROOT%/}/}"
    return
  fi

  if [[ ${DRY_RUN} -eq 1 ]]; then
    printf 'would compress  %s -> %s\n' "${path#${ROOT%/}/}" "${archive_path#${ROOT%/}/}"
    return
  fi

  if [[ "${FORMAT}" == "xz" ]]; then
    xz -T0 -z -f "${path}"
  else
    gzip -f "${path}"
  fi
  printf 'compressed  %s\n' "${archive_path#${ROOT%/}/}"
}

restore_file() {
  local path="$1"
  if [[ ${DRY_RUN} -eq 1 ]]; then
    printf 'would restore  %s\n' "${path#${ROOT%/}/}"
    return
  fi

  case "${path}" in
    *.xz)
      xz -d -f "${path}"
      ;;
    *.gz)
      gzip -d -f "${path}"
      ;;
  esac
  printf 'restored   %s\n' "${path#${ROOT%/}/}"
}

compress_candidates() {
  local count=0
  while IFS=$'\t' read -r _size path; do
    [[ -n "${path:-}" ]] || continue
    compress_file "${path}"
    count=$((count + 1))
  done < <(emit_candidates)
  printf 'Processed %d candidate files\n' "${count}"
}

restore_archives() {
  local count=0
  while IFS=$'\t' read -r _size path; do
    [[ -n "${path:-}" ]] || continue
    restore_file "${path}"
    count=$((count + 1))
  done < <(emit_archives)
  printf 'Processed %d archive files\n' "${count}"
}

printf 'Action: %s\n' "${ACTION}"
printf 'Root: %s\n' "$(cd "${ROOT}" && pwd)"
printf 'Dirs: %s\n' "$(join_by ', ' "${SCAN_DIRS[@]}")"
printf 'Exts: %s\n' "$(join_by ', ' "${TARGET_EXTS[@]}")"
printf 'Min size: %s\n' "$(human_mb "${MIN_SIZE_BYTES}")"
if [[ "${ACTION}" == "compress" ]]; then
  printf 'Format: %s\n' "${FORMAT}"
fi

case "${ACTION}" in
  inventory)
    inventory_candidates
    ;;
  compress)
    compress_candidates
    ;;
  restore)
    restore_archives
    ;;
esac
