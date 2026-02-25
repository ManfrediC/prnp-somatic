#!/usr/bin/env bash
set -euo pipefail

# --- configure ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
JUNCTION_RES_DIR="${REPO_ROOT}/resources/junctions"
cd "$REPO_ROOT"

CONFIG_FILE="${CONFIG_FILE:-${REPO_ROOT}/config/junctions.env}"
if [ -f "$CONFIG_FILE" ]; then
  # shellcheck source=/dev/null
  source "$CONFIG_FILE"
fi

SAMPLEDIR="${SAMPLEDIR:-${REPO_ROOT}/results/final_bam}"

BED="${BED:-${JUNCTION_RES_DIR}/PRNP.pad1kb.hg38.bed}"
JUNC_FA="${JUNC_FA:-${JUNCTION_RES_DIR}/prnp_junctions.fa}"
THREADS="${THREADS:-8}"

# --- output directory ---
OUTDIR="${OUTDIR:-${REPO_ROOT}/results/junctions/junction_align}"
mkdir -p "$OUTDIR"

if [ ! -d "$SAMPLEDIR" ]; then
  echo "Sample BAM directory not found: $SAMPLEDIR" >&2
  exit 1
fi
if [ ! -f "$BED" ]; then
  echo "BED file not found: $BED" >&2
  exit 1
fi
if [ ! -f "$JUNC_FA" ]; then
  echo "Junction FASTA not found: $JUNC_FA" >&2
  exit 1
fi

# 1) Build BWA index for the junction reference (once per JUNC_FA)
if [ ! -e "${JUNC_FA}.bwt" ]; then
  bwa index "$JUNC_FA"
fi

# 2) Loop over all CJD* and Ctrl* BAMs
for BAM in "$SAMPLEDIR"/CJD*.bam "$SAMPLEDIR"/Ctrl*.bam; do
  [ -e "$BAM" ] || continue

  SAMPLE=$(basename "$BAM" .bam)
  echo "Processing sample: $SAMPLE"

  # 2a) Extract PRNP-window reads and name-sort
  samtools view -b -L "$BED" -F 0x400 "$BAM" \
    | samtools sort -n -o "$OUTDIR/${SAMPLE}.PRNP.qname.bam"

  # 2b) Convert to FASTQ
  samtools fastq \
    -1 "$OUTDIR/${SAMPLE}.PRNP.R1.fq.gz" \
    -2 "$OUTDIR/${SAMPLE}.PRNP.R2.fq.gz" \
    -0 /dev/null -s /dev/null -n "$OUTDIR/${SAMPLE}.PRNP.qname.bam"

  # 2c) Align to junction reference and index
  bwa mem -t "$THREADS" "$JUNC_FA" \
    "$OUTDIR/${SAMPLE}.PRNP.R1.fq.gz" "$OUTDIR/${SAMPLE}.PRNP.R2.fq.gz" \
    | samtools sort -o "$OUTDIR/${SAMPLE}.PRNP.toJunc.bam"

  samtools index "$OUTDIR/${SAMPLE}.PRNP.toJunc.bam"
done
