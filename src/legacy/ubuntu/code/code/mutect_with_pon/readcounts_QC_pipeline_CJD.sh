#!/usr/bin/env bash
# ============================================================
# bam-readcount metrics pipeline (with-PoN project, no dilutions)
#   1) Build BED of variant sites from VCFs
#   2) Copy recalibrated BAMs (RG / SM already present)
#   3) Strip PCR duplicates (flag 0x400)
#   4) Run bam-readcount per sample
#   5) Parse readcount-to-metrics TSV  (ALT_COUNT  ALT_BQ  ALT_MQ)
# Each step skips files that already exist, so you can resume.
# ============================================================
set -euo pipefail
shopt -s nullglob

ROOT=/add/results/with_PoN
FASTA=/home/mcarta/databases/chr2_chr4_chr20.fasta

VCF_SAMPLES=$ROOT/annot
BAM_SAMPLES=/add/seq_data/2025-02-07_samples_bwa

# create only the directories we still need (no “…_dil”)
mkdir -p "$ROOT"/{beds,bam_work,bam_nodup,readcounts,metrics}

# Helper: get core sample ID (strip after first dot)
base_id () { basename "$1" | cut -d'.' -f1; }

# ----------------------------------------
# 1) Build BED of variant sites from VCFs
# ----------------------------------------
make_bed() {
  local vcf=$1 bed=$2 tag=$3
  [[ -s $bed ]] && { echo "[1] BED  SKIP  $tag"; return; }
  echo "[1] BED  $tag"
  bcftools query -f '%CHROM\t%POS0\t%POS\t%ALT\n' "$vcf" > "$bed"
}
for vcf in "$VCF_SAMPLES"/*.vcf.gz; do
  id=$(base_id "$vcf")
  make_bed "$vcf" "$ROOT/beds/${id}.bed" "$id"
done

# ----------------------------------------
# 2) Copy BAMs (already have RG/SM tags)
# ----------------------------------------
copy_bam() {
  local src=$1 dest=$2 tag=$3
  [[ -f $dest ]] && { echo "[2] COPY SKIP  $tag"; return; }
  echo "[2] COPY $tag"
  ln -s "$src"       "$dest"
  ln -s "${src}.bai" "${dest}.bai"
}
for bam in "$BAM_SAMPLES"/*.markedDup.recal.bam; do
  id=$(base_id "$bam")
  copy_bam "$bam" "$ROOT/bam_work/${id}.bam" "$id"
done

# ----------------------------------------
# 3) Strip PCR duplicates (flag 0x400)
# ----------------------------------------
strip_dup() {
  local inBam=$1 outBam=$2 tag=$3
  [[ -f $outBam ]] && { echo "[3] DEDUP SKIP $tag"; return; }
  echo "[3] DEDUP $tag"
  samtools view -@4 -b -F 0x400 "$inBam" > "$outBam"
  samtools index -@4 "$outBam"
}
for bam in "$ROOT"/bam_work/*.bam; do
  id=$(base_id "$bam")
  strip_dup "$bam" "$ROOT/bam_nodup/${id}.nodup.bam" "$id"
done

# ----------------------------------------
# 4) Run bam-readcount per sample
# ----------------------------------------
run_brc() {
  local bam=$1 bed=$2 out=$3 tag=$4
  [[ -f $out ]] && { echo "[4] RC   SKIP $tag"; return; }
  echo "[4] RC   $tag"
  bam-readcount -f "$FASTA" -l "$bed" "$bam" > "$out"
}
for bam in "$ROOT"/bam_nodup/*.nodup.bam; do
  id=$(base_id "$bam")
  run_brc "$bam" "$ROOT/beds/${id}.bed" "$ROOT/readcounts/${id}.txt" "$id"
done

# ----------------------------------------
# 5) Parse readcount to full-metrics TSV
# ----------------------------------------
python3 "/add/code/mutect_with_pon/readcount_to_TSV_CJD.py"

echo "DONE readcount pipeline"
