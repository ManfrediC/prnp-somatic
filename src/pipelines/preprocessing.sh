#!/usr/bin/env bash
set -euo pipefail

# -----------------------
# Batch selection (toggle by commenting/uncommenting)
# -----------------------
# This array is the single run-selection switch for preprocessing batches.
BATCHES=(
  CJD_16_samples
  CJD_8_samples
  first_CJD_seq
  sequencing_of_dilutions
)

# -----------------------
# Repo root + config
# -----------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

ENV_FILE="${ENV_FILE:-$REPO_ROOT/config/preprocessing.env}"

# read ENV_FILE, variables defined there become available in this script
source "$ENV_FILE"

# Defaults (if not set in preprocessing.env)
SAMPLES_TSV="${SAMPLES_TSV:-config/preprocessing_samples.tsv}"
FASTQ_DIR="${FASTQ_DIR:-fastq}"
RUNS_DIR="${RUNS_DIR:-runs/preprocessing}"
FINAL_BAM_DIR="${FINAL_BAM_DIR:-results/final_bam}"
THREADS="${THREADS:-8}"
JAVA_MEM_GB="${JAVA_MEM_GB:-8}"
FORCE="${FORCE:-0}"
DRY_RUN="${DRY_RUN:-0}"

# Normalize legacy /run paths to /runs
if [[ "$RUNS_DIR" == run/* ]]; then
  RUNS_DIR="runs/${RUNS_DIR#run/}"
fi

REF_FASTA="${REF_FASTA:-resources/chr2_chr4_chr20.fasta}"
DBSNP_VCF="${DBSNP_VCF:-resources/dbsnp_146.hg38.vcf.gz}"
MILLS_VCF="${MILLS_VCF:-resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz}"
ADAPTERS_FA="${ADAPTERS_FA:-resources/TruSeq3-PE.fa}"

# Convert relative paths to absolute paths
# Keep path normalisation centralised so downstream commands always see absolute inputs.
SAMPLES_TSV="$REPO_ROOT/$SAMPLES_TSV"
RUNS_DIR="$REPO_ROOT/$RUNS_DIR"
FINAL_BAM_DIR="$REPO_ROOT/$FINAL_BAM_DIR"
REF_FASTA="$REPO_ROOT/$REF_FASTA"
DBSNP_VCF="$REPO_ROOT/$DBSNP_VCF"
MILLS_VCF="$REPO_ROOT/$MILLS_VCF"
ADAPTERS_FA="$REPO_ROOT/$ADAPTERS_FA"

# trimming string
TRIMM_STRING="ILLUMINACLIP:${ADAPTERS_FA}:1:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:25 AVGQUAL:30 MINLEN:80"


# make final BAM directory
mkdir -p "$FINAL_BAM_DIR"

# -----------------------
# Minimal helper
# -----------------------
run() {
  # DRY_RUN prints the resolved command exactly as it would execute.
  if [[ "$DRY_RUN" == "1" ]]; then
    echo "+ $*"
  else
    eval "$@"
  fi
}

# -----------------------
# Main loop
# -----------------------

# Read and validate header
# Expected schema is fixed so pipeline parsing remains deterministic.
header="$(grep -v '^#' "$SAMPLES_TSV" | head -n 1)"
[[ "$header" == $'batch_id\tsample_id\tr1\tr2' ]] || {
  echo "ERROR: unexpected header in $SAMPLES_TSV: $header" >&2
  exit 1
}

# Ensure publish directory exists
mkdir -p "$FINAL_BAM_DIR"

# Iterate TSV rows (skip header + comments)
tail -n +2 "$SAMPLES_TSV" | grep -v '^#' | while IFS=$'\t' read -r batch sample r1 r2; do
  [[ -n "${batch// }" ]] || continue #skip blanks/whitespace

  # Filter to selected batches (loop)
  # This enforces BATCHES selection without editing the input TSV.
  selected=0
  for b in "${BATCHES[@]}"; do
    if [[ "$batch" == "$b" ]]; then selected=1; break; fi
  done
  [[ "$selected" == "1" ]] || continue

  # Build FASTQ paths (TSV contains repo-relative paths)
  r1_path="$REPO_ROOT/$r1"
  r2_path="$REPO_ROOT/$r2"
  
  # Build sample path and log directory
  sample_dir="$RUNS_DIR/$batch/$sample"
  log_dir="$sample_dir/logs"
  mkdir -p "$log_dir"

  # Define filenames to be used inside sample_dir
  trim_r1="$sample_dir/${sample}.R1.trimmed.fastq.gz"
  trim_r2="$sample_dir/${sample}.R2.trimmed.fastq.gz"
  unp_r1="$sample_dir/${sample}.R1.unpaired.fastq.gz"
  unp_r2="$sample_dir/${sample}.R2.unpaired.fastq.gz"

  sorted_bam="$sample_dir/${sample}.bwa.sorted.bam"
  sorted_bai="${sorted_bam}.bai"

  rg_bam="$sample_dir/${sample}.bwa.picard.bam"
  md_bam="$sample_dir/${sample}.bwa.picard.markedDup.bam"
  md_metrics="$sample_dir/${sample}.bwa.picard.markedDup.metrics"

  bqsr_table="$sample_dir/${sample}.bwa.recal_data.table"
  final_bam="$sample_dir/${sample}.bwa.picard.markedDup.recal.bam"
  final_bai="${final_bam}.bai"

  publish_bam="$FINAL_BAM_DIR/${sample}.bam"
  publish_bai="$FINAL_BAM_DIR/${sample}.bam.bai"

  # Simple overwrite policy: If final_bam exists and is non-empty, skip (unless FORCE is 1)
  if [[ -s "$final_bam" && "$FORCE" != "1" ]]; then
    echo "SKIP: $batch/$sample (final BAM exists)"
    continue
  fi

  # Write per-sample metadata file (updated each run)
  {
    echo "timestamp: $(date -Is)"
    echo "repo_root: $REPO_ROOT"
    echo "batch_id: $batch"
    echo "sample_id: $sample"
    echo "r1: $r1"
    echo "r2: $r2"
    echo "threads: $THREADS"
    echo "java_mem_gb: $JAVA_MEM_GB"
    echo "git_commit: $(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo NA)"
  } > "$sample_dir/RUN_META.txt"
  
  # Progress message for the terminal
  echo "== $batch / $sample =="

  if [[ "$DRY_RUN" == "1" ]]; then
    echo "  R1: $r1"
    echo "  R2: $r2"
    echo "  out: $final_bam -> $publish_bam"
    continue
  fi
  
  # log directory and creation
  mkdir -p "$log_dir"

  {
    echo "timestamp: $(date -Is)"
    echo "batch_id: $batch"
    echo "sample_id: $sample"
    echo "r1: $r1"
    echo "r2: $r2"
    echo "threads: $THREADS"
    echo "java_mem_gb: $JAVA_MEM_GB"
    echo "git_commit: $(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo NA)"
  } > "$sample_dir/RUN_META.txt"

  # -----------------------
  # Stage 1: Trimmomatic
  # -----------------------
  echo "Starting Trimmomatic for $sample..."
    
  if [[ ! -s "$trim_r1" || ! -s "$trim_r2" ]]; then # checks whether trim_r1 and trim_r2 exist
    trimmomatic PE -threads "$THREADS" -phred33 \
      "$r1_path" "$r2_path" \
      "$trim_r1" "$unp_r1" \
      "$trim_r2" "$unp_r2" \
      $TRIMM_STRING \
      > "$log_dir/trim.log" 2>&1
  else
    echo "  skip trim (exists)"
  fi

  # -----------------------  
  # Stage 2:  bwa mem -> samtools sort (+ index)
  # -----------------------
  echo "Starting bwa mem for $sample..."
  if [[ ! -s "$sorted_bam" ]]; then
    bwa mem -t "$THREADS" "$REF_FASTA" "$trim_r1" "$trim_r2" \
      2> "$log_dir/bwa.stderr.log" \
    | { echo "Starting samtools sort for $sample..." >&2
        samtools sort -@ "$THREADS" -o "$sorted_bam"
      } > "$log_dir/samtools_sort.log" 2>&1
  
    samtools index "$sorted_bam" > "$log_dir/samtools_index.log" 2>&1
  else
    echo "  skip align+sort (exists)"
  fi

  # -----------------------
  # Stage 3: add read groups
  # -----------------------
  echo "AddOrReplaceReadGroups (Picard) for $sample..."
  
  if [[ ! -s "$rg_bam" ]]; then
    picard AddOrReplaceReadGroups \
      I="$sorted_bam" \
      O="$rg_bam" \
      SORT_ORDER=coordinate \
      RGID="$sample" \
      RGLB=Paired_end \
      RGPL=illumina \
      RGSM="$sample" \
      RGPU=illumina \
      USE_JDK_DEFLATER=true \
      USE_JDK_INFLATER=true \
      > "$log_dir/add_rg.log" 2>&1
  else
    echo "  skip addRG (exists)"
  fi

  # -----------------------
  # Stage 4: mark duplicates (placeholder)
  # -----------------------
  echo "MarkDuplicates (Picard) for $sample..."
  
  if [[ ! -s "$md_bam" ]]; then
    picard MarkDuplicates \
      I="$rg_bam" \
      O="$md_bam" \
      M="$md_metrics" \
      USE_JDK_DEFLATER=true \
      USE_JDK_INFLATER=true \
      > "$log_dir/markdup.log" 2>&1
  else
    echo "  skip markdup (exists)"
  fi

  # -----------------------
  # Stage 5: BaseRecalibrator
  # -----------------------
  echo "BaseRecalibrator for $sample..."
  
  if [[ ! -s "$bqsr_table" ]]; then
    gatk --java-options "-Xmx${JAVA_MEM_GB}g" BaseRecalibrator \
      -I "$md_bam" \
      -R "$REF_FASTA" \
      --known-sites "$DBSNP_VCF" \
      --known-sites "$MILLS_VCF" \
      -O "$bqsr_table" \
      > "$log_dir/base_recalibrator.log" 2>&1
  else
    echo "  skip BQSR table (exists)"
  fi
  
  # ApplyBQSR (+ index)
  if [[ ! -s "$final_bam" ]]; then
    gatk --java-options "-Xmx${JAVA_MEM_GB}g" ApplyBQSR \
      -R "$REF_FASTA" \
      -I "$md_bam" \
      --bqsr-recal-file "$bqsr_table" \
      -O "$final_bam" \
      > "$log_dir/apply_bqsr.log" 2>&1
    samtools index "$final_bam" > "$log_dir/index_final.log" 2>&1
  else
    echo "  skip final BAM (exists)"
  fi
  
  # -----------------------
  # Stage 6: publish final BAM
  # -----------------------
 
  echo "Finalising BAM for $sample..."
 
  if [[ -s "$final_bam" && -s "$final_bai" ]]; then
    cp -f "$final_bam" "$publish_bam"
    cp -f "$final_bai" "$publish_bai"
    echo "  published -> $publish_bam"
  else
    echo "ERROR: final outputs missing for $batch/$sample" >&2
    exit 1
  fi
  
done
