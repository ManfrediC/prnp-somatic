#!/bin/bash
# run_MF_loop_original.sh
#
# This script loops over samples in the /add/MosaicForecast/input/ directory,
# modifies the input_positions.tsv file so that the variant “name” field is set to the current sample name,
# runs MosaicForecast (MF) inside a Docker container,
# and then copies the six MF output files to a results directory, renaming them with the sample name.
#
# Changes from the "bwa" version:
#   - Input files now end in _original (e.g. CJD23_original.bam)
#   - RESULTS_DIR is set to /add/results/2025-02-10_MF_loop_original
#   - Sample names include the _original suffix.
#
# Before running this script, verify that:
#   - Your input files (BAM, BAM index, reference fasta, input_positions.tsv, and k24.umap.wg.bw) are in /add/MosaicForecast
#   - The Docker image "yanmei/mosaicforecast:0.0.1" is available (the script will pull it if needed)
#   - You have permission to run Docker (this script uses "sudo docker run" – adjust if needed)

set -euo pipefail

# ***** CONFIGURATION *****
# Base directory where MF files are stored.
MF_DIR="/add/MosaicForecast"
INPUT_DIR="${MF_DIR}/input"
OUTPUT_DIR="${MF_DIR}/output"

# Directory where results from the loop will be saved.
RESULTS_DIR="/add/results/2025-02-10_MF_loop_original"
mkdir -p "$RESULTS_DIR"

# Backup the original input_positions file if not already done.
ORIGINAL_POSITIONS="${INPUT_DIR}/input_positions.tsv.original"
if [ ! -f "$ORIGINAL_POSITIONS" ]; then
    cp "${INPUT_DIR}/input_positions.tsv" "$ORIGINAL_POSITIONS"
fi

# List of expected MF output files (produced in the OUTPUT_DIR by Phase.py).
output_files=("all_2x2table" "all_candidates" "all.merged.inforSNPs.pos" "all.phasing" "all.phasing_2by2" "multiple_inforSNPs.log")

# Pull the Docker image (if not already present)
sudo docker image pull yanmei/mosaicforecast:0.0.1

# ***** MAIN LOOP: Process each sample *****
# Loop over BAM files ending in _original.bam (ignoring the BAI files)
for bam_file in "${INPUT_DIR}"/*_original.bam; do
    # Extract the sample name (e.g. from /add/MosaicForecast/input/CJD23_original.bam, sample becomes "CJD23_original")
    sample=$(basename "$bam_file" .bam)
    echo "-------------------------------------"
    echo "Processing sample: $sample"
    echo "-------------------------------------"

    # Create a modified input_positions.tsv where the 6th column is replaced with the sample name.
    awk -v samp="$sample" 'BEGIN { OFS="\t" } { $NF = samp; print }' "$ORIGINAL_POSITIONS" \
        > "${INPUT_DIR}/input_positions.tsv"

    # Clean the MF output directory before running MF (remove all files and directories).
    rm -rf "${OUTPUT_DIR}"/*

    # ***** RUN MOSAICFORECAST via Docker *****
    sudo docker run \
         -v "${MF_DIR}":/MF \
         -w /MF \
         --rm yanmei/mosaicforecast:0.0.1 \
         python /usr/local/bin/Phase.py \
             input/ output/ \
             input/chr2_chr4_chr20.fasta \
             input/input_positions.tsv \
             20 \
             k24.umap.wg.bw \
             4 \
             bam

    # ***** SAVE RESULTS for the current sample *****
    sample_result_dir="${RESULTS_DIR}/${sample}"
    mkdir -p "$sample_result_dir"

    for file in "${output_files[@]}"; do
        src_file="${OUTPUT_DIR}/${file}"
        if [ -f "$src_file" ]; then
            dest_file="${sample_result_dir}/${sample}_${file}"
            cp "$src_file" "$dest_file"
            echo "Copied $(basename "$src_file") to $(basename "$dest_file")"
        else
            echo "Warning: Expected output file '$file' was not found for sample '$sample'."
        fi
    done

    echo "Finished processing sample: $sample"
    echo ""
done

echo "All samples processed successfully."
