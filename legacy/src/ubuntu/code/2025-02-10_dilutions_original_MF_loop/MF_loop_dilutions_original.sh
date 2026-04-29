#!/bin/bash
# 
# This script loops over samples in the /add/MosaicForecast/input/ directory,
# where sample files are named with the suffix "_original" (e.g. NA100_1to10_original.bam, etc.).
#
# For each sample, the script:
#   1. Modifies the input_positions.tsv file (based on an original backup) so that the 6th column is
#      replaced with the current sample name.
#   2. Clears the MF output directory.
#   3. Runs MosaicForecast (Phase.py) inside a Docker container.
#   4. Copies the expected output files to a results directory, renaming them to include the sample name.
#
# The output files are stored under:
#   /add/results/2025-02-10_MF_loop_original_dilutions/
#
# For example, if the sample name is "NA99A1_undil_original", the output files will be copied as:
#   NA99A1_undil_original_all_2x2table, NA99A1_undil_original_all_candidates, etc.
#
# Before running, ensure:
#   - Your input files (BAM, BAM index, reference FASTA, input_positions.tsv, k24.umap.wg.bw) are in /add/MosaicForecast.
#   - The Docker image "yanmei/mosaicforecast:0.0.1" is available.
#   - You have permission to run Docker (this script uses "sudo docker run").

set -euo pipefail

# ***** CONFIGURATION *****
MF_DIR="/add/MosaicForecast"
INPUT_DIR="${MF_DIR}/input"
OUTPUT_DIR="${MF_DIR}/output"
RESULTS_DIR="/add/results/2025-02-10_MF_loop_original_dilutions"
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
# Loop over BAM files ending with _original.bam
for bam_file in "${INPUT_DIR}"/*_original.bam; do
    # Extract the sample name (e.g. from NA100_1to10_original.bam, sample becomes "NA100_1to10_original")
    sample=$(basename "$bam_file" .bam)
    echo "-------------------------------------"
    echo "Processing sample: $sample"
    echo "-------------------------------------"

    # Create a modified input_positions.tsv where the 6th column is replaced with the sample name.
    awk -v samp="$sample" 'BEGIN { OFS="\t" } { $NF = samp; print }' "$ORIGINAL_POSITIONS" > "${INPUT_DIR}/input_positions.tsv"

    # Clean the MF output directory (remove any files/directories from a previous run)
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
