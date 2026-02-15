#!/bin/bash
# run_mutect_all.sh
# This script sets the working directory to /add/code/2025-02-07_aggressive_pipeline/
# and then runs the following Mutect2 scripts sequentially:
#   2025-02-08_mutect_dilutions_bwa.sh
#   2025-02-08_mutect_all_samples_bwa.sh

# Change to the directory containing the Mutect2 scripts
cd /add/code/mutect_bwa/ || { echo "Failed to change directory"; exit 1; }

# Run the dilution Mutect2 script; if it succeeds, run the sample Mutect2 script.
bash 2025-02-08_mutect_dilutions_bwa.sh && \
bash 2025-02-08_mutect_all_samples_bwa.sh

echo "All Mutect2 scripts executed successfully."
