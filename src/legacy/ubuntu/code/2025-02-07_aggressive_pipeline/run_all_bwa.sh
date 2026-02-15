# Set working directory to where the scripts are located
cd /add/code/2025-02-07_aggressive_pipeline/ || { echo "Failed to change directory"; exit 1; }

# Run the scripts sequentially, stopping if one fails:
bash 2025-02-07_dilutions_bwa.sh && \
bash 2025-02-07_first_CJD_seq.sh && \
bash 2025-02-07_CJD_8_samples_bwa.sh && \
bash 2025-02-07_CJD_16_samples_bwa.sh

echo "All pipeline scripts have been executed."