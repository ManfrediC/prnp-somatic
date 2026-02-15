#!/bin/bash

# Run the first script, show output on-screen and save it
/add/code/mutect_with_pon/mutect_pon_pipeline.sh \
  2>&1 | tee /add/code/mutect_with_pon/mutect_pon_pipeline.log
if [ ${PIPESTATUS[0]} -eq 0 ]; then
  echo "First script completed successfully."

  # Run the second script, show output on-screen and save it
  /add/code/mutect_with_pon/readcounts_QC_pipeline_CJD.sh \
    2>&1 | tee /add/code/mutect_with_pon/readcounts_QC_pipeline_CJD.log
  if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo "Second script completed successfully."
  else
    echo "Second script failed."
  fi
else
  echo "First script failed."
fi
