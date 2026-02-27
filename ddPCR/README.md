# ddPCR raw data

This directory contains raw CSV files generated during the ddPCR experiments in this project (assay runs for `D178N`, `E200K`, and `P102L`).

## Purpose

These files are the original experimental outputs and are treated as raw input data for downstream processing and analysis. Raw CSVs are ignored from version control via [`ddPCR/.gitignore`].

## Sample metadata file

`sample_details.xlsx` contains the ddPCR sample metadata used by the reproducible ddPCR script (`src/ddPCR/create_snv_dataframe.R`) to map sample/well identifiers to participant-level information during dataframe construction.
