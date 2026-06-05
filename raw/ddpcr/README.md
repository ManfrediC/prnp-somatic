# ddPCR Raw Data

This directory contains immutable ddPCR source exports for the `D178N`,
`E200K`, and `P102L` assays.

## Purpose

These files are raw input data for downstream processing and analysis. The
directory is ignored from version control except for this README.

Expected local layout:

- `ddpcr_archive/`: copied `.ddpcr` instrument archives.
- `archive_contents/`: exported droplet-amplitude JSONs and related archive contents.
- `csv_export/`: canonical CSV exports used to validate the raw JSON import.
- `layout_xlsx/`: source layout workbooks used for run/sample mapping.
- `manifests/`: provenance, run, sample, layout, supersession, and QC manifests.

## Sample metadata file

`sample_details.xlsx` contains the ddPCR sample metadata used by
`src/ddpcr/create_snv_dataframe.R` to map sample and well identifiers to
participant-level information during dataframe construction.

## Export Tooling

The local `.ddpcr` archive export/decryption helper is tooling, not raw data.
Keep it untracked at `src/tools/export_ddpcr_archives.py`.
