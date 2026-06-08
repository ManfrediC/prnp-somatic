# ddPCR LoD Raw Data

This directory contains immutable source material for the ddPCR limit-of-
detection (LoD) experiments used for the `D178N`, `E200K`, and `P102L`
assays.

The raw payload is ignored from version control except for this README. Build
or refresh the local copy with:

```bash
python src/tools/build_ddpcr_lod_raw_database.py
```

## Expected Local Layout

The database is grouped by assay:

```text
raw/ddpcr_lod/
  d178n/
  e200k/
  p102l/
```

Each assay directory contains:

- `ddpcr_archive/`: copied `.ddpcr` instrument archive with a canonical name.
- `archive_contents/`: exported archive contents, including droplet-amplitude
  JSONs and related QuantaSoft files.
- `csv_export/`: raw or full QuantaSoft CSV export for the source run.
- `analysis_csv/`: cleaned or merged LoD CSV used by the manuscript analysis.
- `layout_xlsx/`: source layout workbook for the LoD plate.
- `manifests/`: local provenance, file, sample, archive-content, and excluded
  material manifests.

## Provenance Policy

Files are renamed to canonical, identifier-free names such as
`2020-09-11_SNV_D178N.ddpcr`. The untracked manifests link every local file to
its original source path under the USZ ddPCR experiment directory and record
checksums for copied files.
