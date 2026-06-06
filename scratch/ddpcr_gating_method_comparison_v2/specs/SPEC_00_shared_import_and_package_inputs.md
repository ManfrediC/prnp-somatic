# SPEC 00: shared JSON import and package input conversion

## Purpose

Convert the immutable QuantaSoft JSON exports in `raw/ddpcr` into a shared
droplet table and package-specific input files without changing raw data or the
official `src/ddpcr` workflow.

This stage is the prerequisite for every native package workflow. A package is
not considered properly tried unless it receives input in the format expected by
its manual or vignette.

## Evidence Basis

- dMIQE 2020 requires transparent reporting of digital PCR pre-processing,
  partition classification, controls, and software settings.
- `ddPCRclust` vignette expects two-column amplitude CSVs and a run template.
- `dPCP` vignette requires one amplitude file per sample/reference plus a
  nine-column sample table.
- `dpcR::ddpcRquant()` manual expects a directory of raw data.
- The repository's official raw-import helpers already know how to select the
  WT and mutant targets from exported metadata.

## Inputs

- `raw/ddpcr/manifests/runs.csv`
- `raw/ddpcr/manifests/sample_manifest.csv`
- `raw/ddpcr/**/PeakData/<well>.ddpeakjson`
- `raw/ddpcr/**/PeakMetaData/<well>.ddmetajson`
- `src/ddpcr/ddpcr_raw_import_helpers.R`

## Outputs

- `tables/well_manifest.csv`
- `tables/raw_export_audit.csv`
- `tables/missing_raw_amplitude_or_metadata_wells.csv`
- `data/parsed_wells.rds`
- `data/shared_droplets.rds`
- `data/shared_thresholds.csv`
- `inputs/twoddpcr/**`
- `inputs/ddPCRclust/**`
- `inputs/dPCP/**`
- `inputs/definetherain/**`
- `inputs/ddpcRquant/**`
- `tables/package_input_manifest.csv`
- `tables/package_input_validation.csv`

## Data Contract

Every imported droplet must have:

```text
droplet_id
run_id
run_date
assay
well
sample
sample_type
sample_key
droplet_index
ch1_amplitude
ch2_amplitude
Ch1.Amplitude
Ch2.Amplitude
quantaSoft_class
target_class
ch1_role
ch2_role
ch1_target
ch2_target
source_peak_data_json
source_peak_metadata_json
```

`droplet_id` must be reversible and unique:

```text
run_id::assay::well::droplet_index
```

## Package-Specific Conversion Rules

### twoddpcr

Store per-assay and per-run RDS objects:

- a droplet data frame with `Ch1.Amplitude`, `Ch2.Amplitude`;
- a manifest with run, well, assay, sample, and control type;
- a plate-ready object if `ddpcrPlate` construction is reliable.

### ddPCRclust

Write one amplitude CSV per well:

```text
Ch1 Amplitude,Ch2 Amplitude
```

Write one run template per assay/run:

```text
> <run_id>, channel1=<name>, channel2=<name>
Well,Sample type,# of markers,Marker 1,Marker 2,Marker 3,Marker 4
```

The template must pass `ddPCRclust::readTemplate()` and the files must pass
`ddPCRclust::readFiles()`.

### dPCP

Write per-assay/per-run folders to avoid repeated Bio-Rad well IDs across runs.

Amplitude CSVs must initially match the vignette examples:

```text
Ch1 Amplitude,Ch2 Amplitude,Cluster
```

The `Cluster` column may contain imported QuantaSoft metadata labels for
reference compatibility, but native dPCP calls remain the package output, not
the metadata labels.

Sample tables must have exactly the nine dPCP columns:

```text
Sample name
Chip ID/Well ID
No of targets
FAM target
Target 3
Target 4
VIC/HEX target
Reference
Dilution
```

### definetherain

Write single-channel amplitude files and a control manifest. The method is a
one-channel threshold/rain method, so channel combination happens after native
or reimplemented one-channel calls.

### ddpcRquant

Write channel-specific raw-data directories. Each directory contains
two-column Bio-Rad-style amplitude CSVs named `<well>_Amplitude.csv` and one
`summary.csv` with `Well`, `Sample`, `TypeAssay`, and `Assay`. Validate the
exact layout with a small `dpcR::ddpcRquant()` run before full E2E execution.

## E2E Checks

- Every complete manifest well has either a parsed droplet object or an explicit
  missing raw-export row.
- Exported package input row counts equal source JSON droplet counts.
- `ddPCRclust::readFiles()` and `readTemplate()` can read one generated run.
- `dPCP::read_sampleTable()` can read one generated sample table.
- `dPCP::read_sample()` and `read_reference()` can read one generated job.
- `dpcR::ddpcRquant()` can run on a generated control-containing ddpcRquant
  directory, or the failure is logged with the tried layout.

## Failure Handling

Missing `PeakData` amplitude JSONs are not imputed. They are excluded from
amplitude-based classification and reported in
`missing_raw_amplitude_or_metadata_wells.csv`.

Any package-format conversion failure must identify:

- source JSON path;
- intended package;
- generated target path;
- row count expected;
- row count observed;
- package reader error, if any.
