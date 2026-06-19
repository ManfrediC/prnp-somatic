# SPEC 03: ddPCRclust

## Purpose

Run `ddPCRclust` through its native template workflow and low-level clustering
functions, then add a clearly labelled control-projected sensitivity variant.

## Evidence Basis

- `ddPCRclust` vignette documents `readFiles()`, `readTemplate()`, and
  `ddPCRclust(files, template)`.
- The package expects two-column amplitude CSVs and a run template.
- The manual for `runPeaks()` says it finds rain and assigns rain using
  distances to vector lines connecting cluster centres.
- `runPeaks()` exposes `missingClusters`, `similarityParam`, and
  `distanceParam`, which are relevant when sparse clusters are expected.

## Inputs

- `inputs/ddPCRclust/**`
- `models/control_geometry/**`
- `data/parsed_wells.rds`

## Outputs

- `data/droplets/ddPCRclust_plot_droplets.rds`
- `tables/ddPCRclust_well_counts.csv`
- `tables/ddPCRclust_reader_validation.csv`
- `tables/ddPCRclust_parameter_grid.csv`
- `tables/ddPCRclust_run_status.csv`
- `tables/ddPCRclust_low_level_smoke.csv`
- `plots/individual/ddPCRclust/*.svg`

Full package classifications are aggregated to well counts. Persisted
droplet-level output is a deterministic diagnostic sample; full package result
objects remain local/regenerable because native ensemble outputs are bulky.

## Native Template Workflow

For each assay/run:

```r
files <- readFiles(amplitude_files)
template <- readTemplate(template_file)
result <- ddPCRclust(files, template, numOfMarkers = 2)
```

Then:

1. Extract cluster labels and counts.
2. Map package clusters to NN, WT, MUT, DP using known controls and expected
   channel geometry.
3. Convert counts into the common well schema.

Method ID:

- `ddPCRclust_template_native`
- `ddPCRclust_template_fast`

## Low-Level Algorithm Smoke

The exported low-level functions are audited on representative NTC, WT-control,
and positive-control wells from each assay:

- `runPeaks()`;
- `runDensity()`;
- `runSam()`.

The full template workflow already invokes density, SAM, and peaks internally
when `fast = FALSE`; the explicit low-level smoke records which functions fail
on sparse rare-variant wells and avoids duplicating the native ensemble over
hundreds of wells.

These low-level smoke rows are not LoB/LoD candidate methods unless they
produce complete well-count tables in a later run.

## Control-Projected Sensitivity Variant

Use ddPCRclust on high-quality controls to infer cluster geometry, then project
sample droplets onto fixed control geometry using the shared Mahalanobis
acceptance rules.

Method ID:

- `ddPCRclust_control_projected`

This is not a native ddPCRclust workflow and must be labelled as such.

## E2E Checks

- Generated files pass `readFiles()`.
- Generated templates pass `readTemplate()`.
- Native `ddPCRclust()` completes or logs package failures for each assay/run
  and mode.
- All package errors are logged with input file paths.
- The full run produces common-schema well counts.
- Control false-positive and positive-control recovery summaries are written.
- Low-level smoke rows are written for `runDensity()`, `runPeaks()`, and
  `runSam()`.

## Failure Handling

Native package failures are not replaced silently. They are logged, then
control-augmented or control-projected variants are run as separate methods.
