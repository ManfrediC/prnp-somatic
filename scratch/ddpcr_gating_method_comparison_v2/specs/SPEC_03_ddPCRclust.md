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

- `models/ddPCRclust/native_results.rds`
- `models/ddPCRclust/low_level_results.rds`
- `data/droplets/ddPCRclust_*.rds`
- `tables/ddPCRclust_well_counts.csv`
- `tables/ddPCRclust_reader_validation.csv`
- `tables/ddPCRclust_parameter_grid.csv`
- `plots/individual/ddPCRclust/*.svg`

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

## Low-Level Algorithm Grid

For each assay/run and selected parameter grid:

- `runPeaks()`;
- `runDensity()`;
- `runSam()`.

Run modes:

- native per-well;
- control-augmented, where controls provide the expected cluster context and
  final counts are computed only from sample droplets.

Method IDs:

- `ddPCRclust_runPeaks_native`
- `ddPCRclust_runDensity_native`
- `ddPCRclust_runSam_native`
- `ddPCRclust_runPeaks_control_augmented`
- `ddPCRclust_runDensity_control_augmented`
- `ddPCRclust_runSam_control_augmented`

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
- Native `ddPCRclust()` completes on at least one assay/run smoke test.
- All package errors are logged with input file paths.
- The full run produces common-schema well counts.
- Control false-positive and positive-control recovery summaries are written.

## Failure Handling

Native package failures are not replaced silently. They are logged, then
control-augmented or control-projected variants are run as separate methods.

