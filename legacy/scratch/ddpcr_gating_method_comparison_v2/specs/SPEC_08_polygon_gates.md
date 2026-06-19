# SPEC 08: control-anchored polygon gates

## Purpose

Add an explicit polygon-gating comparator to the v2 ddPCR method comparison.
This covers the manual/free-form cluster-gate family that is distinct from
axis-aligned thresholds, k-means, kNN, Mahalanobis rain filters, package-native
clustering, one-channel thresholding, and Bayesian mixture classification.

## Evidence Basis

- dMIQE 2020 defines thresholds, clusters, rain, accepted partitions, and
  requires transparent reporting of controls and analysis choices. It also
  states that positive and negative controls should be used for assay quality
  control and threshold setting.
- Bio-Rad QuantaSoft Analysis Pro documents 2D cluster mode in which square,
  circle, or free-form shapes are drawn around droplet clusters and then mapped
  to target combinations.
- The `ddpcr` package overview describes two-channel ddPCR as four amplitude
  clusters and notes that gating is a key analysis step, with manual gating
  commonly used when automatic gates are inadequate.

References:

- dMIQE 2020: https://academic.oup.com/clinchem/article/66/8/1012/5880117
- QuantaSoft Analysis Pro manual:
  https://www.bio-rad.com/webroot/web/pdf/lsr/literature/QuantaSoft-Analysis-Pro-v1.0-Manual.pdf
- `ddpcr` overview:
  https://rdrr.io/cran/ddpcr/f/inst/doc/overview.Rmd

## Inputs

- `data/shared_droplets.rds`
- `data/parsed_wells.rds`
- `models/control_geometry/control_geometry.rds`

## Gate Construction

For each assay and each biological droplet class (`NN`, `WT`, `MUT`, `DP`):

1. Use only NTC, WT-control, and mutant-positive-control droplets.
2. Start from the class labels exported in the raw Quantasoft JSONs for control
   wells only.
3. Remove peripheral control droplets using the already fitted control
   Mahalanobis gate for that class.
4. Build a convex-hull polygon from the remaining control droplets when at
   least three hull points are available.
5. Expand the polygon slightly around the control-derived centroid so the hull
   is not an exact memorisation of control droplets.
6. Fall back to the control-geometry ellipse sampled as a polygon when a class
   is too sparse or degenerate, which is expected to matter most for sparse
   double-positive control droplets.

The polygon method intentionally stays on the same raw-amplitude coordinate
surface as the other v2 methods. Baseline-shifted polygon gates are a possible
future sensitivity analysis, but they should not be introduced into only this
method during the cross-method comparison.

## Classification

For every droplet in every complete well:

- assign the droplet to a class if it lies inside exactly one class polygon;
- if it lies inside multiple polygons, assign the nearest control-geometry class
  by Mahalanobis distance;
- if it lies outside all class polygons, call it `Rain`;
- reserve `Unclassified` for non-finite amplitudes or invalid inputs.

Method ID:

- `polygon_control_hull_gates`

## Outputs

- `models/polygon_gates/polygon_vertices.csv`
- `data/droplets/polygon_control_hull_gates_plot_droplets.rds`
- `tables/polygon_parameter_grid.csv`
- `tables/polygon_run_status.csv`
- `tables/polygon_well_counts.csv`
- `tables/polygon_control_validation.csv`
- `tables/polygon_plot_manifest.csv`
- `tables/polygon_e2e_checks.csv`
- `plots/individual/polygon_gates/*.svg`
- `plots/individual/polygon_gates/*.pdf`

## E2E Checks

- Four class polygons exist for every assay.
- Every polygon has at least three vertices and finite coordinates.
- Every complete well has one `polygon_control_hull_gates` count row.
- All count rows are marked `ok`.
- Control-validation rows are non-empty.
- Representative plot-droplet rows are written.
- SVG/PDF polygon-gating plots exist and are non-empty.

## Failure Handling

If a class lacks enough usable control data, the fallback polygon must be
labelled in `polygon_vertices.csv`. If any assay cannot produce four finite
polygons after fallback, the method is a failed comparator and must not be
included as an apparently successful LoB/LoD method.
