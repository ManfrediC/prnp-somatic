# SPEC 02: twoddpcr

## Purpose

Run `twoddpcr` using the package's native classifiers and a control-anchored
rare-variant training strategy:

- k-means;
- k-nearest-neighbour classification;
- Mahalanobis rain handling.

## Evidence Basis

- `twoddpcr` manual and vignette document `kmeansClassify()`,
  `knnClassify()`, `mahalanobisRain()`, and `sdRain()`.
- The vignette states that classifying wells together helps avoid judging
  whether clusters in individual wells are empty.
- The vignette uses k-means to clean a low-noise training set, then
  `mahalanobisRain()`, removes rain, and uses `knnClassify()`.
- The `twoddpcr` paper describes two-channel four-class ddPCR classification.

## Inputs

- `inputs/twoddpcr/**`
- `models/control_geometry/**`
- `data/parsed_wells.rds`

## Outputs

- `models/twoddpcr/training_sets.rds`
- `data/droplets/twoddpcr_*.rds`
- `tables/twoddpcr_well_counts.csv`
- `tables/twoddpcr_parameter_grid.csv`
- `tables/twoddpcr_control_validation.csv`
- `plots/individual/twoddpcr/*.svg`

## Variants

### Native k-means

Run `kmeansClassify()`:

- default package centres;
- control-derived initial centres;
- per-assay/run plate-level classification where feasible.

Method IDs:

- `twoddpcr_kmeans_native`
- `twoddpcr_kmeans_control_centres`

### kNN

Build control-only training sets:

1. Use NTC, WT, positive controls, and low-VAF spike-ins if available.
2. Class labels may come from QuantaSoft metadata in controls or a cleaned
   control-only k-means plus Mahalanobis-rain pass.
3. Remove rain from training.
4. Downsample NN/WT if needed so rare MUT/DP neighbours are not overwhelmed.
5. Run `knnClassify()` for `k = 3, 5, 11, 21`.

Method IDs:

- `twoddpcr_knn_k3`
- `twoddpcr_knn_k5`
- `twoddpcr_knn_k11`
- `twoddpcr_knn_k21`

### Mahalanobis Rain

Apply `mahalanobisRain()` to k-means and kNN classifications. Tune
`maxDistances` from controls only.

Method IDs:

- `twoddpcr_kmeans_mah_rain`
- `twoddpcr_knn_mah_rain`

## Parameter Selection

Use leave-control-well-out validation:

- NTCs should remain negative.
- WT controls define false-positive behaviour.
- Positive controls must retain MUT and DP recovery.
- Rain fraction should be plausible and not collapse true positive controls.

## E2E Checks

- Each variant completes for all complete wells or logs a native package error.
- kNN training rows, class balance, and `k` values are recorded.
- `mahalanobisRain()` radii are recorded by class.
- Counts are written to the common method schema.
- Control validation is generated before LoB/LoD application.

## Failure Handling

If `kmeansClassify()` returns empty-cluster errors, the native variant logs
failure. A control-centred or kNN variant may still proceed, but it must not be
reported as native k-means success.

