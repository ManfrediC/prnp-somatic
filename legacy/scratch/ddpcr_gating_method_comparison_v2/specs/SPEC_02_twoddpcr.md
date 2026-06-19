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
- `data/droplets/twoddpcr_plot_droplets.rds`
- `tables/twoddpcr_well_counts.csv`
- `tables/twoddpcr_parameter_grid.csv`
- `tables/twoddpcr_control_validation.csv`
- `plots/individual/twoddpcr/*.svg`

Full classifications are computed for every droplet before aggregation. To
avoid duplicating the 12.5M-row shared droplet cache for every variant, the
persisted droplet-level artefact is a deterministic downsample used for plots
and audit diagnostics. Full per-droplet classifications remain reproducible
from `data/shared_droplets.rds`, `models/control_geometry/**`, and the method
script.

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
- `twoddpcr_knn_k3_mah_rain`
- `twoddpcr_knn_k5_mah_rain`
- `twoddpcr_knn_k11_mah_rain`
- `twoddpcr_knn_k21_mah_rain`

## Parameter Selection

Use leave-control-well-out validation:

- NTCs should remain negative.
- WT controls define false-positive behaviour.
- Positive controls must retain MUT and DP recovery.
- Rain fraction should be plausible and not collapse true positive controls.

## E2E Checks

- Each variant completes for all complete wells or logs a native package error.
- Native k-means failures caused by empty clusters are retained in the status
  table rather than silently replaced by control-centred results.
- kNN training rows, class balance, and `k` values are recorded.
- `mahalanobisRain()` distances are recorded by class.
- Counts are written to the common method schema.
- Control validation is generated before LoB/LoD application.
- Posterior/Bayesian classifiers are not implemented here; they remain reserved
  for `SPEC_06`.

## Failure Handling

If `kmeansClassify()` returns empty-cluster errors, the native variant logs
failure. A control-centred or kNN variant may still proceed, but it must not be
reported as native k-means success.
