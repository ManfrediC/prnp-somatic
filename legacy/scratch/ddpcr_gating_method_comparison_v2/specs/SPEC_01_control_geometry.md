# SPEC 01: shared control-anchored geometry

## Purpose

Fit a control-derived four-partition geometry for rare-variant ddPCR. This is
used by control-anchored method variants and by the Bayesian comparator.

The goal is to avoid training classifiers from rare biological samples alone,
where mutant-only and double-positive partitions may be too sparse to discover.

## Evidence Basis

- dMIQE 2020 requires appropriate controls and transparent analysis rules.
- The 2023 Clinical Chemistry review identifies partition classification,
  baseline shifts, and rain handling as major sources of digital PCR bias.
- `twoddpcr` vignette recommends classifying wells together when individual
  wells have empty clusters and uses Mahalanobis distance for rain.
- `dPCP` vignette recommends high-quality references when few positive
  partitions or rain are expected.

## Inputs

- `data/parsed_wells.rds`
- `data/shared_droplets.rds`
- `tables/well_manifest.csv`

## Outputs

- `models/control_geometry/centroids.csv`
- `models/control_geometry/covariances.rds`
- `models/control_geometry/gate_radii.csv`
- `models/control_geometry/baseline_shifts.csv`
- `models/control_geometry/control_training_manifest.csv`
- `plots/individual/control_geometry/*.svg`

## Model

For each assay, and where possible for each assay/run:

- NN: double-negative cloud.
- WT: reference-only / WT-positive cloud.
- MUT: mutant-only cloud.
- DP: double-positive cloud.

Centres are estimated from controls only. Biological sample droplets must not be
used to tune centres, covariance matrices, rain thresholds, or priors.

## Baseline Alignment

Estimate the NN centre per well. Align sample wells to the matched control NN
centre by a translation:

```text
aligned_x = raw_x - sample_NN + reference_NN
```

The raw coordinates remain stored. Classification methods may use aligned
coordinates only when the method spec explicitly allows baseline correction.

## DP Centre Rule

Use this hierarchy:

1. Observed positive-control DP centre if enough DP droplets exist.
2. Vector-sum prediction:

```text
mu_DP = mu_WT + mu_MUT - mu_NN
```

3. If both fail, mark DP geometry unavailable for that assay/run and do not
force DP calls.

## Covariance Rule

- Use robust/trimmed covariance for abundant classes.
- Use pooled covariance for sparse classes.
- Shrink covariance toward a diagonal pooled estimate.
- Add a small positive diagonal ridge to prevent singular matrices.

## Gate Radius Rule

Estimate Mahalanobis acceptance radii from controls:

- empirical radii when class control droplets are abundant;
- chi-square fallback when empirical estimates are weak;
- final radii must preserve WT-control false-positive behaviour.

## Droplet Assignment Rule

For hard assignments:

1. Compute Mahalanobis distance to NN, WT, MUT, and DP.
2. Pick the nearest class.
3. Accept only if distance is within that class's calibrated radius.
4. MUT and DP assignments must also satisfy channel-compatible signal evidence.
5. Otherwise assign `Rain`.

Never assign a droplet to DP solely because DP is the nearest centroid.

## E2E Checks

- Every assay has NN and WT geometry.
- MUT and DP geometry source is recorded as observed, vector-sum, or unavailable.
- Covariance matrices are finite and positive-definite after shrinkage/ridge.
- WT-control false positives are summarised before sample application.
- Positive-control recovery is summarised before sample application.
- Plots show controls, centres, and ellipses with consistent assay scales.

## Failure Handling

If a control type is missing or too sparse, do not silently borrow from samples.
Record the fallback level:

- assay/run matched;
- assay-wide;
- nearest assay-matched run;
- unavailable.

