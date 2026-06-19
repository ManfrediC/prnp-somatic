# SPEC 06: control-anchored Bayesian mixture comparator

## Purpose

Implement a transparent probabilistic comparator for rare-variant ddPCR that
uses controls to fit class densities and computes posterior-weighted mutant and
WT counts.

This is not a named package workflow. It is a scientific comparator motivated by
the same issues addressed by Umbrella-like and mixture-model partition
classification methods.

## Evidence Basis

- The 2023 Clinical Chemistry review discusses partition-level probabilities,
  baseline correction, rain, and bias from excluding intermediate droplets.
- dMIQE 2020 requires transparent reporting of digital PCR analysis decisions.
- Control-anchored modelling is necessary because rare sample wells may not
  contain enough mutant-positive or double-positive droplets to train a model
  directly.

## Inputs

- `models/control_geometry/**`
- `data/parsed_wells.rds`
- `data/shared_droplets.rds`

## Outputs

- `models/bayesian_mixture/class_parameters.rds`
- `models/bayesian_mixture/priors.csv`
- `data/droplets/bayesian_mixture_plot_droplets.rds`
- `tables/bayesian_well_counts.csv`
- `tables/bayesian_uncertainty.csv`
- `tables/bayesian_e2e_checks.csv`
- `plots/individual/bayesian_mixture/*.svg`

## Model

For each droplet `x_i`:

```text
P(z_i = k | x_i) = pi_k f_k(x_i) / sum_j pi_j f_j(x_i)
```

Classes:

- NN;
- WT;
- MUT;
- DP;
- Rain.

Start with robust Gaussian components for NN, WT, MUT, and DP from
`SPEC_01`. Add a broad explicit Gaussian rain component centred between the
four biological class centres, with inflated covariance. If Gaussian tails
behave poorly, this comparator remains a sensitivity analysis rather than a
candidate final workflow.

## Control Fitting

- Fit centres and covariance from controls only.
- Use baseline-aligned amplitudes if control validation supports it.
- Use shrinkage covariance from `SPEC_01`.
- Priors are estimated from NTC and WT-control class assignments with additive
  smoothing, not from biological sample positivity or mutant-positive controls.
- MUT and DP priors are never zero, so sparse rare-positive droplets can still
  receive posterior mass when the likelihood supports them.

## Counting

Weighted counts:

```text
N_mut = sum(P(MUT) + P(DP))
N_wt  = sum(P(WT) + P(DP))
```

Hard-map counts:

```text
class = argmax_k P(z_i = k | x_i)
```

Both must be reported.

Method IDs:

- `bayesian_control_mixture_weighted`;
- `bayesian_control_mixture_hard_map`.

## E2E Checks

- Posterior rows sum to 1 within numerical tolerance.
- All class likelihoods are finite or explicitly replaced by a safe floor.
- WT controls preserve the false-positive design.
- Positive controls recover MUT and DP signal.
- Weighted counts and hard counts are both exported to common schema.

## Failure Handling

If posterior weights concentrate unrealistically in MUT or DP for WT controls,
the model is not a candidate workflow. It remains a failed sensitivity analysis
with diagnostic plots.
