# SPEC 05: definetherain and ddpcRquant

## Purpose

Run `definetherain` and `ddpcRquant` as one-channel rain/threshold methods,
then combine channel states into two-channel SNV partitions with explicit
limitations.

## Evidence Basis

- `ddpcRquant` manual describes single-channel threshold determination using a
  raw-data directory and default `threshold.int = 0.9995`.
- The `dpcR` manual states `ddpcRquant()` is based on the ddpcRquant web code.
- The definetherain method uses positive-control-based one-channel cluster/rain
  thresholds.
- The 2023 Clinical Chemistry review warns that rain handling can bias digital
  PCR results and discusses methods that model or threshold rain.

## Inputs

- `inputs/definetherain/**`
- `inputs/ddpcRquant/**`
- `data/parsed_wells.rds`

## Outputs

- `models/definetherain/channel_models.csv`
- `models/ddpcRquant/channel_models.rds`
- `data/droplets/definetherain_channel_combined.rds`
- `data/droplets/ddpcRquant_*.rds`
- `tables/one_channel_well_counts.csv`
- `tables/one_channel_thresholds.csv`
- `plots/individual/one_channel/*.svg`

## definetherain

Implement the published method locally if no maintained callable package is
available:

1. Fit two k-means clusters in the positive-control channel amplitudes.
2. Lower cluster is negative; higher cluster is positive.
3. Negative valid region: below `negative_mean + 3 * negative_sd`.
4. Positive valid region: above `positive_mean - 3 * positive_sd`.
5. Intermediate region is rain.

Run separately for WT and mutant channels, then combine:

- WT negative, MUT negative -> NN;
- WT positive, MUT negative -> WT;
- WT negative, MUT positive -> MUT;
- WT positive, MUT positive -> DP;
- either channel rain -> Rain.

Method ID:

- `definetherain_channel_combined`

## ddpcRquant

Run `dpcR::ddpcRquant()` on channel-specific generated raw-data directories.

Test:

- `threshold.int = 0.995`;
- `threshold.int = 0.999`;
- `threshold.int = 0.9995`;
- `reps = 10`;
- `reps = 100` if runtime allows.

Combine channel calls as for definetherain.

Method IDs:

- `ddpcRquant_0995`
- `ddpcRquant_0999`
- `ddpcRquant_09995`

## E2E Checks

- A generated ddpcRquant directory runs through `dpcR::ddpcRquant()` or logs the
  exact accepted-layout failure.
- definetherain control thresholds are finite for each assay/channel.
- Rain-band overlaps are reported.
- Channel thresholds are not interpreted as two-dimensional clusters.
- Common-schema well counts are produced.

## Failure Handling

If a one-channel method cannot produce a finite threshold for a channel, the
affected wells are labelled threshold unavailable. Do not impute thresholds from
biological sample positivity.
