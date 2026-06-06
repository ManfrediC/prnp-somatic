# SPEC 07: LoB/LoD recalculation, plots, panels, and report

## Purpose

Recalculate LoB/LoD for every v2 method using one common schema, then generate
reviewable plots, panels, and a PDF report.

## Evidence Basis

- The current official ddPCR workflow defines the manuscript-facing LoB/LoD
  logic.
- dMIQE 2020 requires transparent reporting of analysis choices, controls, and
  software.
- The 2023 Clinical Chemistry review highlights that classification differences
  can materially change digital PCR estimates, so method comparisons must show
  both droplet geometry and quantitative calls.

## Inputs

- `tables/*_well_counts.csv`
- `data/droplets/*.rds`
- `models/**`
- `raw/ddpcr/sample_details.xlsx`
- official ddPCR helper functions for Poisson and fractional abundance logic.

## Outputs

- `tables/method_well_counts.csv`
- `tables/method_sample_region_results.csv`
- `tables/method_lob_lod_pass_matrix.csv`
- `tables/method_positive_rows.csv`
- `tables/method_summary.csv`
- `tables/method_control_false_positive_summary.csv`
- `tables/method_positive_control_recovery_summary.csv`
- `plots/individual/**/*.svg`
- `plots/individual/**/*.pdf`
- `plots/panels/*.svg`
- `plots/panels/*.pdf`
- `report/ddpcr_gating_method_comparison_v2.pdf`

## LoB/LoD Rules

- Aggregate to the same sample-region unit as the official workflow.
- Use the official Poisson and fractional abundance calculation.
- Use WT controls as blanks.
- Prefer same-plate/run blanks with assay-wide fallback.
- Apply existing assay-specific LoD fractional-abundance thresholds:
  - D178N: 0.056;
  - E200K: 0.067;
  - P102L: 0.13.
- Exclude known germline E200K rows from somatic interpretation but keep them
  visible as classifier checks.
- For probabilistic methods, use expected counts and label the result
  expected-count based.

## Plot Rules

- Individual plots are ggplot SVG/PDF.
- Panel assembly is done by Python from SVG/PDF assets.
- Axis limits are consistent within assay and panel family.
- Numerical counts use all droplets.
- Plot downsampling is allowed only after classification.
- Rain and unclassified droplets are shown explicitly.

## Report Sections

1. Scope and data completeness.
2. Literature and package manual evidence.
3. JSON conversion audit.
4. Control geometry and gate calibration.
5. Method-by-method implementation details.
6. Control validation.
7. LoB/LoD results.
8. Figures and panels.
9. Failure modes and limitations.
10. Recommendation for official workflow.

## E2E Checks

- Every method listed in specs has a row in `method_summary.csv` or a documented
  failure row.
- Positive rows table is generated.
- Control false-positive and positive-control recovery tables are generated.
- PDF report renders without missing figures.
- SVG and PDF panels exist and are non-empty.
- Report states package versions and parameters.

## Failure Handling

Do not report a method as successful if only a fallback or adapter ran. Label
native, control-anchored, one-channel, probabilistic, and failed variants
separately.

