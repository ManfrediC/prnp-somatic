# Experiments

## E1: v2 scaffold

Status: completed

Hypothesis: A separate v2 scratch tree will allow native package workflows and
control-anchored variants to be implemented without disturbing the completed v1
comparison or official workflow files.

Mechanism: Created the `scratch/ddpcr_gating_method_comparison_v2` directory
tree with code, inputs, models, data, tables, plots, report, logs, and
literature folders.

Decision rule: Keep if all generated files remain under ignored `/scratch`.

Result: Scaffold created.

Next: Implement shared JSON import and package input converters.

## E2: shared JSON import and package input conversion

Status: completed

Hypothesis: The exported QuantaSoft JSONs can be converted into shared droplet
tables and method-specific package input files that pass package reader smoke
tests.

Mechanism: Added `specs/SPEC_00_shared_import_and_package_inputs.md`,
`code/lib_v2.R`, `code/00_audit_controls.R`, and
`code/01_export_package_inputs.R`. The converter parses complete active SNV
wells, writes shared droplet/threshold artefacts, exports package inputs for
`twoddpcr`, `ddPCRclust`, `dPCP`, `definetherain`, and `ddpcRquant`, and runs
reader/package smoke validation.

Decision rule: Keep if all complete wells parse, missing amplitude exports are
logged, generated inputs remain under scratch, and `ddPCRclust`, `dPCP`, and
`ddpcRquant` smoke validation pass.

Result: Parsed 824 complete wells and 12,475,437 droplets. Twenty-one D178N
wells lack exported `PeakData` amplitude JSONs and remain logged in
`tables/missing_raw_amplitude_or_metadata_wells.csv`. Generated 1,843 package
input manifest rows. `ddPCRclust::readFiles()`/`readTemplate()`,
`dPCP::read_sampleTable()`, and a small `dpcR::ddpcRquant()` run all pass. The
`ddpcRquant()` smoke call requires `threshold.manual = FALSE` because the
package's documented `NULL` default reaches a length-zero internal `if`.

Next: Fit shared control-anchored geometry.

## E3: control-anchored geometry

Status: completed

Hypothesis: Control wells can provide assay-specific centroids, covariance
matrices, empirical gate radii, and baseline-shift summaries that downstream
classifiers can use without training rare-positive classes from sparse
biological samples alone.

Mechanism: Added `specs/SPEC_01_control_geometry.md` and
`code/02_fit_control_geometry.R`. The script fits NN, WT, MUT, and DP centres
from NTC, WT-control, and positive-control wells; regularises covariance
matrices; derives Mahalanobis acceptance radii; records baseline shifts; writes
control validation tables; and exports one SVG/PDF diagnostic pair per assay.

Decision rule: Keep if the script produces the expected three assays, four
classes per assay, finite centroids and covariance rows, positive gate radii,
non-empty baseline and validation tables, and all diagnostic plot files.

Result: Passed 9 built-in E2E checks. Generated 12 centroid rows, 12 covariance
rows, 12 gate-radius rows, 117 baseline-shift rows, 143 control-validation
rows, and SVG/PDF diagnostic plots for D178N, E200K, and P102L.

Next: Implement `twoddpcr` native and control-anchored variants.

## E4: `twoddpcr` k-means, kNN, and Mahalanobis rain

Status: completed

Hypothesis: `twoddpcr` can be run as both a native four-cluster comparator and
a control-anchored rare-variant comparator when kNN training is built from
control-geometry-cleaned NTC, WT-control, and positive-control droplets.

Mechanism: Updated `specs/SPEC_02_twoddpcr.md` and added
`code/03_run_twoddpcr.R`. The script runs native k-means, control-centred
k-means, native and control-centred Mahalanobis rain variants, kNN for
`k = 3, 5, 11, 21`, and kNN plus Mahalanobis rain for each k. The kNN training
set uses up to 100 balanced control droplets per class and assay after
control-geometry cleaning. Full classifications are aggregated to well-level
counts; a deterministic sampled droplet table is retained locally for plots.

Decision rule: Keep if all methods produce one well-count row for every
complete well, all native package failures are logged with messages, training
contains four classes per assay, control validation is non-empty, run-status
rows cover every method/run combination, and sampled plot droplets are written.

Result: Passed 7 built-in E2E checks. Generated 12 methods, 9,888 well-count
rows, 11,709 control-validation rows, 468 run-status rows, and 1,923,593 sampled
plot-droplet rows. Native package-default k-means failed on 4 of 39 assay/run
jobs because of empty clusters; these failures correspond to 63 well rows and
are retained as native failures. Control-centred k-means and all kNN variants
completed for all 824 wells.

Next: Implement `ddPCRclust` native and control-anchored variants.

## E5: `ddPCRclust` native template and control-projected variants

Status: completed

Hypothesis: The `ddPCRclust` template workflow can be run against the exported
amplitude CSV/template inputs, with package clusters mapped to assay biology by
control-derived geometry, while low-level function failures on sparse wells are
logged explicitly.

Mechanism: Updated `specs/SPEC_03_ddPCRclust.md` and added
`code/04_run_ddpcrclust.R`. The script runs `ddPCRclust()` in full native
ensemble mode and fast density mode for every assay/run, maps package cluster
centres to NN/WT/MUT/DP/Rain using control geometry, writes common well-count
and control-validation tables, and runs `runDensity()`, `runPeaks()`, and
`runSam()` smoke checks on representative NTC, WT-control, and positive-control
wells from each assay. It also writes a labelled
`ddPCRclust_control_projected` sensitivity comparator.

Decision rule: Keep if all generated files pass `readFiles()`/`readTemplate()`,
each candidate method has one well-count row per complete well, run-status rows
cover every method/run combination, package failures remain logged with
messages, low-level smoke rows are present, and plot-sample droplets are
written.

Result: Passed 8 built-in E2E checks. Generated 3 methods, 2,472 well-count
rows, 2,805 control-validation rows, 117 run-status rows, 27 low-level smoke
rows, and 311,739 sampled plot-droplet rows. `ddPCRclust_template_fast` had 97
failed well rows and `ddPCRclust_template_native` had 265 failed well rows after
package/result mapping; both retain failure messages. The control-projected
sensitivity completed all 824 wells. Low-level smoke succeeded for all 9
`runSam()` checks and for 8 of 9 `runDensity()`/`runPeaks()` checks.

Next: Implement the native `dPCP()` workflow and compare it against the v1
adapter-style result.

## E6: native `dPCP()` and control-projected adapter comparison

Status: completed

Hypothesis: Native `dPCP()` can be run end-to-end from the converted Bio-Rad
amplitude files and sample tables; if native clustering fails, the failures can
be compared against a control-projected adapter-style baseline on the same
inputs.

Mechanism: Updated `specs/SPEC_04_dPCP_native.md` and added
`code/05_run_dpcp.R`. The script validates `read_sampleTable()`,
`read_reference()`, and `read_sample()`, then runs native `dPCP()` with a fixed
DBSCAN retry grid: `(eps, minPts) = (200,50), (150,30), (250,50), (300,50),
(200,30)`. Native droplet labels are extracted from `final cluster` when
available and mapped to NN/WT/MUT/DP/Rain. The comparison baseline is
`dPCP_control_projected`, using the same converted dPCP inputs and shared
control geometry.

Decision rule: Keep if reader validation passes, every method has one
well-count row per complete well, native failures are logged with messages,
the adapter comparison has one row per well, and sampled plot droplets are
written.

Result: Passed 8 built-in E2E checks. Generated 2 methods, 1,648 well-count
rows, 1,142 control-validation rows, 226 run-status rows, 824 adapter-comparison
rows, and 137,356 sampled plot-droplet rows. Native `dPCP()` succeeded for 2 of
39 assay/run jobs and 26 of 824 wells. The other 37 jobs failed through all five
pre-specified DBSCAN attempts, mostly with LAPACK `dgecon()` numerical errors;
three attempts failed because the reference had fewer predicted clusters than
expected, and one failed because initial centres were not distinct. The
control-projected baseline completed all 824 wells.

Next: Implement one-channel methods, including `definetherain` and
`ddpcRquant`, then the Bayesian classifiers.

## E7: one-channel `definetherain` and `ddpcRquant`

Status: completed

Hypothesis: One-channel rain/threshold methods can be run as explicit
comparators by classifying WT and mutant channels separately, then combining
channel calls into NN/WT/MUT/DP/Rain with the limitation that channel-wise
thresholds are not two-dimensional clusters.

Mechanism: Added `code/06_run_one_channel_methods.R`. Because no maintained
`definetherain` package is installed in the conda environment, the published
positive-control-based one-channel rain-band logic is implemented locally:
two-cluster k-means per assay/channel, negative valid region below
`negative_mean + 3 * negative_sd`, positive valid region above
`positive_mean - 3 * positive_sd`, and intermediate droplets called rain.
`dpcR::ddpcRquant()` is run for each assay/run/channel with
`threshold.int = 0.995, 0.999, 0.9995`, `reps = 10`, and
`threshold.manual = FALSE`, then channel thresholds are combined into
two-channel calls.

Decision rule: Keep if all methods produce one well-count row per complete
well, definetherain thresholds are finite for each assay/channel, ddpcRquant
thresholds are finite for each channel/run/interval, control validation is
non-empty, and sampled plot droplets are written.

Result: Passed 7 built-in E2E checks. Generated 4 methods, 3,296 well-count
rows, 3,459 control-validation rows, 240 threshold/model rows, 234 ddpcRquant
channel-run status rows, and 400,365 sampled plot-droplet rows. All one-channel
methods completed for all 824 wells.

Next: Implement the Bayesian/probabilistic classifiers.

## E8: Bayesian control-anchored mixture classifiers

Status: completed

Hypothesis: A control-anchored probabilistic classifier can quantify rare
mutant signal without training on rare biological samples alone by using
control-derived class densities, conservative NTC/WT-derived priors, and
posterior-weighted counts.

Mechanism: Updated `specs/SPEC_06_bayesian_mixture.md` and added
`code/07_run_bayesian_mixture.R`. The model uses Gaussian NN, WT, MUT, and DP
components from control geometry, a broad explicit rain component, and priors
estimated from NTC and WT-control geometry assignments with additive smoothing.
It exports both probability-weighted counts
(`bayesian_control_mixture_weighted`) and hard maximum-a-posteriori counts
(`bayesian_control_mixture_hard_map`).

Decision rule: Keep if class parameters are finite, priors sum to one by assay,
posterior probabilities sum to one for every processed droplet, both Bayesian
methods produce one well-count row per complete well, control validation is
non-empty, and sampled plot droplets are written.

Result: Passed 7 built-in E2E checks. Generated 2 methods, 1,648 well-count
rows, 1,037 control-validation rows, 824 uncertainty rows, 39 run-status rows,
and 127,635 sampled plot-droplet rows. MUT/DP priors are non-zero but small
after NTC/WT-only prior fitting: D178N MUT 0.00128 and DP 0.000746; E200K MUT
0.000271 and DP 0.000638; P102L MUT 0.000284 and DP 0.000587.

Next: Implement explicit polygon gates, then recalculate LoB/LoD pass tables
across all implemented methods and build the comparison report.

## E9: control-anchored polygon-gate comparator

Status: completed

Hypothesis: A control-anchored polygon-gate comparator can represent the
manual/free-form gate family without training on sparse rare-positive samples
alone, and can be compared on the same full-droplet count schema as the package
and Bayesian methods.

Mechanism: Added `specs/SPEC_08_polygon_gates.md` and
`code/08_run_polygon_gates.R`. The method builds one convex-hull polygon per
assay/class from NTC, WT-control, and mutant-positive-control droplets after
Mahalanobis trimming against the shared control geometry. Sparse or degenerate
classes fall back to the control-geometry ellipse sampled as a polygon. Droplets
inside one polygon are assigned to that class; droplets inside overlapping
polygons are assigned to the nearest control-geometry class; droplets outside
all polygons are called rain.

Decision rule: Keep if four finite polygons exist per assay, every complete well
has one count row, all count rows are `ok`, control validation is non-empty,
plot-droplet rows are written, and polygon SVG/PDF overlay plots exist.

Result: Passed 7 built-in E2E checks. Generated 326 polygon vertices, 824
well-count rows, 1,079 control-validation rows, 18,108 sampled plot-droplet
rows, 39 run-status rows, and assay-level polygon overlay SVG/PDF plots. The
method ID is `polygon_control_hull_gates`.

Next: Rerun the LoB/LoD synthesis and final report with polygon gates included.

## E10: cross-method LoB/LoD synthesis and report artefacts

Status: completed

Hypothesis: The implemented classifiers can be compared on a shared
sample-region surface by pooling complete well rows per method, recalculating
fractional abundance with the official helper, applying WT-control-derived LoB
and assay LoD thresholds, and rendering vector plots plus a PDF report.

Mechanism: Updated `code/09_summarise_lob_lod_report.R` and renamed the SVG
panel builder to `code/10_build_report_panels.py`. The R synthesis combines the
current exported JSON target-class calls with twoddpcr, polygon gates,
ddPCRclust, dPCP, definetherain, ddpcRquant, and Bayesian result tables; pools
wells by method, assay, and sample-region; recalculates fractional abundance
and confidence intervals; and applies WT-control LoB with same-plate preference
and assay-wide fallback. It exports summary tables, package versions,
downsampled representative gating plots, individual SVG/PDF plots, a PDF panel,
and the final PDF comparison report. The Python panel builder inlines the
R-generated summary SVG plots into a standalone SVG panel with prefixed svglite
clip-path IDs.

Decision rule: Keep if the method summary contains every implemented method,
all requested method families are represented, report and panel artefacts exist
and are non-empty, all method-level E2E checks pass, and the SVG panel builder
passes its own input/output/id-prefix checks.

Result: Passed 53 upstream method E2E rows, 10 final report checks, and 5 SVG
panel E2E checks. Generated 20,600 well-count rows, 16,893 sample-region rows,
25 method summary rows, and 3,831 biological LoB/LoD-positive rows after
excluding germline E200K rows from the biological pass matrix. Method coverage is
`current` 1, `twoddpcr` 12, `polygon_gates` 1, `ddPCRclust` 3, `dPCP` 2,
`definetherain` 1, `ddpcRquant` 3, and `bayesian` 2. The report artefacts are
`report/ddpcr_gating_method_comparison_v2.pdf`,
`plots/panels/method_summary_panel.pdf`, and
`plots/panels/method_summary_panel.svg`.

Next: Run the final completion audit, stage the reviewable code/spec/table
outputs, and commit the polygon plus expanded synthesis milestone.
