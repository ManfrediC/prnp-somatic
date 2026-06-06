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
