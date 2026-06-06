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
