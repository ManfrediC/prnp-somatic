source("scratch/ddpcr_gating_method_comparison_v2/code/lib_v2.R")

setup_v2_dirs()

well_manifest <- build_well_manifest()
write_csv(well_manifest, path_v2("tables", "well_manifest.csv"))

raw_audit <- audit_raw_exports(well_manifest)
write_csv(raw_audit, path_v2("tables", "raw_export_audit.csv"))

complete_manifest <- raw_audit %>%
  filter(has_peak_data, has_peak_metadata)

control_summary <- complete_manifest %>%
  count(assay, sample_type, name = "wells") %>%
  arrange(assay, sample_type)
write_csv(control_summary, path_v2("tables", "control_well_summary.csv"))

missing_raw <- raw_audit %>%
  filter(!has_peak_data | !has_peak_metadata)
write_csv(missing_raw, path_v2("tables", "missing_raw_amplitude_or_metadata_wells.csv"))

cat("well_rows=", nrow(well_manifest), "\n", sep = "")
cat("complete_wells=", nrow(complete_manifest), "\n", sep = "")
cat("missing_raw_wells=", nrow(missing_raw), "\n", sep = "")

