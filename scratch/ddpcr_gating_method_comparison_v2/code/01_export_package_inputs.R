source("scratch/ddpcr_gating_method_comparison_v2/code/lib_v2.R")

setup_v2_dirs()

well_manifest <- build_well_manifest()
raw_audit <- audit_raw_exports(well_manifest)

parsed_path <- path_v2("data", "parsed_wells.rds")
shared_path <- path_v2("data", "shared_droplets.rds")

if (file.exists(parsed_path)) {
  parsed_wells <- readRDS(parsed_path)
} else {
  parsed_wells <- parse_complete_wells(raw_audit)
  saveRDS(parsed_wells, parsed_path)
}

if (file.exists(shared_path)) {
  shared_droplets <- readRDS(shared_path)
} else {
  shared_droplets <- shared_droplet_table(parsed_wells)
  saveRDS(shared_droplets, shared_path)
}

shared_thresholds <- shared_threshold_table(parsed_wells)
write_csv(shared_thresholds, path_v2("data", "shared_thresholds.csv"))

twoddpcr_manifest <- export_twoddpcr_inputs(parsed_wells)
ddpcrclust_manifest <- export_ddpcrclust_inputs(parsed_wells)
dpcp_manifest <- export_dpcp_inputs(parsed_wells)
one_channel_manifest <- export_one_channel_inputs(parsed_wells)

package_manifest <- bind_rows(
  twoddpcr_manifest,
  ddpcrclust_manifest,
  dpcp_manifest,
  one_channel_manifest
)

write_csv(package_manifest, path_v2("tables", "package_input_manifest.csv"))

validation <- validate_package_inputs(package_manifest)
write_csv(validation, path_v2("tables", "package_input_validation.csv"))

cat("parsed_wells=", length(parsed_wells), "\n", sep = "")
cat("shared_droplet_rows=", nrow(shared_droplets), "\n", sep = "")
cat("shared_threshold_rows=", nrow(shared_thresholds), "\n", sep = "")
cat("package_manifest_rows=", nrow(package_manifest), "\n", sep = "")
print(validation)
