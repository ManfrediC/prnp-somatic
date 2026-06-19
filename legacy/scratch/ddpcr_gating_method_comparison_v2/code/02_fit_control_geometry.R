source("scratch/ddpcr_gating_method_comparison_v2/code/lib_v2.R")

setup_v2_dirs()

shared_droplets <- load_shared_droplets()
geometry <- fit_control_geometry(shared_droplets)

geometry_checks <- tibble(
  check = c(
    "expected_assays",
    "four_classes_per_assay",
    "finite_centres",
    "finite_covariances",
    "positive_gate_radii",
    "baseline_shifts_present",
    "control_validation_present"
  ),
  passed = c(
    setequal(unique(geometry$centres$assay), mutation_order),
    all(geometry$centres %>% count(assay) %>% pull(n) == 4L),
    all(is.finite(geometry$centres$centre_ch1)) &&
      all(is.finite(geometry$centres$centre_ch2)),
    all(is.finite(as.matrix(geometry$covariance_rows[, c("cov_11", "cov_12", "cov_21", "cov_22")]))),
    all(is.finite(geometry$gate_radii$radius_d2)) &&
      all(geometry$gate_radii$radius_d2 > 0),
    nrow(geometry$baseline_shifts) > 0L,
    nrow(geometry$control_validation) > 0L
  )
)

saveRDS(geometry, path_v2("models", "control_geometry", "control_geometry.rds"))
write_csv(geometry$centres, path_v2("models", "control_geometry", "centroids.csv"))
write_csv(geometry$covariance_rows, path_v2("models", "control_geometry", "covariances.csv"))
saveRDS(geometry$covariance_rows, path_v2("models", "control_geometry", "covariances.rds"))
write_csv(geometry$gate_radii, path_v2("models", "control_geometry", "gate_radii.csv"))
write_csv(geometry$baseline_shifts, path_v2("models", "control_geometry", "baseline_shifts.csv"))
write_csv(geometry$control_validation, path_v2("tables", "control_geometry_validation.csv"))
write_csv(geometry_checks, path_v2("tables", "control_geometry_e2e_checks.csv"))

plot_manifest <- save_control_geometry_plots(shared_droplets, geometry)
write_csv(plot_manifest, path_v2("tables", "control_geometry_plot_manifest.csv"))

plot_checks <- tibble(
  check = c("all_plot_paths_present", "all_plot_files_exist"),
  passed = c(
    nrow(plot_manifest) == length(mutation_order) &&
      all(c("svg_path", "pdf_path") %in% names(plot_manifest)),
    all(file.exists(c(plot_manifest$svg_path, plot_manifest$pdf_path)))
  )
)
geometry_checks <- bind_rows(geometry_checks, plot_checks)
write_csv(geometry_checks, path_v2("tables", "control_geometry_e2e_checks.csv"))
stopifnot(all(geometry_checks$passed))

cat("centroid_rows=", nrow(geometry$centres), "\n", sep = "")
cat("covariance_rows=", nrow(geometry$covariance_rows), "\n", sep = "")
cat("gate_radius_rows=", nrow(geometry$gate_radii), "\n", sep = "")
cat("baseline_shift_rows=", nrow(geometry$baseline_shifts), "\n", sep = "")
cat("control_validation_rows=", nrow(geometry$control_validation), "\n", sep = "")
cat("e2e_checks_passed=", sum(geometry_checks$passed), "\n", sep = "")
