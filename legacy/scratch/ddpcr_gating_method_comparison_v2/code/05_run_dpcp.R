source("scratch/ddpcr_gating_method_comparison_v2/code/lib_v2.R")

library(dPCP)

setup_v2_dirs()

parse_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }
  as.integer(value)
}

json_text <- function(x) {
  as.character(jsonlite::toJSON(x, auto_unbox = TRUE))
}

run_limit <- parse_int_env("DPCP_RUN_LIMIT", NA_integer_)
plot_sample_per_group <- parse_int_env("DPCP_PLOT_SAMPLE_PER_GROUP", 300L)
random_seed <- 20260606L

set.seed(random_seed)

dpcp_parameter_grid <- tibble(
  attempt = 1:5,
  eps = c(200, 150, 250, 300, 200),
  minPts = c(50, 30, 50, 50, 30)
)

run_captured <- function(expr) {
  output <- character()
  elapsed <- system.time({
    result <- tryCatch({
      output <- capture.output(value <- suppressWarnings(eval.parent(substitute(expr))))
      list(status = "ok", value = value, message = NA_character_)
    }, error = function(err) {
      list(status = "failed", value = NULL, message = conditionMessage(err))
    })
  })
  result$elapsed_seconds <- max(0, unname(elapsed[["elapsed"]]))
  result$output_excerpt <- paste(utils::head(output, 25), collapse = " | ")
  result
}

well_count_rows <- function(droplets, method_id, method_variant, assigned_target_class,
                            status = "ok", message = NA_character_,
                            parameters = list()) {
  droplets %>%
    mutate(assigned_target_class = factor(assigned_target_class, levels = target_class_levels)) %>%
    group_by(run_id, assay, well, sample, sample_type, sample_key, run_date) %>%
    summarise(
      method_id = method_id,
      method_variant = method_variant,
      n_total_droplets = n(),
      n_accepted_droplets = sum(assigned_target_class %in% c("NN", "WT", "MUT", "DP")),
      n_nn = sum(assigned_target_class == "NN"),
      n_wt_only = sum(assigned_target_class == "WT"),
      n_mut_only = sum(assigned_target_class == "MUT"),
      n_double_positive = sum(assigned_target_class == "DP"),
      n_rain = sum(assigned_target_class == "Rain"),
      n_unclassified = sum(assigned_target_class == "Unclassified"),
      n_mut_expected = n_mut_only + n_double_positive,
      n_wt_expected = n_wt_only + n_double_positive,
      classification_status = status,
      classification_message = message,
      method_parameters_json = json_text(parameters),
      .groups = "drop"
    )
}

failure_count_rows <- function(droplets, method_id, method_variant, message,
                               parameters = list()) {
  droplets %>%
    group_by(run_id, assay, well, sample, sample_type, sample_key, run_date) %>%
    summarise(
      method_id = method_id,
      method_variant = method_variant,
      n_total_droplets = n(),
      n_accepted_droplets = NA_real_,
      n_nn = NA_real_,
      n_wt_only = NA_real_,
      n_mut_only = NA_real_,
      n_double_positive = NA_real_,
      n_rain = NA_real_,
      n_unclassified = NA_real_,
      n_mut_expected = NA_real_,
      n_wt_expected = NA_real_,
      classification_status = "failed",
      classification_message = message,
      method_parameters_json = json_text(parameters),
      .groups = "drop"
    )
}

control_validation_rows <- function(classified, method_id, method_variant,
                                    parameters = list()) {
  classified %>%
    filter(sample_type %in% c("NTC", "WT_control", "positive_control")) %>%
    count(
      method_id = method_id,
      method_variant = method_variant,
      assay,
      sample_type,
      quanta_class = target_class,
      assigned_target_class,
      package_cluster,
      name = "droplets"
    ) %>%
    mutate(method_parameters_json = json_text(parameters))
}

diagnostic_sample_rows <- function(classified, method_id, method_variant,
                                   parameters = list()) {
  classified %>%
    mutate(
      method_id = method_id,
      method_variant = method_variant,
      method_parameters_json = json_text(parameters)
    ) %>%
    group_by(method_id, assay, sample_type, assigned_target_class) %>%
    group_modify(function(.x, .y) {
      .x %>% slice_sample(n = min(nrow(.x), plot_sample_per_group))
    }) %>%
    ungroup() %>%
    select(
      method_id, method_variant, assay, run_id, well, sample, sample_type,
      ch1_amplitude, ch2_amplitude, target_class, assigned_target_class,
      package_cluster, method_parameters_json
    )
}

run_status_row <- function(method_id, method_variant, assay_name, run_id,
                           attempt, eps, minPts, status, message,
                           elapsed_seconds, output_excerpt) {
  tibble(
    method_id = method_id,
    method_variant = method_variant,
    assay = assay_name,
    run_id = run_id,
    attempt = attempt,
    eps = eps,
    minPts = minPts,
    status = status,
    message = message,
    elapsed_seconds = elapsed_seconds,
    output_excerpt = output_excerpt
  )
}

dpcp_label_to_target_class <- function(cluster_label, wt_target, mut_target) {
  label <- as.character(cluster_label)
  label[is.na(label)] <- ""
  has_wt <- stringr::str_detect(label, fixed(wt_target))
  has_mut <- stringr::str_detect(label, fixed(mut_target))

  case_when(
    label == "Empty" ~ "NN",
    stringr::str_detect(stringr::str_to_lower(label), "rain") ~ "Rain",
    has_wt & has_mut ~ "DP",
    has_wt ~ "WT",
    has_mut ~ "MUT",
    TRUE ~ "Rain"
  )
}

native_dpcp_run <- function(sample_table_path, run_path, eps, minPts) {
  dPCP(
    file = sample_table_path,
    system = "bio-rad",
    file.location = run_path,
    reference.quality = 0.5,
    sample.quality = 0.5,
    eps = eps,
    minPts = minPts,
    save.template = FALSE,
    rain = TRUE
  )
}

reader_validation_row <- function(manifest_row) {
  sample_table_path <- file.path(manifest_row$path, "sample_table.csv")
  sample_table <- run_captured(read_sampleTable(
    sample_table_path,
    system = "bio-rad",
    file.location = manifest_row$path
  ))

  reference <- if (identical(sample_table$status, "ok")) {
    run_captured(read_reference(
      sample_table$value,
      system = "bio-rad",
      file.location = manifest_row$path,
      eps = 200,
      minPts = 50
    ))
  } else {
    list(status = "failed", message = "sample table read failed", elapsed_seconds = 0, output_excerpt = "")
  }

  sample <- if (identical(sample_table$status, "ok")) {
    run_captured(read_sample(
      sample_table$value,
      system = "bio-rad",
      file.location = manifest_row$path
    ))
  } else {
    list(status = "failed", message = "sample table read failed", elapsed_seconds = 0, output_excerpt = "")
  }

  tibble(
    assay = manifest_row$assay,
    run_id = manifest_row$run_id,
    path = manifest_row$path,
    sample_table_status = sample_table$status,
    sample_table_message = sample_table$message,
    reference_status = reference$status,
    reference_message = reference$message,
    sample_status = sample$status,
    sample_message = sample$message,
    sample_table_rows = if (identical(sample_table$status, "ok")) nrow(sample_table$value) else NA_integer_,
    reference_files = if (identical(reference$status, "ok")) length(reference$value) else NA_integer_,
    sample_files = if (identical(sample$status, "ok")) length(sample$value) else NA_integer_
  )
}

classified_from_dpcp_result <- function(dpcp_result, run_droplets, sample_table) {
  sample_table <- as_tibble(unclass(sample_table))
  classified_wells <- list()
  failed_wells <- list()

  for (i in seq_len(nrow(sample_table))) {
    sample_name <- sample_table$Sample.name[[i]]
    well_name <- sample_table$Chip.ID.Well.ID[[i]]
    result_name <- paste(sample_name, well_name, sep = "_")
    well_droplets <- run_droplets %>% filter(well == well_name)
    sample_result <- dpcp_result$samples[[result_name]]

    if (is.null(sample_result) || is.null(sample_result$data)) {
      failed_wells <- append(failed_wells, list(list(
        well = well_name,
        message = paste("Missing dPCP sample result:", result_name)
      )))
      next
    }

    data <- sample_result$data
    cluster_col <- if ("final cluster" %in% names(data)) {
      "final cluster"
    } else if ("cmeans cluster" %in% names(data)) {
      "cmeans cluster"
    } else {
      NA_character_
    }
    if (is.na(cluster_col)) {
      failed_wells <- append(failed_wells, list(list(
        well = well_name,
        message = paste("Missing dPCP cluster column:", result_name)
      )))
      next
    }
    if (nrow(data) != nrow(well_droplets)) {
      failed_wells <- append(failed_wells, list(list(
        well = well_name,
        message = paste("dPCP row count mismatch:", result_name)
      )))
      next
    }

    wt_target <- sample_table$FAM.target[[i]]
    mut_target <- sample_table$VIC.HEX.target[[i]]
    assigned <- dpcp_label_to_target_class(data[[cluster_col]], wt_target, mut_target)

    classified_wells <- append(classified_wells, list(well_droplets %>%
      mutate(
        package_cluster = as.character(data[[cluster_col]]),
        assigned_target_class = factor(assigned, levels = target_class_levels)
      )))
  }

  list(
    classified = bind_rows(classified_wells),
    failed = failed_wells
  )
}

control_projected_classified <- function(run_droplets, geometry, assay_name) {
  classify_by_control_geometry(
    droplets = run_droplets,
    centres = geometry$centres %>% filter(assay == assay_name),
    covariance_rows = geometry$covariance_rows %>% filter(assay == assay_name),
    gate_radii = geometry$gate_radii %>% filter(assay == assay_name)
  ) %>%
    mutate(
      package_cluster = assigned_class,
      assigned_target_class = factor(assigned_class, levels = target_class_levels)
    )
}

geometry <- readRDS(path_v2("models", "control_geometry", "control_geometry.rds"))
package_manifest <- read_csv(path_v2("tables", "package_input_manifest.csv"), show_col_types = FALSE) %>%
  filter(package == "dPCP") %>%
  arrange(assay, run_id)
if (!is.na(run_limit)) {
  package_manifest <- package_manifest %>% slice_head(n = run_limit)
}

write_csv(dpcp_parameter_grid, path_v2("tables", "dpcp_parameter_grid.csv"))
reference_selection <- package_manifest %>%
  select(
    assay, run_id, path, sample_table_path, reference_well,
    reference_file, reference_sample_type
  )
write_csv(reference_selection, path_v2("tables", "dpcp_reference_selection.csv"))

reader_validation <- bind_rows(lapply(seq_len(nrow(package_manifest)), function(i) {
  reader_validation_row(package_manifest[i, ])
}))
write_csv(reader_validation, path_v2("tables", "dpcp_reader_validation.csv"))

all_count_rows <- list()
all_validation_rows <- list()
all_sample_rows <- list()
all_status_rows <- list()

for (i in seq_len(nrow(package_manifest))) {
  manifest_row <- package_manifest[i, ]
  assay_name <- manifest_row$assay
  run_id <- manifest_row$run_id
  run_path <- manifest_row$path
  sample_table_path <- file.path(run_path, "sample_table.csv")
  twoddpcr_input <- path_v2("inputs", "twoddpcr", assay_name, safe_file_component(run_id), "droplets.rds")
  run_droplets <- readRDS(twoddpcr_input)
  sample_table <- read_sampleTable(sample_table_path, system = "bio-rad", file.location = run_path)

  native_result <- NULL
  native_parameters <- NULL
  for (attempt_i in seq_len(nrow(dpcp_parameter_grid))) {
    params <- dpcp_parameter_grid[attempt_i, ]
    result <- run_captured(native_dpcp_run(
      sample_table_path,
      run_path,
      eps = params$eps,
      minPts = params$minPts
    ))
    all_status_rows <- append(all_status_rows, list(run_status_row(
      "dPCP_native",
      "native dPCP() with fixed DBSCAN retry grid",
      assay_name,
      run_id,
      params$attempt,
      params$eps,
      params$minPts,
      result$status,
      result$message,
      result$elapsed_seconds,
      result$output_excerpt
    )))

    if (identical(result$status, "ok")) {
      native_result <- result$value
      native_parameters <- list(attempt = params$attempt, eps = params$eps, minPts = params$minPts)
      break
    }
  }

  if (!is.null(native_result)) {
    native_classified <- classified_from_dpcp_result(native_result, run_droplets, sample_table)
    if (nrow(native_classified$classified) > 0L) {
      all_count_rows <- append(all_count_rows, list(well_count_rows(
        native_classified$classified,
        "dPCP_native",
        "native dPCP() with fixed DBSCAN retry grid",
        native_classified$classified$assigned_target_class,
        parameters = native_parameters
      )))
      all_validation_rows <- append(all_validation_rows, list(control_validation_rows(
        native_classified$classified,
        "dPCP_native",
        "native dPCP() with fixed DBSCAN retry grid",
        native_parameters
      )))
      all_sample_rows <- append(all_sample_rows, list(diagnostic_sample_rows(
        native_classified$classified,
        "dPCP_native",
        "native dPCP() with fixed DBSCAN retry grid",
        native_parameters
      )))
    }
    if (length(native_classified$failed) > 0L) {
      for (failed in native_classified$failed) {
        all_count_rows <- append(all_count_rows, list(failure_count_rows(
          run_droplets %>% filter(well == failed$well),
          "dPCP_native",
          "native dPCP() with fixed DBSCAN retry grid",
          failed$message,
          native_parameters
        )))
      }
    }
  } else {
    last_status <- all_status_rows[[length(all_status_rows)]]
    all_count_rows <- append(all_count_rows, list(failure_count_rows(
      run_droplets,
      "dPCP_native",
      "native dPCP() with fixed DBSCAN retry grid",
      last_status$message,
      list(attempts = dpcp_parameter_grid)
    )))
  }

  control_parameters <- list(source = "control_geometry", package_input = "dPCP")
  control_classified <- control_projected_classified(run_droplets, geometry, assay_name)
  all_status_rows <- append(all_status_rows, list(run_status_row(
    "dPCP_control_projected",
    "control-projected adapter-style baseline using dPCP inputs",
    assay_name,
    run_id,
    NA_integer_,
    NA_real_,
    NA_real_,
    "ok",
    NA_character_,
    0,
    NA_character_
  )))
  all_count_rows <- append(all_count_rows, list(well_count_rows(
    control_classified,
    "dPCP_control_projected",
    "control-projected adapter-style baseline using dPCP inputs",
    control_classified$assigned_target_class,
    parameters = control_parameters
  )))
  all_validation_rows <- append(all_validation_rows, list(control_validation_rows(
    control_classified,
    "dPCP_control_projected",
    "control-projected adapter-style baseline using dPCP inputs",
    control_parameters
  )))
  all_sample_rows <- append(all_sample_rows, list(diagnostic_sample_rows(
    control_classified,
    "dPCP_control_projected",
    "control-projected adapter-style baseline using dPCP inputs",
    control_parameters
  )))

  cat("processed ", assay_name, " ", run_id, "\n", sep = "")
}

well_counts <- bind_rows(all_count_rows) %>%
  arrange(method_id, assay, run_id, well)
control_validation <- bind_rows(all_validation_rows) %>%
  arrange(method_id, assay, sample_type, quanta_class, assigned_target_class)
plot_droplets <- bind_rows(all_sample_rows)
run_status <- bind_rows(all_status_rows) %>%
  arrange(method_id, assay, run_id, attempt)

adapter_comparison <- well_counts %>%
  select(
    method_id, assay, run_id, well, sample, sample_type,
    classification_status, n_mut_expected, n_wt_expected, n_rain
  ) %>%
  pivot_wider(
    names_from = method_id,
    values_from = c(classification_status, n_mut_expected, n_wt_expected, n_rain)
  ) %>%
  mutate(
    delta_native_minus_control_mut =
      n_mut_expected_dPCP_native - n_mut_expected_dPCP_control_projected,
    delta_native_minus_control_wt =
      n_wt_expected_dPCP_native - n_wt_expected_dPCP_control_projected
  )

write_csv(well_counts, path_v2("tables", "dpcp_well_counts.csv"))
write_csv(control_validation, path_v2("tables", "dpcp_control_validation.csv"))
write_csv(run_status, path_v2("tables", "dpcp_run_status.csv"))
write_csv(adapter_comparison, path_v2("tables", "dpcp_adapter_comparison.csv"))
dir.create(path_v2("data", "droplets"), recursive = TRUE, showWarnings = FALSE)
saveRDS(plot_droplets, path_v2("data", "droplets", "dPCP_plot_droplets.rds"))

expected_wells <- bind_rows(lapply(seq_len(nrow(package_manifest)), function(i) {
  read_csv(package_manifest$sample_table_path[[i]], show_col_types = FALSE) %>%
    transmute(
      run_id = package_manifest$run_id[[i]],
      assay = package_manifest$assay[[i]],
      well = `Chip ID/Well ID`
    )
})) %>%
  distinct(run_id, assay, well)
expected_methods <- c("dPCP_native", "dPCP_control_projected")

e2e_checks <- tibble(
  check = c(
    "all_methods_present",
    "well_count_rows_complete",
    "reader_validation_passed",
    "native_status_present",
    "control_validation_present",
    "adapter_comparison_complete",
    "failed_rows_have_messages",
    "plot_sample_present"
  ),
  passed = c(
    setequal(unique(well_counts$method_id), expected_methods),
    nrow(well_counts) == nrow(expected_wells) * length(expected_methods),
    all(reader_validation$sample_table_status == "ok") &&
      all(reader_validation$reference_status == "ok") &&
      all(reader_validation$sample_status == "ok"),
    all(package_manifest$run_id %in% run_status$run_id[run_status$method_id == "dPCP_native"]),
    nrow(control_validation) > 0L,
    nrow(adapter_comparison) == nrow(expected_wells),
    all(!is.na(well_counts$classification_message[well_counts$classification_status == "failed"])),
    nrow(plot_droplets) > 0L
  )
)
write_csv(e2e_checks, path_v2("tables", "dpcp_e2e_checks.csv"))
stopifnot(all(e2e_checks$passed))

cat("dpcp_methods=", length(expected_methods), "\n", sep = "")
cat("dpcp_well_count_rows=", nrow(well_counts), "\n", sep = "")
cat("dpcp_control_validation_rows=", nrow(control_validation), "\n", sep = "")
cat("dpcp_run_status_rows=", nrow(run_status), "\n", sep = "")
cat("dpcp_adapter_comparison_rows=", nrow(adapter_comparison), "\n", sep = "")
cat("dpcp_plot_sample_rows=", nrow(plot_droplets), "\n", sep = "")
