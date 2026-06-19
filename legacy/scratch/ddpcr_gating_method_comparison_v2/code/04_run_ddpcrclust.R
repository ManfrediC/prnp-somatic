source("scratch/ddpcr_gating_method_comparison_v2/code/lib_v2.R")

library(ddPCRclust)

setup_v2_dirs()

parse_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }
  as.integer(value)
}

parse_modes_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }
  strsplit(value, ",", fixed = TRUE)[[1]]
}

json_text <- function(x) {
  as.character(jsonlite::toJSON(x, auto_unbox = TRUE))
}

template_modes <- parse_modes_env("DDPCRCLUST_TEMPLATE_MODES", c("full", "fast"))
run_limit <- parse_int_env("DDPCRCLUST_RUN_LIMIT", NA_integer_)
plot_sample_per_group <- parse_int_env("DDPCRCLUST_PLOT_SAMPLE_PER_GROUP", 300L)
random_seed <- 20260606L

set.seed(random_seed)

mode_to_method_id <- function(mode) {
  if (identical(mode, "full")) {
    "ddPCRclust_template_native"
  } else if (identical(mode, "fast")) {
    "ddPCRclust_template_fast"
  } else {
    stop("Unsupported ddPCRclust template mode: ", mode)
  }
}

mode_to_variant <- function(mode) {
  if (identical(mode, "full")) {
    "native template workflow, full ensemble"
  } else {
    "native template workflow, fast density mode"
  }
}

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
                           status, message, elapsed_seconds, output_excerpt,
                           parameters = list()) {
  tibble(
    method_id = method_id,
    method_variant = method_variant,
    assay = assay_name,
    run_id = run_id,
    status = status,
    message = message,
    elapsed_seconds = elapsed_seconds,
    output_excerpt = output_excerpt,
    method_parameters_json = json_text(parameters)
  )
}

read_ddpcrclust_input <- function(run_path) {
  old_wd <- setwd(run_path)
  on.exit(setwd(old_wd), add = TRUE)
  amplitude_files <- list.files(pattern = "_Amplitude.csv$")
  files <- readFiles(amplitude_files)
  template <- readTemplate("ddPCR_run_template.csv")
  list(files = files, template = template, amplitude_files = amplitude_files)
}

reader_validation_row <- function(manifest_row) {
  result <- run_captured(read_ddpcrclust_input(manifest_row$path))
  if (identical(result$status, "ok")) {
    ids <- result$value$files$ids
    template_rows <- nrow(result$value$template$template)
  } else {
    ids <- character()
    template_rows <- NA_integer_
  }

  tibble(
    assay = manifest_row$assay,
    run_id = manifest_row$run_id,
    path = manifest_row$path,
    status = result$status,
    message = result$message,
    amplitude_files = length(ids),
    template_rows = template_rows,
    output_excerpt = result$output_excerpt
  )
}

cluster_labels_from_geometry <- function(result_data, geometry, assay_name) {
  cluster_centres <- result_data %>%
    mutate(package_cluster = as.character(Cluster)) %>%
    group_by(package_cluster) %>%
    summarise(
      ch1_amplitude = median(Ch1.Amplitude, na.rm = TRUE),
      ch2_amplitude = median(Ch2.Amplitude, na.rm = TRUE),
      .groups = "drop"
    )

  centre_class <- classify_by_control_geometry(
    droplets = cluster_centres,
    centres = geometry$centres %>% filter(assay == assay_name),
    covariance_rows = geometry$covariance_rows %>% filter(assay == assay_name),
    gate_radii = geometry$gate_radii %>% filter(assay == assay_name)
  ) %>%
    select(package_cluster, assigned_class)

  centre_class$assigned_class[match(as.character(result_data$Cluster), centre_class$package_cluster)]
}

classified_from_package_result <- function(well_droplets, well_result, geometry, assay_name) {
  result_data <- well_result$data
  if (!all(c("Ch1.Amplitude", "Ch2.Amplitude", "Cluster") %in% names(result_data))) {
    stop("ddPCRclust result lacks expected Ch1.Amplitude, Ch2.Amplitude, and Cluster columns")
  }
  if (nrow(result_data) != nrow(well_droplets)) {
    stop("ddPCRclust result row count does not match source droplet count")
  }

  assigned <- cluster_labels_from_geometry(result_data, geometry, assay_name)
  well_droplets %>%
    mutate(
      package_cluster = as.character(result_data$Cluster),
      assigned_target_class = factor(assigned, levels = target_class_levels)
    )
}

run_template_mode <- function(run_path, mode) {
  input <- read_ddpcrclust_input(run_path)
  ddPCRclust(
    files = input$files,
    template = input$template,
    numOfMarkers = 2,
    fast = identical(mode, "fast"),
    multithread = FALSE
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

representative_low_level_inputs <- function(manifest) {
  run_inputs <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i, ]
    read_csv(file.path(row$path, "ddPCR_run_template.csv"), skip = 1, show_col_types = FALSE) %>%
      transmute(
        assay = row$assay,
        run_id = row$run_id,
        path = row$path,
        well = Well,
        sample_type = `Sample type`,
        amplitude_path = file.path(row$path, paste0(Well, "_Amplitude.csv"))
      )
  }))

  run_inputs %>%
    filter(sample_type %in% c("NTC", "WT_control", "positive_control")) %>%
    group_by(assay, sample_type) %>%
    slice_head(n = 1) %>%
    ungroup()
}

run_low_level_smoke <- function(smoke_inputs) {
  methods <- c("runDensity", "runPeaks", "runSam")

  bind_rows(lapply(seq_len(nrow(smoke_inputs)), function(i) {
    input <- smoke_inputs[i, ]
    file_data <- read.csv(input$amplitude_path)
    bind_rows(lapply(methods, function(method_name) {
      result <- run_captured(do.call(
        method_name,
        list(file = file_data, sensitivity = 1, numOfMarkers = 2)
      ))
      tibble(
        assay = input$assay,
        run_id = input$run_id,
        well = input$well,
        sample_type = input$sample_type,
        method = method_name,
        status = result$status,
        message = result$message,
        elapsed_seconds = result$elapsed_seconds,
        output_excerpt = result$output_excerpt,
        result_rows = ifelse(
          identical(result$status, "ok") && !is.null(result$value$data),
          nrow(result$value$data),
          NA_integer_
        )
      )
    }))
  }))
}

geometry <- readRDS(path_v2("models", "control_geometry", "control_geometry.rds"))
package_manifest <- read_csv(path_v2("tables", "package_input_manifest.csv"), show_col_types = FALSE) %>%
  filter(package == "ddPCRclust") %>%
  arrange(assay, run_id)
if (!is.na(run_limit)) {
  package_manifest <- package_manifest %>% slice_head(n = run_limit)
}

method_grid <- bind_rows(
  tibble(
    method_id = vapply(template_modes, mode_to_method_id, character(1)),
    method_variant = vapply(template_modes, mode_to_variant, character(1)),
    classifier = "ddPCRclust_template",
    template_mode = template_modes,
    random_seed = random_seed
  ),
  tibble(
    method_id = "ddPCRclust_control_projected",
    method_variant = "control-projected sensitivity using ddPCRclust input geometry",
    classifier = "control_projected",
    template_mode = NA_character_,
    random_seed = random_seed
  )
)
write_csv(method_grid, path_v2("tables", "ddPCRclust_parameter_grid.csv"))

reader_validation <- bind_rows(lapply(seq_len(nrow(package_manifest)), function(i) {
  reader_validation_row(package_manifest[i, ])
}))
write_csv(reader_validation, path_v2("tables", "ddPCRclust_reader_validation.csv"))

all_count_rows <- list()
all_validation_rows <- list()
all_sample_rows <- list()
all_status_rows <- list()

for (i in seq_len(nrow(package_manifest))) {
  manifest_row <- package_manifest[i, ]
  assay_name <- manifest_row$assay
  run_id <- manifest_row$run_id
  twoddpcr_input <- path_v2("inputs", "twoddpcr", assay_name, safe_file_component(run_id), "droplets.rds")
  run_droplets <- readRDS(twoddpcr_input)

  for (mode in template_modes) {
    method_id <- mode_to_method_id(mode)
    method_variant <- mode_to_variant(mode)
    parameters <- list(mode = mode, numOfMarkers = 2, multithread = FALSE)
    package_result <- run_captured(run_template_mode(manifest_row$path, mode))

    all_status_rows <- append(all_status_rows, list(run_status_row(
      method_id,
      method_variant,
      assay_name,
      run_id,
      package_result$status,
      package_result$message,
      package_result$elapsed_seconds,
      package_result$output_excerpt,
      parameters
    )))

    if (identical(package_result$status, "ok")) {
      classified_wells <- list()
      failed_wells <- list()

      for (well_name in names(package_result$value)) {
        well_droplets <- run_droplets %>% filter(well == well_name)
        well_result <- package_result$value[[well_name]]
        classified <- tryCatch(
          classified_from_package_result(well_droplets, well_result, geometry, assay_name),
          error = function(err) err
        )

        if (inherits(classified, "error")) {
          failed_wells <- append(failed_wells, list(failure_count_rows(
            well_droplets,
            method_id,
            method_variant,
            conditionMessage(classified),
            parameters
          )))
        } else {
          classified_wells <- append(classified_wells, list(classified))
        }
      }

      if (length(classified_wells) > 0L) {
        classified_run <- bind_rows(classified_wells)
        all_count_rows <- append(all_count_rows, list(well_count_rows(
          classified_run,
          method_id,
          method_variant,
          classified_run$assigned_target_class,
          parameters = parameters
        )))
        all_validation_rows <- append(all_validation_rows, list(control_validation_rows(
          classified_run,
          method_id,
          method_variant,
          parameters
        )))
        all_sample_rows <- append(all_sample_rows, list(diagnostic_sample_rows(
          classified_run,
          method_id,
          method_variant,
          parameters
        )))
      }
      all_count_rows <- c(all_count_rows, failed_wells)
    } else {
      all_count_rows <- append(all_count_rows, list(failure_count_rows(
        run_droplets,
        method_id,
        method_variant,
        package_result$message,
        parameters
      )))
    }
  }

  control_parameters <- list(source = "control_geometry", package_input = "ddPCRclust")
  control_classified <- control_projected_classified(run_droplets, geometry, assay_name)
  all_status_rows <- append(all_status_rows, list(run_status_row(
    "ddPCRclust_control_projected",
    "control-projected sensitivity using ddPCRclust input geometry",
    assay_name,
    run_id,
    "ok",
    NA_character_,
    0,
    NA_character_,
    control_parameters
  )))
  all_count_rows <- append(all_count_rows, list(well_count_rows(
    control_classified,
    "ddPCRclust_control_projected",
    "control-projected sensitivity using ddPCRclust input geometry",
    control_classified$assigned_target_class,
    parameters = control_parameters
  )))
  all_validation_rows <- append(all_validation_rows, list(control_validation_rows(
    control_classified,
    "ddPCRclust_control_projected",
    "control-projected sensitivity using ddPCRclust input geometry",
    control_parameters
  )))
  all_sample_rows <- append(all_sample_rows, list(diagnostic_sample_rows(
    control_classified,
    "ddPCRclust_control_projected",
    "control-projected sensitivity using ddPCRclust input geometry",
    control_parameters
  )))

  cat("processed ", assay_name, " ", run_id, "\n", sep = "")
}

low_level_smoke <- run_low_level_smoke(representative_low_level_inputs(package_manifest))

well_counts <- bind_rows(all_count_rows) %>%
  arrange(method_id, assay, run_id, well)
control_validation <- bind_rows(all_validation_rows) %>%
  arrange(method_id, assay, sample_type, quanta_class, assigned_target_class)
plot_droplets <- bind_rows(all_sample_rows)
run_status <- bind_rows(all_status_rows) %>%
  arrange(method_id, assay, run_id)

write_csv(well_counts, path_v2("tables", "ddPCRclust_well_counts.csv"))
write_csv(control_validation, path_v2("tables", "ddPCRclust_control_validation.csv"))
write_csv(run_status, path_v2("tables", "ddPCRclust_run_status.csv"))
write_csv(low_level_smoke, path_v2("tables", "ddPCRclust_low_level_smoke.csv"))
dir.create(path_v2("data", "droplets"), recursive = TRUE, showWarnings = FALSE)
saveRDS(plot_droplets, path_v2("data", "droplets", "ddPCRclust_plot_droplets.rds"))

expected_methods <- method_grid$method_id
expected_wells <- bind_rows(lapply(seq_len(nrow(package_manifest)), function(i) {
  read_csv(
    file.path(package_manifest$path[[i]], "ddPCR_run_template.csv"),
    skip = 1,
    show_col_types = FALSE
  ) %>%
    transmute(
      run_id = package_manifest$run_id[[i]],
      assay = package_manifest$assay[[i]],
      well = Well
    )
})) %>%
  distinct(run_id, assay, well)

e2e_checks <- tibble(
  check = c(
    "all_methods_present",
    "well_count_rows_complete",
    "reader_validation_passed",
    "run_status_rows_complete",
    "control_validation_present",
    "low_level_smoke_complete",
    "failed_rows_have_messages",
    "plot_sample_present"
  ),
  passed = c(
    setequal(unique(well_counts$method_id), expected_methods),
    nrow(well_counts) == nrow(expected_wells) * length(expected_methods),
    all(reader_validation$status == "ok"),
    nrow(run_status) == nrow(package_manifest) * length(expected_methods),
    nrow(control_validation) > 0L,
    nrow(low_level_smoke) == length(unique(low_level_smoke$assay)) * 3L * 3L,
    all(!is.na(well_counts$classification_message[well_counts$classification_status == "failed"])),
    nrow(plot_droplets) > 0L
  )
)
write_csv(e2e_checks, path_v2("tables", "ddPCRclust_e2e_checks.csv"))
stopifnot(all(e2e_checks$passed))

cat("ddPCRclust_methods=", length(expected_methods), "\n", sep = "")
cat("ddPCRclust_well_count_rows=", nrow(well_counts), "\n", sep = "")
cat("ddPCRclust_control_validation_rows=", nrow(control_validation), "\n", sep = "")
cat("ddPCRclust_run_status_rows=", nrow(run_status), "\n", sep = "")
cat("ddPCRclust_low_level_smoke_rows=", nrow(low_level_smoke), "\n", sep = "")
cat("ddPCRclust_plot_sample_rows=", nrow(plot_droplets), "\n", sep = "")
