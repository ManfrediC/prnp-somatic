source("scratch/ddpcr_gating_method_comparison_v2/code/lib_v2.R")

library(twoddpcr)

setup_v2_dirs()

twoddpcr_classes <- c("NN", "NP", "PN", "PP")
parse_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }
  as.integer(value)
}

parse_int_vector_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

knn_values <- parse_int_vector_env("TWODDPCR_K_VALUES", c(3L, 5L, 11L, 21L))
max_training_per_class <- parse_int_env("TWODDPCR_MAX_TRAINING_PER_CLASS", 100L)
knn_chunk_size <- parse_int_env("TWODDPCR_KNN_CHUNK_SIZE", 75000L)
plot_sample_per_group <- parse_int_env("TWODDPCR_PLOT_SAMPLE_PER_GROUP", 400L)
run_limit <- parse_int_env("TWODDPCR_RUN_LIMIT", NA_integer_)
random_seed <- 20260606L

set.seed(random_seed)

json_text <- function(x) {
  as.character(jsonlite::toJSON(x, auto_unbox = TRUE))
}

control_centres_matrix <- function(geometry, assay_name) {
  rows <- geometry$centres %>%
    filter(assay == assay_name) %>%
    arrange(match(twoddpcr_class, twoddpcr_classes))
  centres <- as.matrix(rows[, c("centre_ch1", "centre_ch2")])
  rownames(centres) <- rows$twoddpcr_class
  centres
}

control_rain_distances <- function(geometry, assay_name) {
  geometry$gate_radii %>%
    filter(assay == assay_name) %>%
    left_join(
      geometry$centres %>% select(assay, target_class, twoddpcr_class),
      by = c("assay", "target_class")
    ) %>%
    arrange(match(twoddpcr_class, twoddpcr_classes)) %>%
    { as.list(setNames(.$radius_d2, .$twoddpcr_class)) }
}

as_twoddpcr_target_class <- function(channel_class, droplets) {
  target_class_from_twoddpcr(
    channel_class = channel_class,
    ch1_role = droplets$ch1_role[[1]],
    ch2_role = droplets$ch2_role[[1]]
  )
}

well_count_rows <- function(droplets, method_id, method_variant, channel_class,
                            status = "ok", message = NA_character_,
                            parameters = list()) {
  assigned_target_class <- as_twoddpcr_target_class(channel_class, droplets)

  droplets %>%
    mutate(
      assigned_twoddpcr_class = as.character(channel_class),
      assigned_target_class = factor(assigned_target_class, levels = target_class_levels)
    ) %>%
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

control_validation_rows <- function(classified, method_id, method_variant, parameters = list()) {
  classified %>%
    filter(sample_type %in% c("NTC", "WT_control", "positive_control")) %>%
    count(
      method_id = method_id,
      method_variant = method_variant,
      assay,
      sample_type,
      quanta_class = target_class,
      assigned_target_class,
      assigned_twoddpcr_class,
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
      assigned_twoddpcr_class, method_parameters_json
    )
}

classified_table <- function(droplets, channel_class) {
  droplets %>%
    mutate(
      assigned_twoddpcr_class = as.character(channel_class),
      assigned_target_class = factor(
        as_twoddpcr_target_class(channel_class, droplets),
        levels = target_class_levels
      )
    )
}

run_kmeans_variant <- function(droplets, centres = NULL) {
  if (is.null(centres)) {
    kmeansClassify(droplets, fullTable = TRUE)
  } else {
    kmeansClassify(droplets, centres = centres, fullTable = TRUE)
  }
}

run_mahalanobis_rain_variant <- function(droplets, channel_class, max_distances) {
  rain_input <- data.frame(
    Ch1.Amplitude = droplets$Ch1.Amplitude,
    Ch2.Amplitude = droplets$Ch2.Amplitude,
    class = factor(as.character(channel_class), levels = c(twoddpcr_classes, "Rain", "N/A"))
  )

  mahalanobisRain(
    rain_input,
    cMethod = "class",
    maxDistances = max_distances,
    fullTable = TRUE
  )$class
}

run_knn_variant <- function(droplets, training, k) {
  chunks <- split(seq_len(nrow(droplets)), ceiling(seq_len(nrow(droplets)) / knn_chunk_size))
  class_out <- vector("character", nrow(droplets))

  for (chunk_indices in chunks) {
    class_out[chunk_indices] <- as.character(knnClassify(
      droplets = droplets[chunk_indices, c("Ch1.Amplitude", "Ch2.Amplitude")],
      trainData = training$data,
      cl = training$labels,
      k = k,
      fullTable = FALSE
    ))
  }

  factor(class_out, levels = c(twoddpcr_classes, "Rain", "N/A"))
}

build_twoddpcr_training_sets <- function(shared_droplets, geometry) {
  controls <- shared_droplets %>%
    filter(sample_type %in% c("NTC", "WT_control", "positive_control"))

  cleaned_controls <- bind_rows(lapply(mutation_order, function(assay_name) {
    assay_controls <- controls %>% filter(assay == assay_name)
    classify_by_control_geometry(
      droplets = assay_controls,
      centres = geometry$centres %>% filter(assay == assay_name),
      covariance_rows = geometry$covariance_rows %>% filter(assay == assay_name),
      gate_radii = geometry$gate_radii %>% filter(assay == assay_name)
    )
  }))

  labelled_controls <- cleaned_controls %>%
    filter(assigned_class %in% c("NN", "WT", "MUT", "DP")) %>%
    group_by(assay) %>%
    group_modify(function(.x, .y) {
      .x %>%
        mutate(
          twoddpcr_class = twoddpcr_class_from_target(
            assigned_class,
            ch1_role = ch1_role[[1]],
            ch2_role = ch2_role[[1]]
          )
        )
    }) %>%
    ungroup()

  training_table <- labelled_controls %>%
    group_by(assay, twoddpcr_class) %>%
    group_modify(function(.x, .y) {
      .x %>% slice_sample(n = min(nrow(.x), max_training_per_class))
    }) %>%
    ungroup() %>%
    select(
      assay, run_id, well, sample, sample_type, assigned_class,
      twoddpcr_class, Ch1.Amplitude, Ch2.Amplitude
    )

  training_summary <- training_table %>%
    count(assay, twoddpcr_class, assigned_class, name = "training_rows") %>%
    mutate(
      label_source = "control_geometry_cleaned_controls",
      max_training_per_class = max_training_per_class,
      random_seed = random_seed
    )

  training_sets <- lapply(split(training_table, training_table$assay), function(assay_training) {
    list(
      data = assay_training[, c("Ch1.Amplitude", "Ch2.Amplitude")],
      labels = factor(assay_training$twoddpcr_class, levels = twoddpcr_classes),
      rows = assay_training
    )
  })

  list(
    training_table = training_table,
    training_summary = training_summary,
    training_sets = training_sets
  )
}

record_success <- function(droplets, method_id, method_variant, channel_class,
                           parameters, count_rows, validation_rows, sample_rows) {
  classified <- classified_table(droplets, channel_class)

  list(
    counts = append(
      count_rows,
      list(well_count_rows(droplets, method_id, method_variant, channel_class,
                           parameters = parameters))
    ),
    validation = append(
      validation_rows,
      list(control_validation_rows(classified, method_id, method_variant, parameters))
    ),
    samples = append(
      sample_rows,
      list(diagnostic_sample_rows(classified, method_id, method_variant, parameters))
    )
  )
}

run_status_row <- function(method_id, method_variant, assay_name, run_id,
                           status, message, elapsed_seconds, parameters) {
  tibble(
    method_id = method_id,
    method_variant = method_variant,
    assay = assay_name,
    run_id = run_id,
    status = status,
    message = message,
    elapsed_seconds = elapsed_seconds,
    method_parameters_json = json_text(parameters)
  )
}

run_with_status <- function(expr) {
  elapsed <- system.time({
    result <- tryCatch(
      list(status = "ok", value = eval.parent(substitute(expr)), message = NA_character_),
      error = function(err) {
        list(status = "failed", value = NULL, message = conditionMessage(err))
      }
    )
  })
  result$elapsed_seconds <- max(0, unname(elapsed[["elapsed"]]))
  result
}

geometry <- readRDS(path_v2("models", "control_geometry", "control_geometry.rds"))
shared_droplets <- load_shared_droplets()
training <- build_twoddpcr_training_sets(shared_droplets, geometry)
saveRDS(training, path_v2("models", "twoddpcr", "training_sets.rds"))
write_csv(training$training_summary, path_v2("tables", "twoddpcr_training_summary.csv"))

package_manifest <- read_csv(path_v2("tables", "package_input_manifest.csv"), show_col_types = FALSE) %>%
  filter(package == "twoddpcr") %>%
  arrange(assay, run_id)
if (!is.na(run_limit)) {
  package_manifest <- package_manifest %>% slice_head(n = run_limit)
}

method_grid <- bind_rows(
  tibble(
    method_id = c(
      "twoddpcr_kmeans_native",
      "twoddpcr_kmeans_control_centres",
      "twoddpcr_kmeans_native_mah_rain",
      "twoddpcr_kmeans_control_mah_rain"
    ),
    method_variant = c(
      "native k-means with package default centres",
      "k-means initialised from control geometry",
      "native k-means plus control-derived Mahalanobis rain distances",
      "control-centred k-means plus control-derived Mahalanobis rain distances"
    ),
    classifier = c("kmeans", "kmeans", "kmeans_mahalanobis_rain", "kmeans_mahalanobis_rain"),
    k = NA_integer_,
    max_training_per_class = NA_integer_,
    random_seed = random_seed
  ),
  tibble(
    method_id = paste0("twoddpcr_knn_k", knn_values),
    method_variant = paste0("control-geometry-cleaned kNN, k=", knn_values),
    classifier = "knn",
    k = knn_values,
    max_training_per_class = max_training_per_class,
    random_seed = random_seed
  ),
  tibble(
    method_id = paste0("twoddpcr_knn_k", knn_values, "_mah_rain"),
    method_variant = paste0("control-geometry-cleaned kNN plus Mahalanobis rain, k=", knn_values),
    classifier = "knn_mahalanobis_rain",
    k = knn_values,
    max_training_per_class = max_training_per_class,
    random_seed = random_seed
  )
)
write_csv(method_grid, path_v2("tables", "twoddpcr_parameter_grid.csv"))

all_count_rows <- list()
all_validation_rows <- list()
all_sample_rows <- list()
all_status_rows <- list()

for (i in seq_len(nrow(package_manifest))) {
  manifest_row <- package_manifest[i, ]
  assay_name <- manifest_row$assay
  run_id <- manifest_row$run_id
  droplets <- readRDS(file.path(manifest_row$path, "droplets.rds"))
  centres <- control_centres_matrix(geometry, assay_name)
  max_distances <- control_rain_distances(geometry, assay_name)
  assay_training <- training$training_sets[[assay_name]]

  kmeans_native <- run_with_status(run_kmeans_variant(droplets))
  native_parameters <- list(centres = "package_default")
  all_status_rows <- append(all_status_rows, list(run_status_row(
    "twoddpcr_kmeans_native",
    "native k-means with package default centres",
    assay_name,
    run_id,
    kmeans_native$status,
    kmeans_native$message,
    kmeans_native$elapsed_seconds,
    native_parameters
  )))

  if (identical(kmeans_native$status, "ok")) {
    native_record <- record_success(
      droplets = droplets,
      method_id = "twoddpcr_kmeans_native",
      method_variant = "native k-means with package default centres",
      channel_class = kmeans_native$value$data$class,
      parameters = native_parameters,
      count_rows = all_count_rows,
      validation_rows = all_validation_rows,
      sample_rows = all_sample_rows
    )
    all_count_rows <- native_record$counts
    all_validation_rows <- native_record$validation
    all_sample_rows <- native_record$samples

    native_rain_parameters <- list(
      base_classifier = "twoddpcr_kmeans_native",
      max_distances = max_distances
    )
    native_rain <- run_with_status(run_mahalanobis_rain_variant(
      droplets,
      kmeans_native$value$data$class,
      max_distances
    ))
    all_status_rows <- append(all_status_rows, list(run_status_row(
      "twoddpcr_kmeans_native_mah_rain",
      "native k-means plus control-derived Mahalanobis rain distances",
      assay_name,
      run_id,
      native_rain$status,
      native_rain$message,
      native_rain$elapsed_seconds,
      native_rain_parameters
    )))

    if (identical(native_rain$status, "ok")) {
      native_rain_record <- record_success(
        droplets = droplets,
        method_id = "twoddpcr_kmeans_native_mah_rain",
        method_variant = "native k-means plus control-derived Mahalanobis rain distances",
        channel_class = native_rain$value,
        parameters = native_rain_parameters,
        count_rows = all_count_rows,
        validation_rows = all_validation_rows,
        sample_rows = all_sample_rows
      )
      all_count_rows <- native_rain_record$counts
      all_validation_rows <- native_rain_record$validation
      all_sample_rows <- native_rain_record$samples
    } else {
      all_count_rows <- append(all_count_rows, list(failure_count_rows(
        droplets,
        "twoddpcr_kmeans_native_mah_rain",
        "native k-means plus control-derived Mahalanobis rain distances",
        native_rain$message,
        native_rain_parameters
      )))
    }
  } else {
    all_count_rows <- append(all_count_rows, list(failure_count_rows(
      droplets,
      "twoddpcr_kmeans_native",
      "native k-means with package default centres",
      kmeans_native$message,
      native_parameters
    )))
    all_status_rows <- append(all_status_rows, list(run_status_row(
      "twoddpcr_kmeans_native_mah_rain",
      "native k-means plus control-derived Mahalanobis rain distances",
      assay_name,
      run_id,
      "failed",
      paste("Base classifier failed:", kmeans_native$message),
      0,
      list(base_classifier = "twoddpcr_kmeans_native", max_distances = max_distances)
    )))
    all_count_rows <- append(all_count_rows, list(failure_count_rows(
      droplets,
      "twoddpcr_kmeans_native_mah_rain",
      "native k-means plus control-derived Mahalanobis rain distances",
      paste("Base classifier failed:", kmeans_native$message),
      list(base_classifier = "twoddpcr_kmeans_native", max_distances = max_distances)
    )))
  }

  control_parameters <- list(centres = "control_geometry")
  kmeans_control <- run_with_status(run_kmeans_variant(droplets, centres = centres))
  all_status_rows <- append(all_status_rows, list(run_status_row(
    "twoddpcr_kmeans_control_centres",
    "k-means initialised from control geometry",
    assay_name,
    run_id,
    kmeans_control$status,
    kmeans_control$message,
    kmeans_control$elapsed_seconds,
    control_parameters
  )))

  if (identical(kmeans_control$status, "ok")) {
    control_record <- record_success(
      droplets = droplets,
      method_id = "twoddpcr_kmeans_control_centres",
      method_variant = "k-means initialised from control geometry",
      channel_class = kmeans_control$value$data$class,
      parameters = control_parameters,
      count_rows = all_count_rows,
      validation_rows = all_validation_rows,
      sample_rows = all_sample_rows
    )
    all_count_rows <- control_record$counts
    all_validation_rows <- control_record$validation
    all_sample_rows <- control_record$samples

    control_rain_parameters <- list(
      base_classifier = "twoddpcr_kmeans_control_centres",
      max_distances = max_distances
    )
    control_rain <- run_with_status(run_mahalanobis_rain_variant(
      droplets,
      kmeans_control$value$data$class,
      max_distances
    ))
    all_status_rows <- append(all_status_rows, list(run_status_row(
      "twoddpcr_kmeans_control_mah_rain",
      "control-centred k-means plus control-derived Mahalanobis rain distances",
      assay_name,
      run_id,
      control_rain$status,
      control_rain$message,
      control_rain$elapsed_seconds,
      control_rain_parameters
    )))

    if (identical(control_rain$status, "ok")) {
      control_rain_record <- record_success(
        droplets = droplets,
        method_id = "twoddpcr_kmeans_control_mah_rain",
        method_variant = "control-centred k-means plus control-derived Mahalanobis rain distances",
        channel_class = control_rain$value,
        parameters = control_rain_parameters,
        count_rows = all_count_rows,
        validation_rows = all_validation_rows,
        sample_rows = all_sample_rows
      )
      all_count_rows <- control_rain_record$counts
      all_validation_rows <- control_rain_record$validation
      all_sample_rows <- control_rain_record$samples
    } else {
      all_count_rows <- append(all_count_rows, list(failure_count_rows(
        droplets,
        "twoddpcr_kmeans_control_mah_rain",
        "control-centred k-means plus control-derived Mahalanobis rain distances",
        control_rain$message,
        control_rain_parameters
      )))
    }
  } else {
    all_count_rows <- append(all_count_rows, list(failure_count_rows(
      droplets,
      "twoddpcr_kmeans_control_centres",
      "k-means initialised from control geometry",
      kmeans_control$message,
      control_parameters
    )))
    all_status_rows <- append(all_status_rows, list(run_status_row(
      "twoddpcr_kmeans_control_mah_rain",
      "control-centred k-means plus control-derived Mahalanobis rain distances",
      assay_name,
      run_id,
      "failed",
      paste("Base classifier failed:", kmeans_control$message),
      0,
      list(base_classifier = "twoddpcr_kmeans_control_centres", max_distances = max_distances)
    )))
    all_count_rows <- append(all_count_rows, list(failure_count_rows(
      droplets,
      "twoddpcr_kmeans_control_mah_rain",
      "control-centred k-means plus control-derived Mahalanobis rain distances",
      paste("Base classifier failed:", kmeans_control$message),
      list(base_classifier = "twoddpcr_kmeans_control_centres", max_distances = max_distances)
    )))
  }

  for (k_value in knn_values) {
    knn_method_id <- paste0("twoddpcr_knn_k", k_value)
    knn_variant <- paste0("control-geometry-cleaned kNN, k=", k_value)
    knn_parameters <- list(
      k = k_value,
      training_source = "control_geometry_cleaned_controls",
      max_training_per_class = max_training_per_class,
      random_seed = random_seed
    )
    knn_result <- run_with_status(run_knn_variant(droplets, assay_training, k_value))
    all_status_rows <- append(all_status_rows, list(run_status_row(
      knn_method_id,
      knn_variant,
      assay_name,
      run_id,
      knn_result$status,
      knn_result$message,
      knn_result$elapsed_seconds,
      knn_parameters
    )))

    if (identical(knn_result$status, "ok")) {
      knn_record <- record_success(
        droplets = droplets,
        method_id = knn_method_id,
        method_variant = knn_variant,
        channel_class = knn_result$value,
        parameters = knn_parameters,
        count_rows = all_count_rows,
        validation_rows = all_validation_rows,
        sample_rows = all_sample_rows
      )
      all_count_rows <- knn_record$counts
      all_validation_rows <- knn_record$validation
      all_sample_rows <- knn_record$samples

      knn_rain_method_id <- paste0(knn_method_id, "_mah_rain")
      knn_rain_variant <- paste0(
        "control-geometry-cleaned kNN plus Mahalanobis rain, k=",
        k_value
      )
      knn_rain_parameters <- c(knn_parameters, list(max_distances = max_distances))
      knn_rain_result <- run_with_status(run_mahalanobis_rain_variant(
        droplets,
        knn_result$value,
        max_distances
      ))
      all_status_rows <- append(all_status_rows, list(run_status_row(
        knn_rain_method_id,
        knn_rain_variant,
        assay_name,
        run_id,
        knn_rain_result$status,
        knn_rain_result$message,
        knn_rain_result$elapsed_seconds,
        knn_rain_parameters
      )))

      if (identical(knn_rain_result$status, "ok")) {
        knn_rain_record <- record_success(
          droplets = droplets,
          method_id = knn_rain_method_id,
          method_variant = knn_rain_variant,
          channel_class = knn_rain_result$value,
          parameters = knn_rain_parameters,
          count_rows = all_count_rows,
          validation_rows = all_validation_rows,
          sample_rows = all_sample_rows
        )
        all_count_rows <- knn_rain_record$counts
        all_validation_rows <- knn_rain_record$validation
        all_sample_rows <- knn_rain_record$samples
      } else {
        all_count_rows <- append(all_count_rows, list(failure_count_rows(
          droplets,
          knn_rain_method_id,
          knn_rain_variant,
          knn_rain_result$message,
          knn_rain_parameters
        )))
      }
    } else {
      all_count_rows <- append(all_count_rows, list(failure_count_rows(
        droplets,
        knn_method_id,
        knn_variant,
        knn_result$message,
        knn_parameters
      )))
      all_status_rows <- append(all_status_rows, list(run_status_row(
        paste0(knn_method_id, "_mah_rain"),
        paste0("control-geometry-cleaned kNN plus Mahalanobis rain, k=", k_value),
        assay_name,
        run_id,
        "failed",
        paste("Base classifier failed:", knn_result$message),
        0,
        c(knn_parameters, list(max_distances = max_distances))
      )))
      all_count_rows <- append(all_count_rows, list(failure_count_rows(
        droplets,
        paste0(knn_method_id, "_mah_rain"),
        paste0("control-geometry-cleaned kNN plus Mahalanobis rain, k=", k_value),
        paste("Base classifier failed:", knn_result$message),
        c(knn_parameters, list(max_distances = max_distances))
      )))
    }
  }

  cat("processed ", assay_name, " ", run_id, "\n", sep = "")
}

well_counts <- bind_rows(all_count_rows) %>%
  arrange(method_id, assay, run_id, well)
control_validation <- bind_rows(all_validation_rows) %>%
  arrange(method_id, assay, sample_type, quanta_class, assigned_target_class)
plot_droplets <- bind_rows(all_sample_rows)
run_status <- bind_rows(all_status_rows) %>%
  arrange(method_id, assay, run_id)

write_csv(well_counts, path_v2("tables", "twoddpcr_well_counts.csv"))
write_csv(control_validation, path_v2("tables", "twoddpcr_control_validation.csv"))
write_csv(run_status, path_v2("tables", "twoddpcr_run_status.csv"))
dir.create(path_v2("data", "droplets"), recursive = TRUE, showWarnings = FALSE)
saveRDS(plot_droplets, path_v2("data", "droplets", "twoddpcr_plot_droplets.rds"))

expected_methods <- method_grid$method_id
expected_wells <- bind_rows(lapply(package_manifest$path, function(run_path) {
  read_csv(file.path(run_path, "well_manifest.csv"), show_col_types = FALSE)
})) %>%
  distinct(run_id, assay, well)
training_class_balance <- training$training_summary %>%
  count(assay, twoddpcr_class, name = "label_rows") %>%
  count(assay, name = "classes_present")

e2e_checks <- tibble(
  check = c(
    "all_methods_present",
    "well_count_rows_complete",
    "run_status_rows_present",
    "training_has_four_classes_per_assay",
    "control_validation_present",
    "failed_rows_have_messages",
    "plot_sample_present"
  ),
  passed = c(
    setequal(unique(well_counts$method_id), expected_methods),
    nrow(well_counts) == nrow(expected_wells) * length(expected_methods),
    nrow(run_status) >= nrow(package_manifest) * length(expected_methods),
    all(training_class_balance$classes_present == length(twoddpcr_classes)),
    nrow(control_validation) > 0L,
    all(!is.na(well_counts$classification_message[well_counts$classification_status == "failed"])),
    nrow(plot_droplets) > 0L
  )
)
write_csv(e2e_checks, path_v2("tables", "twoddpcr_e2e_checks.csv"))
stopifnot(all(e2e_checks$passed))

cat("twoddpcr_methods=", length(expected_methods), "\n", sep = "")
cat("twoddpcr_well_count_rows=", nrow(well_counts), "\n", sep = "")
cat("twoddpcr_control_validation_rows=", nrow(control_validation), "\n", sep = "")
cat("twoddpcr_run_status_rows=", nrow(run_status), "\n", sep = "")
cat("twoddpcr_plot_sample_rows=", nrow(plot_droplets), "\n", sep = "")
