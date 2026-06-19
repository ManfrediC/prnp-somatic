source("scratch/ddpcr_gating_method_comparison_v2/code/lib_v2.R")

library(dpcR)

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

run_limit <- parse_int_env("ONE_CHANNEL_RUN_LIMIT", NA_integer_)
plot_sample_per_group <- parse_int_env("ONE_CHANNEL_PLOT_SAMPLE_PER_GROUP", 300L)
random_seed <- 20260606L
threshold_intervals <- c(0.995, 0.999, 0.9995)

set.seed(random_seed)

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

channel_call_definetherain <- function(amplitude, negative_upper, positive_lower) {
  case_when(
    amplitude <= negative_upper & amplitude < positive_lower ~ "negative",
    amplitude >= positive_lower & amplitude > negative_upper ~ "positive",
    TRUE ~ "rain"
  )
}

combine_channel_calls <- function(ch1_call, ch2_call, ch1_role, ch2_role) {
  ref_call <- if (ch1_role == "reference") ch1_call else ch2_call
  mut_call <- if (ch1_role == "mutant") ch1_call else ch2_call

  case_when(
    ref_call == "rain" | mut_call == "rain" ~ "Rain",
    ref_call == "positive" & mut_call == "positive" ~ "DP",
    ref_call == "positive" & mut_call == "negative" ~ "WT",
    ref_call == "negative" & mut_call == "positive" ~ "MUT",
    ref_call == "negative" & mut_call == "negative" ~ "NN",
    TRUE ~ "Unclassified"
  )
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
      ch1_call, ch2_call, method_parameters_json
    )
}

fit_definetherain_thresholds <- function(shared_droplets) {
  controls <- shared_droplets %>% filter(sample_type == "positive_control")
  bind_rows(lapply(mutation_order, function(assay_name) {
    assay_controls <- controls %>% filter(assay == assay_name)
    bind_rows(lapply(c("ch1", "ch2"), function(channel_name) {
      amplitude_col <- paste0(channel_name, "_amplitude")
      role_col <- paste0(channel_name, "_role")
      x <- assay_controls[[amplitude_col]]
      x <- x[is.finite(x)]
      if (length(x) < 10L) {
        return(tibble(
          method_id = "definetherain_channel_combined",
          assay = assay_name,
          channel = channel_name,
          channel_role = assay_controls[[role_col]][[1]],
          status = "failed",
          message = "fewer than 10 positive-control amplitudes",
          negative_mean = NA_real_,
          negative_sd = NA_real_,
          positive_mean = NA_real_,
          positive_sd = NA_real_,
          negative_upper = NA_real_,
          positive_lower = NA_real_,
          rain_band_overlap = NA
        ))
      }
      km <- kmeans(x, centers = 2)
      centres <- sort(km$centers[, 1])
      low_cluster <- which.min(km$centers[, 1])
      high_cluster <- which.max(km$centers[, 1])
      low <- x[km$cluster == low_cluster]
      high <- x[km$cluster == high_cluster]
      negative_mean <- mean(low)
      negative_sd <- sd(low)
      positive_mean <- mean(high)
      positive_sd <- sd(high)
      negative_upper <- negative_mean + 3 * negative_sd
      positive_lower <- positive_mean - 3 * positive_sd
      tibble(
        method_id = "definetherain_channel_combined",
        assay = assay_name,
        channel = channel_name,
        channel_role = assay_controls[[role_col]][[1]],
        status = "ok",
        message = NA_character_,
        negative_mean = negative_mean,
        negative_sd = negative_sd,
        positive_mean = positive_mean,
        positive_sd = positive_sd,
        negative_upper = negative_upper,
        positive_lower = positive_lower,
        rain_band_overlap = negative_upper >= positive_lower
      )
    }))
  }))
}

classify_definetherain_run <- function(run_droplets, thresholds) {
  assay_name <- run_droplets$assay[[1]]
  ch1 <- thresholds %>% filter(assay == assay_name, channel == "ch1")
  ch2 <- thresholds %>% filter(assay == assay_name, channel == "ch2")
  if (nrow(ch1) != 1L || nrow(ch2) != 1L || any(c(ch1$status, ch2$status) != "ok")) {
    stop("definetherain threshold unavailable")
  }

  run_droplets %>%
    mutate(
      ch1_call = channel_call_definetherain(ch1_amplitude, ch1$negative_upper, ch1$positive_lower),
      ch2_call = channel_call_definetherain(ch2_amplitude, ch2$negative_upper, ch2$positive_lower),
      assigned_target_class = factor(
        combine_channel_calls(ch1_call, ch2_call, ch1_role[[1]], ch2_role[[1]]),
        levels = target_class_levels
      )
    )
}

classify_threshold_run <- function(run_droplets, ch1_threshold, ch2_threshold) {
  run_droplets %>%
    mutate(
      ch1_call = ifelse(ch1_amplitude >= ch1_threshold, "positive", "negative"),
      ch2_call = ifelse(ch2_amplitude >= ch2_threshold, "positive", "negative"),
      assigned_target_class = factor(
        combine_channel_calls(ch1_call, ch2_call, ch1_role[[1]], ch2_role[[1]]),
        levels = target_class_levels
      )
    )
}

run_ddpcrquant_channel <- function(path, threshold_interval) {
  ddpcRquant(
    path,
    threshold.int = threshold_interval,
    reps = 10,
    threshold.manual = FALSE
  )
}

method_grid <- bind_rows(
  tibble(
    method_id = "definetherain_channel_combined",
    method_variant = "local definetherain one-channel rain bands combined into two-channel calls",
    classifier = "definetherain",
    threshold_interval = NA_real_,
    reps = NA_integer_
  ),
  tibble(
    method_id = paste0("ddpcRquant_", c("0995", "0999", "09995")),
    method_variant = paste0("ddpcRquant threshold.int=", threshold_intervals, ", reps=10"),
    classifier = "ddpcRquant",
    threshold_interval = threshold_intervals,
    reps = 10L
  )
)
write_csv(method_grid, path_v2("tables", "one_channel_parameter_grid.csv"))

shared_droplets <- load_shared_droplets()
definetherain_thresholds <- fit_definetherain_thresholds(shared_droplets)

package_manifest <- read_csv(path_v2("tables", "package_input_manifest.csv"), show_col_types = FALSE) %>%
  filter(package == "ddpcRquant") %>%
  arrange(assay, run_id, channel)
run_manifest <- package_manifest %>%
  distinct(assay, run_id)
if (!is.na(run_limit)) {
  run_manifest <- run_manifest %>% slice_head(n = run_limit)
  package_manifest <- package_manifest %>%
    semi_join(run_manifest, by = c("assay", "run_id"))
}

all_count_rows <- list()
all_validation_rows <- list()
all_sample_rows <- list()
all_threshold_rows <- list(definetherain_thresholds)
all_status_rows <- list()

for (i in seq_len(nrow(run_manifest))) {
  run_row <- run_manifest[i, ]
  assay_name <- run_row$assay
  run_id <- run_row$run_id
  twoddpcr_input <- path_v2("inputs", "twoddpcr", assay_name, safe_file_component(run_id), "droplets.rds")
  run_droplets <- readRDS(twoddpcr_input)

  def_parameters <- list(source = "positive_control", sd_multiplier = 3)
  def_result <- tryCatch(
    classify_definetherain_run(run_droplets, definetherain_thresholds),
    error = function(err) err
  )
  if (inherits(def_result, "error")) {
    all_count_rows <- append(all_count_rows, list(failure_count_rows(
      run_droplets,
      "definetherain_channel_combined",
      "local definetherain one-channel rain bands combined into two-channel calls",
      conditionMessage(def_result),
      def_parameters
    )))
  } else {
    all_count_rows <- append(all_count_rows, list(well_count_rows(
      def_result,
      "definetherain_channel_combined",
      "local definetherain one-channel rain bands combined into two-channel calls",
      def_result$assigned_target_class,
      parameters = def_parameters
    )))
    all_validation_rows <- append(all_validation_rows, list(control_validation_rows(
      def_result,
      "definetherain_channel_combined",
      "local definetherain one-channel rain bands combined into two-channel calls",
      def_parameters
    )))
    all_sample_rows <- append(all_sample_rows, list(diagnostic_sample_rows(
      def_result,
      "definetherain_channel_combined",
      "local definetherain one-channel rain bands combined into two-channel calls",
      def_parameters
    )))
  }

  for (threshold_interval in threshold_intervals) {
    suffix <- c("0.995" = "0995", "0.999" = "0999", "0.9995" = "09995")[[as.character(threshold_interval)]]
    method_id <- paste0("ddpcRquant_", suffix)
    method_variant <- paste0("ddpcRquant threshold.int=", threshold_interval, ", reps=10")
    channel_thresholds <- list()
    channel_status <- list()

    for (channel_name in c("ch1", "ch2")) {
      input_row <- package_manifest %>%
        filter(assay == assay_name, run_id == run_id, channel == channel_name) %>%
        slice(1)
      result <- run_captured(run_ddpcrquant_channel(input_row$path, threshold_interval))
      threshold_value <- if (identical(result$status, "ok")) {
        methods::slot(result$value, "threshold")
      } else {
        NA_real_
      }
      channel_thresholds[[channel_name]] <- threshold_value
      channel_status[[channel_name]] <- result

      all_threshold_rows <- append(all_threshold_rows, list(tibble(
        method_id = method_id,
        assay = assay_name,
        run_id = run_id,
        channel = channel_name,
        channel_role = run_droplets[[paste0(channel_name, "_role")]][[1]],
        status = result$status,
        message = result$message,
        threshold_interval = threshold_interval,
        reps = 10L,
        threshold = threshold_value,
        elapsed_seconds = result$elapsed_seconds
      )))
      all_status_rows <- append(all_status_rows, list(tibble(
        method_id = method_id,
        method_variant = method_variant,
        assay = assay_name,
        run_id = run_id,
        channel = channel_name,
        status = result$status,
        message = result$message,
        elapsed_seconds = result$elapsed_seconds,
        output_excerpt = result$output_excerpt
      )))
    }

    parameters <- list(threshold_interval = threshold_interval, reps = 10)
    if (any(!is.finite(unlist(channel_thresholds)))) {
      messages <- paste(vapply(channel_status, function(x) x$message, character(1)), collapse = "; ")
      all_count_rows <- append(all_count_rows, list(failure_count_rows(
        run_droplets,
        method_id,
        method_variant,
        messages,
        parameters
      )))
    } else {
      classified <- classify_threshold_run(
        run_droplets,
        ch1_threshold = channel_thresholds$ch1,
        ch2_threshold = channel_thresholds$ch2
      )
      all_count_rows <- append(all_count_rows, list(well_count_rows(
        classified,
        method_id,
        method_variant,
        classified$assigned_target_class,
        parameters = parameters
      )))
      all_validation_rows <- append(all_validation_rows, list(control_validation_rows(
        classified,
        method_id,
        method_variant,
        parameters
      )))
      all_sample_rows <- append(all_sample_rows, list(diagnostic_sample_rows(
        classified,
        method_id,
        method_variant,
        parameters
      )))
    }
  }

  cat("processed ", assay_name, " ", run_id, "\n", sep = "")
}

well_counts <- bind_rows(all_count_rows) %>%
  arrange(method_id, assay, run_id, well)
control_validation <- bind_rows(all_validation_rows) %>%
  arrange(method_id, assay, sample_type, quanta_class, assigned_target_class)
thresholds <- bind_rows(all_threshold_rows) %>%
  arrange(method_id, assay, run_id, channel)
run_status <- bind_rows(all_status_rows) %>%
  arrange(method_id, assay, run_id, channel)
plot_droplets <- bind_rows(all_sample_rows)

write_csv(well_counts, path_v2("tables", "one_channel_well_counts.csv"))
write_csv(control_validation, path_v2("tables", "one_channel_control_validation.csv"))
write_csv(thresholds, path_v2("tables", "one_channel_thresholds.csv"))
write_csv(run_status, path_v2("tables", "one_channel_run_status.csv"))
write_csv(definetherain_thresholds, path_v2("models", "definetherain", "channel_models.csv"))
dir.create(path_v2("data", "droplets"), recursive = TRUE, showWarnings = FALSE)
saveRDS(plot_droplets, path_v2("data", "droplets", "one_channel_plot_droplets.rds"))

expected_methods <- method_grid$method_id
expected_wells <- bind_rows(lapply(seq_len(nrow(run_manifest)), function(i) {
  run_row <- run_manifest[i, ]
  read_csv(
    path_v2("inputs", "twoddpcr", run_row$assay, safe_file_component(run_row$run_id), "well_manifest.csv"),
    show_col_types = FALSE
  ) %>%
    transmute(run_id, assay, well)
})) %>%
  distinct(run_id, assay, well)

e2e_checks <- tibble(
  check = c(
    "all_methods_present",
    "well_count_rows_complete",
    "definetherain_thresholds_finite",
    "ddpcRquant_thresholds_present",
    "control_validation_present",
    "failed_rows_have_messages",
    "plot_sample_present"
  ),
  passed = c(
    setequal(unique(well_counts$method_id), expected_methods),
    nrow(well_counts) == nrow(expected_wells) * length(expected_methods),
    all(definetherain_thresholds$status == "ok") &&
      all(is.finite(definetherain_thresholds$negative_upper)) &&
      all(is.finite(definetherain_thresholds$positive_lower)),
    all(is.finite(thresholds$threshold[thresholds$method_id != "definetherain_channel_combined"])),
    nrow(control_validation) > 0L,
    all(!is.na(well_counts$classification_message[well_counts$classification_status == "failed"])),
    nrow(plot_droplets) > 0L
  )
)
write_csv(e2e_checks, path_v2("tables", "one_channel_e2e_checks.csv"))
stopifnot(all(e2e_checks$passed))

cat("one_channel_methods=", length(expected_methods), "\n", sep = "")
cat("one_channel_well_count_rows=", nrow(well_counts), "\n", sep = "")
cat("one_channel_control_validation_rows=", nrow(control_validation), "\n", sep = "")
cat("one_channel_threshold_rows=", nrow(thresholds), "\n", sep = "")
cat("one_channel_run_status_rows=", nrow(run_status), "\n", sep = "")
cat("one_channel_plot_sample_rows=", nrow(plot_droplets), "\n", sep = "")
