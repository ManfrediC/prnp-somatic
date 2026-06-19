source("scratch/ddpcr_gating_method_comparison_v2/code/lib_v2.R")

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

run_limit <- parse_int_env("BAYESIAN_RUN_LIMIT", NA_integer_)
plot_sample_per_group <- parse_int_env("BAYESIAN_PLOT_SAMPLE_PER_GROUP", 300L)
prior_smoothing <- 100
rain_covariance_multiplier <- 25
random_seed <- 20260606L
bayes_classes <- c("NN", "WT", "MUT", "DP", "Rain")

set.seed(random_seed)

covariance_matrix <- function(covariance_rows, assay_name, class_name) {
  row <- covariance_rows %>% filter(assay == assay_name, target_class == class_name) %>% slice(1)
  matrix(c(row$cov_11, row$cov_12, row$cov_21, row$cov_22), nrow = 2, byrow = TRUE)
}

log_dmvnorm <- function(x, mean, sigma) {
  chol_sigma <- chol(sigma)
  solved <- backsolve(chol_sigma, t(x - matrix(mean, nrow = nrow(x), ncol = 2, byrow = TRUE)))
  quad <- colSums(solved^2)
  log_det <- 2 * sum(log(diag(chol_sigma)))
  -0.5 * (2 * log(2 * pi) + log_det + quad)
}

build_class_parameters <- function(geometry, shared_droplets) {
  controls <- shared_droplets %>%
    filter(sample_type %in% c("NTC", "WT_control", "positive_control"))

  bind_rows(lapply(mutation_order, function(assay_name) {
    centres <- geometry$centres %>% filter(assay == assay_name)
    cov_rows <- geometry$covariance_rows %>% filter(assay == assay_name)
    biological <- bind_rows(lapply(c("NN", "WT", "MUT", "DP"), function(class_name) {
      cov_mat <- covariance_matrix(cov_rows, assay_name, class_name)
      centre <- centres %>% filter(target_class == class_name)
      tibble(
        assay = assay_name,
        target_class = class_name,
        mean_ch1 = centre$centre_ch1,
        mean_ch2 = centre$centre_ch2,
        cov_11 = cov_mat[1, 1],
        cov_12 = cov_mat[1, 2],
        cov_21 = cov_mat[2, 1],
        cov_22 = cov_mat[2, 2],
        component_source = "control_geometry"
      )
    }))

    pooled_cov <- Reduce("+", lapply(c("NN", "WT", "MUT", "DP"), function(class_name) {
      covariance_matrix(cov_rows, assay_name, class_name)
    })) / 4
    centre_mat <- as.matrix(centres[, c("centre_ch1", "centre_ch2")])
    rain_cov <- pooled_cov * rain_covariance_multiplier + cov(centre_mat) + diag(c(1000, 1000), 2)
    rain_mean <- colMeans(centre_mat)

    bind_rows(
      biological,
      tibble(
        assay = assay_name,
        target_class = "Rain",
        mean_ch1 = rain_mean[1],
        mean_ch2 = rain_mean[2],
        cov_11 = rain_cov[1, 1],
        cov_12 = rain_cov[1, 2],
        cov_21 = rain_cov[2, 1],
        cov_22 = rain_cov[2, 2],
        component_source = "broad_control_geometry_rain"
      )
    )
  }))
}

build_priors <- function(geometry, shared_droplets) {
  controls <- shared_droplets %>%
    filter(sample_type %in% c("NTC", "WT_control"))

  bind_rows(lapply(mutation_order, function(assay_name) {
    assigned <- classify_by_control_geometry(
      droplets = controls %>% filter(assay == assay_name),
      centres = geometry$centres %>% filter(assay == assay_name),
      covariance_rows = geometry$covariance_rows %>% filter(assay == assay_name),
      gate_radii = geometry$gate_radii %>% filter(assay == assay_name)
    )
    counts <- tibble(target_class = bayes_classes) %>%
      left_join(
        assigned %>% count(target_class = assigned_class, name = "observed"),
        by = "target_class"
      ) %>%
      mutate(
        observed = replace_na(observed, 0),
        smoothed = observed + prior_smoothing,
        prior = smoothed / sum(smoothed),
        assay = assay_name,
        prior_source = "ntc_wt_control_geometry_assignments_smoothed"
      )
    counts
  })) %>%
    select(assay, target_class, observed, smoothed, prior, prior_source)
}

component_covariance <- function(parameters, assay_name, class_name) {
  row <- parameters %>% filter(assay == assay_name, target_class == class_name) %>% slice(1)
  matrix(c(row$cov_11, row$cov_12, row$cov_21, row$cov_22), nrow = 2, byrow = TRUE)
}

posterior_probabilities <- function(droplets, parameters, priors, assay_name) {
  x <- as.matrix(droplets[, c("ch1_amplitude", "ch2_amplitude")])
  log_lik <- sapply(bayes_classes, function(class_name) {
    row <- parameters %>% filter(assay == assay_name, target_class == class_name) %>% slice(1)
    log_dmvnorm(
      x,
      mean = c(row$mean_ch1, row$mean_ch2),
      sigma = component_covariance(parameters, assay_name, class_name)
    )
  })
  colnames(log_lik) <- bayes_classes
  log_prior <- priors$prior[match(bayes_classes, priors$target_class)]
  log_weight <- sweep(log_lik, 2, log(log_prior), "+")
  row_max <- apply(log_weight, 1, max)
  weight <- exp(log_weight - row_max)
  posterior <- weight / rowSums(weight)
  colnames(posterior) <- bayes_classes
  posterior
}

weighted_count_rows <- function(droplets, posterior, parameters_json) {
  posterior_df <- as_tibble(posterior)
  droplets %>%
    bind_cols(posterior_df) %>%
    group_by(run_id, assay, well, sample, sample_type, sample_key, run_date) %>%
    summarise(
      method_id = "bayesian_control_mixture_weighted",
      method_variant = "control-anchored Gaussian mixture with probability-weighted counts",
      n_total_droplets = n(),
      n_accepted_droplets = sum(NN + WT + MUT + DP),
      n_nn = sum(NN),
      n_wt_only = sum(WT),
      n_mut_only = sum(MUT),
      n_double_positive = sum(DP),
      n_rain = sum(Rain),
      n_unclassified = 0,
      n_mut_expected = sum(MUT + DP),
      n_wt_expected = sum(WT + DP),
      classification_status = "ok",
      classification_message = NA_character_,
      method_parameters_json = parameters_json,
      .groups = "drop"
    )
}

hard_count_rows <- function(droplets, assigned_target_class, parameters_json) {
  droplets %>%
    mutate(assigned_target_class = factor(assigned_target_class, levels = target_class_levels)) %>%
    group_by(run_id, assay, well, sample, sample_type, sample_key, run_date) %>%
    summarise(
      method_id = "bayesian_control_mixture_hard_map",
      method_variant = "control-anchored Gaussian mixture with hard MAP calls",
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
      classification_status = "ok",
      classification_message = NA_character_,
      method_parameters_json = parameters_json,
      .groups = "drop"
    )
}

control_validation_rows <- function(classified, method_id, method_variant, parameters_json) {
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
    mutate(method_parameters_json = parameters_json)
}

diagnostic_sample_rows <- function(classified, method_id, method_variant, parameters_json) {
  classified %>%
    mutate(
      method_id = method_id,
      method_variant = method_variant,
      method_parameters_json = parameters_json
    ) %>%
    group_by(method_id, assay, sample_type, assigned_target_class) %>%
    group_modify(function(.x, .y) {
      .x %>% slice_sample(n = min(nrow(.x), plot_sample_per_group))
    }) %>%
    ungroup() %>%
    select(
      method_id, method_variant, assay, run_id, well, sample, sample_type,
      ch1_amplitude, ch2_amplitude, target_class, assigned_target_class,
      max_posterior, posterior_entropy, method_parameters_json
    )
}

geometry <- readRDS(path_v2("models", "control_geometry", "control_geometry.rds"))
shared_droplets <- load_shared_droplets()
class_parameters <- build_class_parameters(geometry, shared_droplets)
priors <- build_priors(geometry, shared_droplets)

saveRDS(class_parameters, path_v2("models", "bayesian_mixture", "class_parameters.rds"))
write_csv(class_parameters, path_v2("models", "bayesian_mixture", "class_parameters.csv"))
write_csv(priors, path_v2("models", "bayesian_mixture", "priors.csv"))

method_parameters <- json_text(list(
  class_parameters = "control_geometry",
  rain_covariance_multiplier = rain_covariance_multiplier,
  prior_smoothing = prior_smoothing
))

package_manifest <- read_csv(path_v2("tables", "package_input_manifest.csv"), show_col_types = FALSE) %>%
  filter(package == "twoddpcr") %>%
  arrange(assay, run_id)
if (!is.na(run_limit)) {
  package_manifest <- package_manifest %>% slice_head(n = run_limit)
}

all_count_rows <- list()
all_validation_rows <- list()
all_sample_rows <- list()
all_uncertainty_rows <- list()
all_status_rows <- list()

for (i in seq_len(nrow(package_manifest))) {
  manifest_row <- package_manifest[i, ]
  assay_name <- manifest_row$assay
  run_id <- manifest_row$run_id
  droplets <- readRDS(file.path(manifest_row$path, "droplets.rds"))
  assay_priors <- priors %>% filter(assay == assay_name)

  elapsed <- system.time({
    posterior <- posterior_probabilities(droplets, class_parameters, assay_priors, assay_name)
  })
  posterior_row_sum <- rowSums(posterior)
  max_posterior <- apply(posterior, 1, max)
  assigned <- bayes_classes[max.col(posterior, ties.method = "first")]
  entropy <- -rowSums(posterior * log(pmax(posterior, .Machine$double.xmin)))

  all_count_rows <- append(all_count_rows, list(weighted_count_rows(
    droplets,
    posterior,
    method_parameters
  )))
  all_count_rows <- append(all_count_rows, list(hard_count_rows(
    droplets,
    assigned,
    method_parameters
  )))

  hard_classified <- droplets %>%
    mutate(
      assigned_target_class = factor(assigned, levels = target_class_levels),
      max_posterior = max_posterior,
      posterior_entropy = entropy
    )

  all_validation_rows <- append(all_validation_rows, list(control_validation_rows(
    hard_classified,
    "bayesian_control_mixture_hard_map",
    "control-anchored Gaussian mixture with hard MAP calls",
    method_parameters
  )))
  all_sample_rows <- append(all_sample_rows, list(diagnostic_sample_rows(
    hard_classified,
    "bayesian_control_mixture_hard_map",
    "control-anchored Gaussian mixture with hard MAP calls",
    method_parameters
  )))
  all_uncertainty_rows <- append(all_uncertainty_rows, list(hard_classified %>%
    group_by(run_id, assay, well, sample, sample_type, sample_key) %>%
    summarise(
      mean_max_posterior = mean(max_posterior),
      median_max_posterior = median(max_posterior),
      mean_entropy = mean(posterior_entropy),
      posterior_rowsum_min = min(posterior_row_sum),
      posterior_rowsum_max = max(posterior_row_sum),
      low_confidence_droplets = sum(max_posterior < 0.8),
      .groups = "drop"
    )))
  all_status_rows <- append(all_status_rows, list(tibble(
    assay = assay_name,
    run_id = run_id,
    status = "ok",
    message = NA_character_,
    elapsed_seconds = max(0, unname(elapsed[["elapsed"]])),
    posterior_rowsum_min = min(posterior_row_sum),
    posterior_rowsum_max = max(posterior_row_sum)
  )))

  cat("processed ", assay_name, " ", run_id, "\n", sep = "")
}

well_counts <- bind_rows(all_count_rows) %>%
  arrange(method_id, assay, run_id, well)
control_validation <- bind_rows(all_validation_rows) %>%
  arrange(method_id, assay, sample_type, quanta_class, assigned_target_class)
plot_droplets <- bind_rows(all_sample_rows)
uncertainty <- bind_rows(all_uncertainty_rows) %>%
  arrange(assay, run_id, well)
run_status <- bind_rows(all_status_rows) %>%
  arrange(assay, run_id)

write_csv(well_counts, path_v2("tables", "bayesian_well_counts.csv"))
write_csv(control_validation, path_v2("tables", "bayesian_control_validation.csv"))
write_csv(uncertainty, path_v2("tables", "bayesian_uncertainty.csv"))
write_csv(run_status, path_v2("tables", "bayesian_run_status.csv"))
dir.create(path_v2("data", "droplets"), recursive = TRUE, showWarnings = FALSE)
saveRDS(plot_droplets, path_v2("data", "droplets", "bayesian_mixture_plot_droplets.rds"))

expected_wells <- bind_rows(lapply(seq_len(nrow(package_manifest)), function(i) {
  read_csv(file.path(package_manifest$path[[i]], "well_manifest.csv"), show_col_types = FALSE) %>%
    transmute(run_id, assay, well)
})) %>%
  distinct(run_id, assay, well)
expected_methods <- c("bayesian_control_mixture_weighted", "bayesian_control_mixture_hard_map")

e2e_checks <- tibble(
  check = c(
    "all_methods_present",
    "well_count_rows_complete",
    "finite_class_parameters",
    "priors_sum_to_one",
    "posterior_rows_sum_to_one",
    "control_validation_present",
    "plot_sample_present"
  ),
  passed = c(
    setequal(unique(well_counts$method_id), expected_methods),
    nrow(well_counts) == nrow(expected_wells) * length(expected_methods),
    all(is.finite(class_parameters$mean_ch1)) &&
      all(is.finite(class_parameters$mean_ch2)) &&
      all(is.finite(as.matrix(class_parameters[, c("cov_11", "cov_12", "cov_21", "cov_22")]))),
    all(abs((priors %>% group_by(assay) %>% summarise(total = sum(prior), .groups = "drop"))$total - 1) < 1e-9),
    all(abs(run_status$posterior_rowsum_min - 1) < 1e-8) &&
      all(abs(run_status$posterior_rowsum_max - 1) < 1e-8),
    nrow(control_validation) > 0L,
    nrow(plot_droplets) > 0L
  )
)
write_csv(e2e_checks, path_v2("tables", "bayesian_e2e_checks.csv"))
stopifnot(all(e2e_checks$passed))

cat("bayesian_methods=", length(expected_methods), "\n", sep = "")
cat("bayesian_well_count_rows=", nrow(well_counts), "\n", sep = "")
cat("bayesian_control_validation_rows=", nrow(control_validation), "\n", sep = "")
cat("bayesian_uncertainty_rows=", nrow(uncertainty), "\n", sep = "")
cat("bayesian_run_status_rows=", nrow(run_status), "\n", sep = "")
cat("bayesian_plot_sample_rows=", nrow(plot_droplets), "\n", sep = "")
