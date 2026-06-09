library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(jsonlite)

# ---- paths and inputs ----

find_project_root <- function(start_dir = getwd()) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current_dir, "src", "ddpcr")) &&
        dir.exists(file.path(current_dir, "raw", "ddpcr"))) {
      return(current_dir)
    }
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Could not locate project root from: ", start_dir)
    }
    current_dir <- parent_dir
  }
}

project_root <- find_project_root()
raw_root <- file.path(project_root, "raw", "ddpcr")
validation_dir <- file.path(project_root, "results", "ddPCR", "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "src", "ddpcr", "ddpcr_raw_import_helpers.R"))

default_run_ids <- c(
  "2020-10-12_SNV_E200K",
  "2021-01-26_SNV_E200K",
  "2020-12-01_SNV_P102L",
  "2021-03-29_SNV_D178N",
  "2021-03-29_SNV_E200K",
  "2021-03-29_SNV_P102L",
  "2021-04-06_SNV_E200K",
  "2021-05-10_SNV_E200K"
)

args <- commandArgs(trailingOnly = TRUE)
target_run_ids <- if (length(args) > 0L) args else default_run_ids

# Bio-Rad concentration exports use 0.85 nL partitions, expressed as uL here.
droplet_volume_ul <- 0.00085
confidence_68 <- 0.682689492137

# ---- formatting helpers ----

first_matching_column <- function(columns, pattern) {
  matches <- grep(pattern, columns, value = TRUE)
  if (length(matches) == 0L) {
    return(NA_character_)
  }
  matches[[1]]
}

path_from_manifest_value <- function(value) {
  value <- as.character(value[[1]])
  if (str_detect(value, "^C:\\\\") && dir.exists("/mnt/c")) {
    return(file.path("/mnt/c", str_replace_all(str_sub(value, 4), "\\\\", "/")))
  }
  value
}

format_decimal <- function(x, digits) {
  if (is.na(x) || !is.finite(x)) {
    return("")
  }
  formatC(x, format = "f", digits = digits)
}

format_integer <- function(x) {
  if (is.na(x) || !is.finite(x)) {
    return("")
  }
  as.character(as.integer(round(x)))
}

should_replace_value <- function(old_value, new_value, tolerance = 1e-9) {
  old_value <- old_value %||% ""
  new_value <- new_value %||% ""
  if (identical(old_value, new_value)) {
    return(FALSE)
  }
  old_number <- suppressWarnings(as.numeric(old_value))
  new_number <- suppressWarnings(as.numeric(new_value))
  if (!is.na(old_number) && !is.na(new_number) &&
      abs(old_number - new_number) <= tolerance) {
    return(FALSE)
  }
  TRUE
}

replace_preserving_equivalent <- function(old_values, raw_values, formatter, tolerance = 1e-9) {
  purrr::map2_chr(old_values, raw_values, function(old_value, raw_value) {
    new_value <- formatter(raw_value)
    if (should_replace_value(old_value, new_value, tolerance)) {
      new_value
    } else {
      old_value
    }
  })
}

# ---- archive-derived values ----

read_well_amplitude_rows <- function(raw_root, run_row, well_manifest) {
  extract_path <- file.path(raw_root, run_row$archive_contents_relative_dir)
  assay <- normalise_assay(run_row$assay)
  well <- well_manifest$well[[1]]
  peak_path <- file.path(extract_path, "PeakData", paste0(well, ".ddpeakjson"))
  metadata_path <- file.path(extract_path, "PeakMetaData", paste0(well, ".ddmetajson"))

  if (!file.exists(peak_path)) {
    stop("Missing peak data: ", peak_path)
  }
  if (!file.exists(metadata_path)) {
    stop("Missing peak metadata: ", metadata_path)
  }

  peak <- jsonlite::fromJSON(peak_path, simplifyVector = FALSE)
  metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)
  targets <- metadata_targets(metadata)
  selected <- selected_target_indices(targets, assay)
  if (is.null(targets) || is.null(selected)) {
    stop("Could not select targets for ", run_row$run_id, " ", well)
  }

  target_names <- vapply(targets, target_name, character(1))
  target_clean <- clean_ddpcr_target(target_names)
  channels <- vapply(targets, target_channel, integer(1))
  selected_indices <- c(selected$ref, selected$mut)

  amplitudes <- peak$PeakInfo$Amplitudes
  if (is.null(amplitudes) || length(amplitudes) < 2L) {
    stop("Peak data does not contain two amplitude channels: ", peak_path)
  }
  amplitude_by_channel <- list(
    `1` = as.numeric(amplitudes[[1]]),
    `2` = as.numeric(amplitudes[[2]])
  )

  total_values <- setNames(vector("list", length(selected_indices)), selected_indices)
  positive_values <- setNames(vector("list", length(selected_indices)), selected_indices)
  negative_values <- setNames(vector("list", length(selected_indices)), selected_indices)
  for (idx in as.character(selected_indices)) {
    total_values[[idx]] <- numeric()
    positive_values[[idx]] <- numeric()
    negative_values[[idx]] <- numeric()
  }

  for (cluster in metadata$Clusters %||% list()) {
    droplet_indices <- as.integer(cluster$Droplets %||% integer(0)) + 1L
    if (length(droplet_indices) == 0L) {
      next
    }

    results <- as.character(unlist(cluster$Results, use.names = FALSE))
    if (length(results) < max(selected_indices) || is_true(cluster$Unassigned)) {
      next
    }

    selected_results <- results[selected_indices]
    if (!all(selected_results %in% c("Negative", "Positive"))) {
      next
    }

    for (target_index in selected_indices) {
      channel <- channels[[target_index]]
      values <- amplitude_by_channel[[as.character(channel)]]
      row_indices <- droplet_indices[droplet_indices >= 1L & droplet_indices <= length(values)]
      if (length(row_indices) == 0L) {
        next
      }

      idx <- as.character(target_index)
      channel_values <- values[row_indices]
      total_values[[idx]] <- c(total_values[[idx]], channel_values)
      if (results[[target_index]] == "Positive") {
        positive_values[[idx]] <- c(positive_values[[idx]], channel_values)
      } else {
        negative_values[[idx]] <- c(negative_values[[idx]], channel_values)
      }
    }
  }

  mean_or_missing <- function(x) {
    if (length(x) == 0L) {
      return(NA_real_)
    }
    mean(x, na.rm = TRUE)
  }

  purrr::map_dfr(selected_indices, function(target_index) {
    idx <- as.character(target_index)
    tibble(
      run_id = run_row$run_id,
      Well = well,
      target_clean = target_clean[[target_index]],
      MeanAmplitudeOfPositives_raw = mean_or_missing(positive_values[[idx]]),
      MeanAmplitudeOfNegatives_raw = mean_or_missing(negative_values[[idx]]),
      MeanAmplitudeTotal_raw = mean_or_missing(total_values[[idx]]),
      raw_peak_path = peak_path
    )
  })
}

target_concentration_rows <- function(raw_rows) {
  intervals_95 <- purrr::map2(
    raw_rows$Positives,
    raw_rows$`Accepted Droplets`,
    ~ binomial_lambda_interval(.x, .y, conf.level = 0.95)
  )
  intervals_68 <- purrr::map2(
    raw_rows$Positives,
    raw_rows$`Accepted Droplets`,
    ~ binomial_lambda_interval(.x, .y, conf.level = confidence_68)
  )

  raw_rows %>%
    mutate(
      lambda_target = purrr::map_dbl(intervals_95, "lambda"),
      lambda_target_low_95 = purrr::map_dbl(intervals_95, "lower"),
      lambda_target_high_95 = purrr::map_dbl(intervals_95, "upper"),
      lambda_target_low_68 = purrr::map_dbl(intervals_68, "lower"),
      lambda_target_high_68 = purrr::map_dbl(intervals_68, "upper"),
      concentration_raw = lambda_target / droplet_volume_ul,
      concentration_low_95_raw = lambda_target_low_95 / droplet_volume_ul,
      concentration_high_95_raw = lambda_target_high_95 / droplet_volume_ul,
      concentration_low_68_raw = lambda_target_low_68 / droplet_volume_ul,
      concentration_high_68_raw = lambda_target_high_68 / droplet_volume_ul,
      copies_20ul_raw = concentration_raw * 20
    )
}

well_fractional_abundance_rows <- function(raw_rows) {
  raw_rows %>%
    group_by(run_id, Well) %>%
    summarise(
      accepted = first(`Accepted Droplets`),
      ref_positive = Positives[raw_role == "reference"][[1]],
      mut_positive = Positives[raw_role == "mutant"][[1]],
      ref_negative = Negatives[raw_role == "reference"][[1]],
      mut_negative = Negatives[raw_role == "mutant"][[1]],
      .groups = "drop"
    ) %>%
    mutate(
      fa_95 = pmap(
        list(ref_positive, ref_negative, mut_positive, mut_negative, accepted),
        ~ fractional_abundance(..1, ..2, ..3, ..4, ..5, conf.level = 0.95)
      ),
      fa_68 = pmap(
        list(ref_positive, ref_negative, mut_positive, mut_negative, accepted),
        ~ fractional_abundance(..1, ..2, ..3, ..4, ..5, conf.level = confidence_68)
      )
    ) %>%
    transmute(
      run_id,
      Well,
      linkage_raw = purrr::map_dbl(fa_95, "lambda_mut") / droplet_volume_ul,
      ratio_raw = purrr::map_dbl(fa_95, "ratio"),
      ratio_low_95_raw = purrr::map_dbl(fa_95, "ratio_low"),
      ratio_high_95_raw = purrr::map_dbl(fa_95, "ratio_high"),
      fractional_abundance_raw = purrr::map_dbl(fa_95, "fractional_abundance"),
      fractional_abundance_low_95_raw = purrr::map_dbl(fa_95, "ci_low"),
      fractional_abundance_high_95_raw = purrr::map_dbl(fa_95, "ci_high"),
      ratio_low_68_raw = purrr::map_dbl(fa_68, "ratio_low"),
      ratio_high_68_raw = purrr::map_dbl(fa_68, "ratio_high"),
      fractional_abundance_low_68_raw = purrr::map_dbl(fa_68, "ci_low"),
      fractional_abundance_high_68_raw = purrr::map_dbl(fa_68, "ci_high")
    )
}

build_run_updates <- function(raw_root, run_row, sample_manifest) {
  run_manifest <- sample_manifest %>%
    filter(run_id == run_row$run_id) %>%
    arrange(well_order, target_order)

  raw_rows <- read_archive_rows_from_database(raw_root, run_row, sample_manifest) %>%
    mutate(
      Positives = as.numeric(Positives),
      Negatives = as.numeric(Negatives),
      `Accepted Droplets` = as.numeric(`Accepted Droplets`)
    )

  amplitude_rows <- purrr::map_dfr(unique(run_manifest$well), function(well) {
    read_well_amplitude_rows(
      raw_root = raw_root,
      run_row = run_row,
      well_manifest = run_manifest %>% filter(.data$well == !!well)
    )
  })

  raw_rows %>%
    target_concentration_rows() %>%
    left_join(well_fractional_abundance_rows(raw_rows), by = c("run_id", "Well")) %>%
    left_join(amplitude_rows, by = c("run_id", "Well", "target_clean")) %>%
    mutate(
      ratio_raw = if_else(raw_role == "mutant", ratio_raw, NA_real_),
      ratio_low_95_raw = if_else(raw_role == "mutant", ratio_low_95_raw, NA_real_),
      ratio_high_95_raw = if_else(raw_role == "mutant", ratio_high_95_raw, NA_real_),
      ratio_low_68_raw = if_else(raw_role == "mutant", ratio_low_68_raw, NA_real_),
      ratio_high_68_raw = if_else(raw_role == "mutant", ratio_high_68_raw, NA_real_),
      fractional_abundance_raw = if_else(raw_role == "mutant", fractional_abundance_raw, NA_real_),
      fractional_abundance_low_95_raw = if_else(raw_role == "mutant", fractional_abundance_low_95_raw, NA_real_),
      fractional_abundance_high_95_raw = if_else(raw_role == "mutant", fractional_abundance_high_95_raw, NA_real_),
      fractional_abundance_low_68_raw = if_else(raw_role == "mutant", fractional_abundance_low_68_raw, NA_real_),
      fractional_abundance_high_68_raw = if_else(raw_role == "mutant", fractional_abundance_high_68_raw, NA_real_)
    )
}

# ---- CSV patching ----

patch_one_csv <- function(run_row, updates) {
  csv_path <- path_from_manifest_value(run_row$csv_database_path)
  if (!file.exists(csv_path)) {
    stop("Missing CSV export: ", csv_path)
  }

  csv_data <- readr::read_csv(
    csv_path,
    col_types = readr::cols(.default = readr::col_character()),
    na = character(),
    show_col_types = FALSE
  )
  original_columns <- names(csv_data)
  conc_col <- first_matching_column(original_columns, "^Conc\\(copies/")
  copies_col <- first_matching_column(original_columns, "^Copies/20")

  indexed <- csv_data %>%
    mutate(
      .row_id = row_number(),
      .target_clean = clean_ddpcr_target(.data$Target)
    ) %>%
    left_join(
      updates %>%
        transmute(
          Well,
          target_clean,
          run_id_raw = run_id,
          raw_role,
          accepted_droplets_raw = `Accepted Droplets`,
          positives_raw = Positives,
          negatives_raw = Negatives,
          ch1_ch2_pos_raw = `Ch1+Ch2+`,
          ch1_pos_ch2_neg_raw = `Ch1+Ch2-`,
          ch1_neg_ch2_pos_raw = `Ch1-Ch2+`,
          ch1_ch2_neg_raw = `Ch1-Ch2-`,
          concentration_raw,
          copies_20ul_raw,
          linkage_raw,
          concentration_high_95_raw,
          concentration_low_95_raw,
          concentration_high_68_raw,
          concentration_low_68_raw,
          ratio_raw,
          ratio_high_95_raw,
          ratio_low_95_raw,
          ratio_high_68_raw,
          ratio_low_68_raw,
          fractional_abundance_raw,
          fractional_abundance_high_95_raw,
          fractional_abundance_low_95_raw,
          fractional_abundance_high_68_raw,
          fractional_abundance_low_68_raw,
          MeanAmplitudeOfPositives_raw,
          MeanAmplitudeOfNegatives_raw,
          MeanAmplitudeTotal_raw,
          raw_metadata_path,
          raw_peak_path
        ),
      by = c("Well" = "Well", ".target_clean" = "target_clean"),
      relationship = "one-to-one"
    )

  missing_updates <- indexed %>% filter(is.na(run_id_raw))
  if (nrow(missing_updates) > 0L) {
    stop(
      "Could not map CSV rows back to archive-derived rows for ",
      run_row$run_id, ":\n",
      paste(utils::capture.output(print(missing_updates %>% select(Well, Target, Sample))), collapse = "\n")
    )
  }

  update_specs <- list(
    list(column = conc_col, field = "concentration_raw", formatter = function(x) format_decimal(x, 8), tolerance = 1e-7),
    list(column = copies_col, field = "copies_20ul_raw", formatter = function(x) format_decimal(x, 6), tolerance = 1e-6),
    list(column = "Accepted Droplets", field = "accepted_droplets_raw", formatter = format_integer, tolerance = 0),
    list(column = "Positives", field = "positives_raw", formatter = format_integer, tolerance = 0),
    list(column = "Negatives", field = "negatives_raw", formatter = format_integer, tolerance = 0),
    list(column = "Ch1+Ch2+", field = "ch1_ch2_pos_raw", formatter = format_integer, tolerance = 0),
    list(column = "Ch1+Ch2-", field = "ch1_pos_ch2_neg_raw", formatter = format_integer, tolerance = 0),
    list(column = "Ch1-Ch2+", field = "ch1_neg_ch2_pos_raw", formatter = format_integer, tolerance = 0),
    list(column = "Ch1-Ch2-", field = "ch1_ch2_neg_raw", formatter = format_integer, tolerance = 0),
    list(column = "Linkage", field = "linkage_raw", formatter = function(x) format_decimal(x, 12), tolerance = 1e-9),
    list(column = "PoissonConfMax", field = "concentration_high_95_raw", formatter = function(x) format_decimal(x, 9), tolerance = 1e-7),
    list(column = "PoissonConfMin", field = "concentration_low_95_raw", formatter = function(x) format_decimal(x, 9), tolerance = 1e-7),
    list(column = "PoissonConfidenceMax68", field = "concentration_high_68_raw", formatter = function(x) format_decimal(x, 9), tolerance = 1e-7),
    list(column = "PoissonConfidenceMin68", field = "concentration_low_68_raw", formatter = function(x) format_decimal(x, 9), tolerance = 1e-7),
    list(column = "Ratio", field = "ratio_raw", formatter = function(x) format_decimal(x, 12), tolerance = 1e-12),
    list(column = "PoissonRatioMax", field = "ratio_high_95_raw", formatter = function(x) format_decimal(x, 12), tolerance = 1e-12),
    list(column = "PoissonRatioMin", field = "ratio_low_95_raw", formatter = function(x) format_decimal(x, 12), tolerance = 1e-12),
    list(column = "PoissonRatioMax68", field = "ratio_high_68_raw", formatter = function(x) format_decimal(x, 12), tolerance = 1e-12),
    list(column = "PoissonRatioMin68", field = "ratio_low_68_raw", formatter = function(x) format_decimal(x, 12), tolerance = 1e-12),
    list(column = "Fractional Abundance", field = "fractional_abundance_raw", formatter = function(x) format_decimal(x, 10), tolerance = 1e-8),
    list(column = "PoissonFractionalAbundanceMax", field = "fractional_abundance_high_95_raw", formatter = function(x) format_decimal(x, 10), tolerance = 1e-8),
    list(column = "PoissonFractionalAbundanceMin", field = "fractional_abundance_low_95_raw", formatter = function(x) format_decimal(x, 10), tolerance = 1e-8),
    list(column = "PoissonFractionalAbundanceMax68", field = "fractional_abundance_high_68_raw", formatter = function(x) format_decimal(x, 10), tolerance = 1e-8),
    list(column = "PoissonFractionalAbundanceMin68", field = "fractional_abundance_low_68_raw", formatter = function(x) format_decimal(x, 10), tolerance = 1e-8),
    list(column = "MeanAmplitudeOfPositives", field = "MeanAmplitudeOfPositives_raw", formatter = function(x) format_decimal(x, 12), tolerance = 1e-8),
    list(column = "MeanAmplitudeOfNegatives", field = "MeanAmplitudeOfNegatives_raw", formatter = function(x) format_decimal(x, 12), tolerance = 1e-8),
    list(column = "MeanAmplitudeTotal", field = "MeanAmplitudeTotal_raw", formatter = function(x) format_decimal(x, 12), tolerance = 1e-8)
  )

  audit_rows <- list()
  patched <- indexed
  for (spec in update_specs) {
    column <- spec$column
    field <- spec$field
    if (is.na(column) || !column %in% names(patched) || !field %in% names(patched)) {
      next
    }

    old_values <- patched[[column]]
    new_values <- replace_preserving_equivalent(
      old_values = old_values,
      raw_values = patched[[field]],
      formatter = spec$formatter,
      tolerance = spec$tolerance
    )
    status <- ifelse(old_values == new_values, "unchanged", "changed")

    audit_rows[[length(audit_rows) + 1L]] <- tibble(
      run_id = run_row$run_id,
      run_date = as.Date(run_row$run_date),
      assay = run_row$assay,
      csv_path = as.character(csv_path),
      row_id = patched$.row_id,
      well = patched$Well,
      sample = patched$Sample,
      target = patched$Target,
      target_clean = patched$.target_clean,
      column = column,
      old_value = old_values,
      new_value = new_values,
      status = status,
      raw_metadata_path = patched$raw_metadata_path,
      raw_peak_path = patched$raw_peak_path
    )

    patched[[column]] <- new_values
  }

  readr::write_csv(
    patched %>% select(all_of(original_columns)),
    csv_path,
    na = ""
  )

  bind_rows(audit_rows)
}

# ---- run amendment ----

runs <- readr::read_csv(file.path(raw_root, "manifests", "runs.csv"), show_col_types = FALSE) %>%
  filter(run_id %in% target_run_ids) %>%
  mutate(
    run_date = as.Date(run_date),
    assay = normalise_assay(assay)
  ) %>%
  arrange(run_date, assay)

missing_run_ids <- setdiff(target_run_ids, runs$run_id)
if (length(missing_run_ids) > 0L) {
  stop("Requested run IDs are absent from raw/ddpcr/manifests/runs.csv: ", paste(missing_run_ids, collapse = ", "))
}

sample_manifest <- readr::read_csv(file.path(raw_root, "manifests", "sample_manifest.csv"), show_col_types = FALSE) %>%
  filter(experiment == "SNV", assay %in% mutation_list_raw_import) %>%
  mutate(
    run_date = as.Date(run_date),
    assay = normalise_assay(assay),
    target_clean = clean_ddpcr_target(target_clean)
  )

audit <- purrr::map_dfr(seq_len(nrow(runs)), function(i) {
  run_row <- runs[i, ]
  message("Amending CSV export from archive-derived gates: ", run_row$run_id)
  updates <- build_run_updates(raw_root, run_row, sample_manifest)
  patch_one_csv(run_row, updates)
})

readr::write_csv(audit, file.path(validation_dir, "csv_amendment_manifest.csv"))

summary <- audit %>%
  count(run_id, status, name = "cell_count") %>%
  arrange(run_id, status)
readr::write_csv(summary, file.path(validation_dir, "csv_amendment_summary.csv"))

changed_summary <- audit %>%
  filter(status == "changed") %>%
  count(run_id, column, name = "changed_cells") %>%
  arrange(run_id, column)
readr::write_csv(changed_summary, file.path(validation_dir, "csv_amendment_changed_columns.csv"))

message("Wrote CSV amendment audit: ", file.path(validation_dir, "csv_amendment_manifest.csv"))
