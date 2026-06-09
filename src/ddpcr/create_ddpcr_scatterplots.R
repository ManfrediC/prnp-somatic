library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(jsonlite)
library(ggplot2)

# ---- relative paths ----
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
raw_root <- file.path(project_root, "raw", "ddpcr")
out_root <- file.path(project_root, "results", "ddPCR", "scatterplots")
positive_out_dir <- file.path(out_root, "lob_lod_positive")
positive_individual_dir <- file.path(positive_out_dir, "individual_wells")
positive_merged_dir <- file.path(positive_out_dir, "merged_samples")
legacy_pooled_dir <- file.path(positive_out_dir, "pooled_samples")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(positive_individual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(positive_merged_dir, recursive = TRUE, showWarnings = FALSE)

# ---- import helper functions ----
source(file.path(project_root, "src", "ddpcr", "ddpcr_raw_import_helpers.R"))

# ---- shared scatterplot helpers ----
source(file.path(project_root, "src", "ddpcr", "ddpcr_scatterplot_helpers.R"))

# ---- well plotter ----

# Render one well and return the manifest row describing either the output file
# or the reason this well could not be plotted.
plot_well <- function(well_row) {
  contributing_well <- paste(well_row$run_id, well_row$well, sep = "_")

  # Keep only manuscript-facing identifiers in the filename.
  file_stem <- paste(
    safe_file_component(as.character(well_row$run_date)),
    well_row$well,
    safe_file_component(well_row$participant_display),
    safe_file_component(well_row$brain_region_display),
    safe_file_component(well_row$mutation),
    sep = "_"
  )
  output_paths <- plot_output_paths(positive_individual_dir, file_stem)

  extract_path <- file.path(raw_root, well_row$archive_contents_relative_dir)
  peak_path <- file.path(extract_path, "PeakData", paste0(well_row$well, ".ddpeakjson"))
  metadata_path <- file.path(extract_path, "PeakMetaData", paste0(well_row$well, ".ddmetajson"))

  # Record missing archive components in the manifest rather than aborting the
  # whole batch.
  missing_status <- case_when(
    !file.exists(peak_path) & !file.exists(metadata_path) ~ "missing_peak_data_and_metadata",
    !file.exists(peak_path) ~ "missing_peak_data",
    !file.exists(metadata_path) ~ "missing_peak_metadata",
    TRUE ~ NA_character_
  )

  if (!is.na(missing_status)) {
    return(bind_cols(
      tibble(
        plot_kind = "individual_well",
        positive_id = well_row$positive_id,
        participant = well_row$participant,
        participant_display = well_row$participant_display,
        code = well_row$code,
        brain_region = well_row$brain_region,
        brain_region_display = well_row$brain_region_display,
        mutation = well_row$mutation,
        replicate_index = well_row$replicate_index,
        replicate_label = well_row$replicate_label,
        sample_id = well_row$sample_id,
        sample_key = well_row$sample_key,
        is_pooled = well_row$is_pooled,
        well_index = well_row$well_index,
        run_id = well_row$run_id,
        run_date = as.Date(well_row$run_date),
        assay = well_row$assay,
        well = well_row$well,
        sample = well_row$sample,
        status = missing_status,
        droplet_count = NA_integer_,
        plotted_droplet_count = NA_integer_,
        accepted_clustered_droplets = NA_integer_,
        gated_or_unassigned_droplets = NA_integer_,
        rejected_or_unassigned_droplets = NA_integer_,
        rejected_droplet_count = NA_integer_,
        x_threshold = NA_real_,
        y_threshold = NA_real_,
        x_threshold_source = NA_character_,
        y_threshold_source = NA_character_,
        threshold_details = NA_character_,
        peak_path = peak_path,
        metadata_path = metadata_path,
        n_wells = 1L,
        contributing_wells = contributing_well
      ),
      plot_path_columns(output_paths)
    ))
  }

  parsed <- read_well_droplets(extract_path, well_row$run_id, well_row$well, well_row$assay)

  # Fix legend order even when a well lacks one or more quadrants.
  droplets <- parsed$droplets %>%
    mutate(gate = factor(
      gate,
      levels = c(plot_gate_levels, rejected_gate_level)
    ))

  x_threshold <- parsed$thresholds$threshold[parsed$thresholds$axis == "x"][[1]]
  y_threshold <- parsed$thresholds$threshold[parsed$thresholds$axis == "y"][[1]]
  x_threshold_source <- parsed$thresholds$threshold_source[parsed$thresholds$axis == "x"][[1]]
  y_threshold_source <- parsed$thresholds$threshold_source[parsed$thresholds$axis == "y"][[1]]

  axis_limits <- tibble(
    x_min = well_row$x_min,
    x_max = well_row$x_max,
    y_min = well_row$y_min,
    y_max = well_row$y_max
  )

  plot <- build_scatterplot(
    droplets = droplets,
    thresholds = parsed$thresholds,
    mutation = well_row$mutation,
    title = paste(well_row$participant_display, well_row$brain_region_display, well_row$mutation, well_row$replicate_label),
    axis_limits = axis_limits
  )

  save_plot_outputs(plot, output_paths, width = 6, height = 5, dpi = 220)

  # The manifest doubles as an audit table for droplet counts, thresholds, and
  # source files used for each scatterplot.
  bind_cols(
    tibble(
      plot_kind = "individual_well",
      positive_id = well_row$positive_id,
      participant = well_row$participant,
      participant_display = well_row$participant_display,
      code = well_row$code,
      brain_region = well_row$brain_region,
      brain_region_display = well_row$brain_region_display,
      mutation = well_row$mutation,
      replicate_index = well_row$replicate_index,
      replicate_label = well_row$replicate_label,
      sample_id = well_row$sample_id,
      sample_key = well_row$sample_key,
      is_pooled = well_row$is_pooled,
      well_index = well_row$well_index,
      run_id = well_row$run_id,
      run_date = as.Date(well_row$run_date),
      assay = well_row$assay,
      well = well_row$well,
      sample = well_row$sample,
      status = "written",
      droplet_count = nrow(droplets),
      plotted_droplet_count = sum(droplets$gate %in% called_gate_levels),
      accepted_clustered_droplets = sum(droplets$gate %in% called_gate_levels),
      gated_or_unassigned_droplets = sum(droplets$gate == "gated_unassigned"),
      rejected_or_unassigned_droplets = sum(droplets$gate == rejected_gate_level),
      rejected_droplet_count = parsed$rejected_droplet_count,
      x_threshold = x_threshold,
      y_threshold = y_threshold,
      x_threshold_source = x_threshold_source,
      y_threshold_source = y_threshold_source,
      threshold_details = threshold_details(parsed$thresholds),
      peak_path = peak_path,
      metadata_path = metadata_path,
      n_wells = 1L,
      contributing_wells = contributing_well
    ),
    plot_path_columns(output_paths)
  )
}

# ---- merged sample plotter ----

# Read all wells for one positive sample and pool the droplets into a single
# scatterplot.
plot_merged_sample <- function(sample_rows) {
  row <- sample_rows[1, ]
  output_paths <- plot_output_paths(positive_merged_dir, paste0("merged_", row$positive_id))

  parsed_wells <- purrr::map(seq_len(nrow(sample_rows)), function(i) {
    well_row <- sample_rows[i, ]
    extract_path <- file.path(raw_root, well_row$archive_contents_relative_dir)
    parsed <- read_well_droplets(extract_path, well_row$run_id, well_row$well, well_row$assay)
    list(row = well_row, parsed = parsed)
  })

  droplets <- purrr::map_dfr(parsed_wells, function(entry) {
    entry$parsed$droplets %>%
      mutate(
        gate = factor(gate, levels = c(plot_gate_levels, rejected_gate_level)),
        well_plot_id = paste(entry$row$run_id, entry$row$well, sep = "_"),
        run_id = entry$row$run_id,
        well = entry$row$well,
        sample = entry$row$sample
      )
  })

  thresholds <- purrr::map_dfr(parsed_wells, function(entry) {
    entry$parsed$thresholds %>%
      mutate(
        well_plot_id = paste(entry$row$run_id, entry$row$well, sep = "_"),
        run_id = entry$row$run_id,
        well = entry$row$well,
        sample = entry$row$sample
      )
  })

  axis_limits <- tibble(
    x_min = row$x_min,
    x_max = row$x_max,
    y_min = row$y_min,
    y_max = row$y_max
  )

  plot <- build_scatterplot(
    droplets = droplets,
    thresholds = thresholds,
    mutation = row$mutation,
    title = paste(row$participant_display, row$brain_region_display, row$mutation),
    subtitle = paste0(n_distinct(droplets$well_plot_id), " well(s); ", sum(droplets$gate %in% called_gate_levels), " plotted droplets"),
    axis_limits = axis_limits,
    called_point_size = 0.30,
    called_alpha = 0.55,
    unassigned_point_size = 0.24,
    unassigned_alpha = 0.12
  )

  save_plot_outputs(plot, output_paths, width = 6, height = 5, dpi = 220)
  status <- "written"

  bind_cols(
    tibble(
      plot_kind = "merged_sample",
      positive_id = row$positive_id,
      participant = row$participant,
      participant_display = row$participant_display,
      code = row$code,
      brain_region = row$brain_region,
      brain_region_display = row$brain_region_display,
      mutation = row$mutation,
      replicate_index = NA_integer_,
      replicate_label = NA_character_,
      sample_id = row$sample_id,
      sample_key = row$sample_key,
      is_pooled = row$is_pooled,
      well_index = NA_integer_,
      run_id = NA_character_,
      run_date = as.Date(NA),
      assay = row$assay,
      well = NA_character_,
      sample = NA_character_,
      status = status,
      droplet_count = nrow(droplets),
      plotted_droplet_count = sum(droplets$gate %in% called_gate_levels),
      accepted_clustered_droplets = sum(droplets$gate %in% called_gate_levels),
      gated_or_unassigned_droplets = sum(droplets$gate == "gated_unassigned"),
      rejected_or_unassigned_droplets = sum(droplets$gate == rejected_gate_level),
      rejected_droplet_count = sum(purrr::map_int(parsed_wells, function(entry) entry$parsed$rejected_droplet_count %||% 0L)),
      x_threshold = NA_real_,
      y_threshold = NA_real_,
      x_threshold_source = "multiple_wells",
      y_threshold_source = "multiple_wells",
      threshold_details = threshold_details(thresholds),
      peak_path = NA_character_,
      metadata_path = NA_character_,
      n_wells = n_distinct(droplets$well_plot_id),
      contributing_wells = paste(sort(unique(droplets$well_plot_id)), collapse = ";")
    ),
    plot_path_columns(output_paths)
  )
}

# ---- active well set ----

# Only active SNV runs are eligible for positive-sample scatterplot generation.
runs <- read_ddpcr_manifest_table(raw_root, "runs.csv") %>%
  filter(status == "active", experiment == "SNV") %>%
  select(run_id, archive_contents_relative_dir)

# Build a stable well index before selecting positive samples.
active_wells <- read_ddpcr_manifest_table(raw_root, "sample_manifest.csv") %>%
  filter(experiment == "SNV", assay %in% mutation_list_raw_import) %>%
  distinct(run_id, run_date, assay, well, sample) %>%
  left_join(runs, by = "run_id") %>%
  arrange(run_date, assay, well) %>%
  mutate(
    run_date = as.Date(run_date),
    sample_key = normalise_sample_key(sample),
    well_index = row_number()
  )

active_well_total <- nrow(active_wells)

# ---- positive sample selection ----

sample_details <- readxl::read_excel(file.path(raw_root, "sample_details.xlsx")) %>%
  as_tibble() %>%
  transmute(
    code = as.character(code),
    participant_display = participant_display_label(new_name)
  )

# SNV_data_final.xlsx is the source of truth for LoB/LoD positivity.
positive_samples_all <- readxl::read_excel(file.path(project_root, "results", "ddPCR", "SNV_data_final.xlsx")) %>%
  as_tibble() %>%
  mutate(
    participant = as.character(participant),
    code = as.character(code),
    detected_above_LoB = as_true_flag(detected_above_LoB),
    detected_above_LoD = as_true_flag(detected_above_LoD),
    is_pooled = as_true_flag(is_pooled),
    sample_id = paste(coalesce(as.character(code), as.character(participant)), brain_region, sep = "_"),
    sample_key = normalise_sample_key(sample_id)
  ) %>%
  left_join(sample_details, by = "code") %>%
  mutate(
    participant_display = coalesce(participant_display, participant_display_label(participant)),
    brain_region_display = brain_region_display_label(brain_region),
    positive_id = safe_file_component(paste(participant_display, brain_region_display, mutation, sep = "_"))
  ) %>%
  filter(
    detected_above_LoB,
    detected_above_LoD
  )

positive_sample_total_before_exclusion <- nrow(positive_samples_all)

# CJD30 E200K is heterozygous and is excluded from these positive-sample QC plots
positive_samples <- positive_samples_all %>%
  filter(!(participant == "CJD30" & mutation == "E200K")) %>%
  arrange(mutation, participant, brain_region)

# Map each positive sample-region result back to its raw physical well(s).
selected_wells <- positive_samples %>%
  select(
    positive_id, participant, participant_display, code, brain_region,
    brain_region_display, mutation, sample_id,
    sample_key, is_pooled, fractional_abundance, ci_low, ci_high,
    detected_above_LoB, detected_above_LoD
  ) %>%
  inner_join(
    active_wells,
    by = c("mutation" = "assay", "sample_key" = "sample_key"),
    relationship = "many-to-many"
  ) %>%
  mutate(assay = mutation) %>%
  arrange(mutation, participant_display, brain_region_display, run_date, well)

unmapped_positive_samples <- positive_samples %>%
  anti_join(
    selected_wells %>% distinct(positive_id),
    by = "positive_id"
  )
if (nrow(unmapped_positive_samples) > 0L) {
  stop(
    "Could not map LoB+/LoD+ sample rows back to raw wells:\n",
    paste(utils::capture.output(print(unmapped_positive_samples)), collapse = "\n")
  )
}

replicate_dates <- selected_wells %>%
  distinct(positive_id, run_date) %>%
  arrange(positive_id, run_date) %>%
  group_by(positive_id) %>%
  mutate(
    replicate_index = row_number(),
    replicate_label = paste("replicate", replicate_index)
  ) %>%
  ungroup()

selected_wells <- selected_wells %>%
  left_join(replicate_dates, by = c("positive_id", "run_date"))

compute_sample_axis_limits <- function(sample_rows) {
  droplets <- purrr::map_dfr(seq_len(nrow(sample_rows)), function(i) {
    well_row <- sample_rows[i, ]
    extract_path <- file.path(raw_root, well_row$archive_contents_relative_dir)
    tryCatch(
      read_well_droplets(extract_path, well_row$run_id, well_row$well, well_row$assay)$droplets,
      error = function(e) tibble(x_amplitude = numeric(), y_amplitude = numeric(), gate = character())
    )
  }) %>%
    filter(gate != rejected_gate_level)

  x_limits <- axis_limits_for_values(droplets$x_amplitude)
  y_limits <- axis_limits_for_values(droplets$y_amplitude)
  tibble(
    positive_id = sample_rows$positive_id[[1]],
    x_min = x_limits[[1]],
    x_max = x_limits[[2]],
    y_min = y_limits[[1]],
    y_max = y_limits[[2]]
  )
}

sample_axis_limits <- if (nrow(selected_wells) > 0L) {
  selected_wells %>%
    group_by(positive_id) %>%
    group_split() %>%
    purrr::map_dfr(compute_sample_axis_limits)
} else {
  tibble(positive_id = character(), x_min = numeric(), x_max = numeric(), y_min = numeric(), y_max = numeric())
}

selected_wells <- selected_wells %>%
  left_join(sample_axis_limits, by = "positive_id")

# ---- render batch ----

if (dir.exists(legacy_pooled_dir)) {
  legacy_pooled_dir_norm <- normalizePath(legacy_pooled_dir, winslash = "/", mustWork = TRUE)
  expected_legacy_dir_norm <- normalizePath(
    file.path(positive_out_dir, "pooled_samples"),
    winslash = "/",
    mustWork = TRUE
  )
  if (legacy_pooled_dir_norm != expected_legacy_dir_norm) {
    stop("Refusing to remove unexpected legacy pooled output directory: ", legacy_pooled_dir_norm)
  }
  unlink(legacy_pooled_dir_norm, recursive = TRUE, force = TRUE)
  if (dir.exists(legacy_pooled_dir_norm)) {
    stop("Failed to remove legacy pooled output directory: ", legacy_pooled_dir_norm)
  }
}

# These are derived QC outputs; clear stale graphics so each run reflects the
# current code and data.
clean_plot_output_dir(positive_individual_dir)
clean_plot_output_dir(positive_merged_dir)

empty_plot_manifest <- tibble(
  plot_kind = character(),
  positive_id = character(),
  participant = character(),
  participant_display = character(),
  code = character(),
  brain_region = character(),
  brain_region_display = character(),
  mutation = character(),
  replicate_index = integer(),
  replicate_label = character(),
  sample_id = character(),
  sample_key = character(),
  is_pooled = logical(),
  well_index = integer(),
  run_id = character(),
  run_date = as.Date(character()),
  assay = character(),
  well = character(),
  sample = character(),
  status = character(),
  droplet_count = integer(),
  plotted_droplet_count = integer(),
  accepted_clustered_droplets = integer(),
  gated_or_unassigned_droplets = integer(),
  rejected_or_unassigned_droplets = integer(),
  rejected_droplet_count = integer(),
  x_threshold = numeric(),
  y_threshold = numeric(),
  x_threshold_source = character(),
  y_threshold_source = character(),
  threshold_details = character(),
  peak_path = character(),
  metadata_path = character(),
  n_wells = integer(),
  contributing_wells = character(),
  output_path = character(),
  output_png_path = character(),
  output_svg_path = character(),
  output_pdf_path = character()
)

# Plot each selected physical well.
individual_manifest <- if (nrow(selected_wells) > 0L) {
  purrr::map_dfr(seq_len(nrow(selected_wells)), function(i) {
    plot_well(selected_wells[i, ])
  })
} else {
  empty_plot_manifest
}

# Pool repeat wells for each positive sample-region result.
merged_manifest <- if (nrow(selected_wells) > 0L) {
  selected_wells %>%
    group_by(positive_id) %>%
    group_split() %>%
    purrr::map_dfr(plot_merged_sample)
} else {
  empty_plot_manifest
}

plot_manifest <- bind_rows(individual_manifest, merged_manifest) %>%
  arrange(positive_id, plot_kind, run_date, well)

# ---- manifest update ----

readr::write_csv(
  selected_wells %>%
    select(
      positive_id, participant, participant_display, code, brain_region,
      brain_region_display, mutation, replicate_index, replicate_label, sample_id,
      sample_key, is_pooled,
      run_id, run_date, assay, well, sample, archive_contents_relative_dir,
      well_index, fractional_abundance, ci_low, ci_high,
      detected_above_LoB, detected_above_LoD,
      x_min, x_max, y_min, y_max
    ),
  file.path(positive_out_dir, "selected_wells.csv")
)
readr::write_csv(plot_manifest, file.path(positive_out_dir, "plot_manifest.csv"))

# ---- run summary ----

plot_file_paths <- unlist(plot_manifest %>% select(output_png_path, output_svg_path, output_pdf_path), use.names = FALSE)

# Write a short text summary for quick terminal and file-based review.
summary_lines <- c(
  paste0("Active wells in full set: ", active_well_total),
  paste0("LoB+/LoD+ sample rows before CJD30 E200K exclusion: ", positive_sample_total_before_exclusion),
  paste0("LoB+/LoD+ sample rows after CJD30 E200K exclusion: ", nrow(positive_samples)),
  paste0("Selected positive wells: ", nrow(selected_wells)),
  paste0("Individual well scatterplots written this run: ", sum(plot_manifest$plot_kind == "individual_well" & plot_manifest$status == "written")),
  paste0("Merged sample scatterplots written this run: ", sum(plot_manifest$plot_kind == "merged_sample" & plot_manifest$status == "written")),
  paste0("Plot manifest rows: ", nrow(plot_manifest)),
  paste0("Scatterplot files present in manifest: ", sum(!is.na(plot_file_paths) & file.exists(plot_file_paths)))
)
writeLines(summary_lines, file.path(positive_out_dir, "scatterplot_summary.txt"))
cat(paste(summary_lines, collapse = "\n"))
cat("\n")
