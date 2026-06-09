library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)
library(jsonlite)
library(ggplot2)

# ---- relative paths ----
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
lod_raw_root <- file.path(project_root, "raw", "ddpcr_lod")
out_root <- file.path(project_root, "results", "ddPCR", "scatterplots")
lod_out_dir <- file.path(out_root, "lod_experiments")
lod_individual_dir <- file.path(lod_out_dir, "individual_wells")
dir.create(lod_individual_dir, recursive = TRUE, showWarnings = FALSE)

# ---- import helper functions ----
source(file.path(project_root, "src", "ddpcr", "ddpcr_raw_import_helpers.R"))
source(file.path(project_root, "src", "ddpcr", "ddpcr_scatterplot_helpers.R"))

# ---- LoD raw database ----

read_lod_assay_rows <- function(assay_dir) {
  manifests_dir <- file.path(assay_dir, "manifests")
  run_path <- file.path(manifests_dir, "run.csv")
  sample_manifest_path <- file.path(manifests_dir, "sample_manifest.csv")

  if (!file.exists(run_path)) {
    stop("Missing LoD run manifest: ", run_path)
  }
  if (!file.exists(sample_manifest_path)) {
    stop("Missing LoD sample manifest: ", sample_manifest_path)
  }

  run_row <- readr::read_csv(run_path, show_col_types = FALSE)
  if (nrow(run_row) != 1L) {
    stop("Expected one run row in ", run_path, ", found ", nrow(run_row))
  }

  assay <- normalise_assay(run_row$assay[[1]])
  archive_contents_dir <- file.path(assay_dir, "archive_contents")
  layout_path <- list.files(file.path(assay_dir, "layout_xlsx"), pattern = "\\.xlsx$", full.names = TRUE)
  analysis_path <- list.files(file.path(assay_dir, "analysis_csv"), pattern = "\\.csv$", full.names = TRUE)

  readr::read_csv(sample_manifest_path, show_col_types = FALSE) %>%
    mutate(
      assay = normalise_assay(assay),
      target_normalised = clean_ddpcr_target(target_normalised),
      run_date = as.Date(run_row$run_date[[1]]),
      archive_contents_dir = archive_contents_dir,
      source_layout = if (length(layout_path) > 0L) layout_path[[1]] else NA_character_,
      source_analysis_csv = if (length(analysis_path) > 0L) analysis_path[[1]] else NA_character_,
      source_ddpcr = file.path(assay_dir, "ddpcr_archive", basename(run_row$ddpcr_archive[[1]])),
      original_source_layout = run_row$source_layout_xlsx[[1]],
      original_source_analysis_csv = run_row$source_analysis_csv[[1]]
    ) %>%
    filter(target_normalised == assay) %>%
    separate_rows(merged_wells, sep = ",\\s*") %>%
    mutate(
      physical_well = str_squish(merged_wells),
      lod_merged_well = well,
      lod_sample_label = sample,
      lod_fraction = as.character(sample_normalised)
    ) %>%
    filter(physical_well != "")
}

assay_dirs <- list.dirs(lod_raw_root, recursive = FALSE, full.names = TRUE)
assay_dirs <- assay_dirs[file.exists(file.path(assay_dirs, "manifests", "sample_manifest.csv"))]
if (length(assay_dirs) == 0L) {
  stop("No LoD assay manifests found under ", lod_raw_root)
}

selected_wells <- assay_dirs %>%
  sort() %>%
  purrr::map_dfr(read_lod_assay_rows) %>%
  arrange(run_date, assay, lod_merged_well, physical_well) %>%
  mutate(well_index = row_number())

if (nrow(selected_wells) == 0L) {
  stop("No LoD physical wells were selected from ", lod_raw_root)
}

compute_lod_axis_limits <- function(run_rows) {
  droplets <- purrr::map_dfr(seq_len(nrow(run_rows)), function(i) {
    well_row <- run_rows[i, ]
    tryCatch(
      read_well_droplets(
        extract_path = well_row$archive_contents_dir,
        run_id = well_row$run_id,
        well = well_row$physical_well,
        assay = well_row$assay
      )$droplets,
      error = function(e) tibble(x_amplitude = numeric(), y_amplitude = numeric(), gate = character())
    )
  }) %>%
    filter(gate != rejected_gate_level)

  x_limits <- axis_limits_for_values(droplets$x_amplitude)
  y_limits <- axis_limits_for_values(droplets$y_amplitude)
  tibble(
    run_id = run_rows$run_id[[1]],
    x_min = x_limits[[1]],
    x_max = x_limits[[2]],
    y_min = y_limits[[1]],
    y_max = y_limits[[2]]
  )
}

run_axis_limits <- selected_wells %>%
  group_by(run_id) %>%
  group_split() %>%
  purrr::map_dfr(compute_lod_axis_limits)

selected_wells <- selected_wells %>%
  left_join(run_axis_limits, by = "run_id")

# ---- well plotter ----

plot_lod_well <- function(well_row) {
  contributing_well <- paste(well_row$run_id, well_row$physical_well, sep = "_")
  file_stem <- paste(
    safe_file_component(as.character(well_row$run_date)),
    safe_file_component(well_row$assay),
    well_row$physical_well,
    safe_file_component(well_row$lod_sample_label),
    sep = "_"
  )
  output_paths <- plot_output_paths(lod_individual_dir, file_stem)

  peak_path <- file.path(well_row$archive_contents_dir, "PeakData", paste0(well_row$physical_well, ".ddpeakjson"))
  metadata_path <- file.path(well_row$archive_contents_dir, "PeakMetaData", paste0(well_row$physical_well, ".ddmetajson"))

  missing_status <- case_when(
    !file.exists(peak_path) & !file.exists(metadata_path) ~ "missing_peak_data_and_metadata",
    !file.exists(peak_path) ~ "missing_peak_data",
    !file.exists(metadata_path) ~ "missing_peak_metadata",
    TRUE ~ NA_character_
  )

  if (!is.na(missing_status)) {
    return(bind_cols(
      tibble(
        plot_kind = "lod_individual_well",
        lod_id = safe_file_component(contributing_well),
        run_id = well_row$run_id,
        run_date = as.Date(well_row$run_date),
        assay = well_row$assay,
        lod_merged_well = well_row$lod_merged_well,
        physical_well = well_row$physical_well,
        sample = well_row$lod_sample_label,
        lod_fraction = well_row$lod_fraction,
        accepted_droplets = well_row$accepted_droplets,
        positives = well_row$positives,
        negatives = well_row$negatives,
        fractional_abundance = well_row$fractional_abundance,
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
        source_layout = well_row$source_layout,
        source_analysis_csv = well_row$source_analysis_csv,
        n_wells = 1L,
        contributing_wells = contributing_well
      ),
      plot_path_columns(output_paths)
    ))
  }

  parsed <- read_well_droplets(
    extract_path = well_row$archive_contents_dir,
    run_id = well_row$run_id,
    well = well_row$physical_well,
    assay = well_row$assay
  )
  droplets <- parsed$droplets %>%
    mutate(gate = factor(gate, levels = c(plot_gate_levels, rejected_gate_level)))

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
    mutation = well_row$assay,
    title = paste(well_row$run_date, well_row$assay, well_row$physical_well, well_row$lod_sample_label),
    subtitle = paste0("LoD fraction ", well_row$lod_fraction, "; merged sample ", well_row$lod_merged_well),
    axis_limits = axis_limits
  )

  save_plot_outputs(plot, output_paths, width = 6, height = 5, dpi = 220)

  bind_cols(
    tibble(
      plot_kind = "lod_individual_well",
      lod_id = safe_file_component(contributing_well),
      run_id = well_row$run_id,
      run_date = as.Date(well_row$run_date),
      assay = well_row$assay,
      lod_merged_well = well_row$lod_merged_well,
      physical_well = well_row$physical_well,
      sample = well_row$lod_sample_label,
      lod_fraction = well_row$lod_fraction,
      accepted_droplets = well_row$accepted_droplets,
      positives = well_row$positives,
      negatives = well_row$negatives,
      fractional_abundance = well_row$fractional_abundance,
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
      source_layout = well_row$source_layout,
      source_analysis_csv = well_row$source_analysis_csv,
      n_wells = 1L,
      contributing_wells = contributing_well
    ),
    plot_path_columns(output_paths)
  )
}

# ---- render batch ----

clean_plot_output_dir(lod_individual_dir)

plot_manifest <- purrr::map_dfr(seq_len(nrow(selected_wells)), function(i) {
  plot_lod_well(selected_wells[i, ])
}) %>%
  arrange(run_date, assay, lod_merged_well, physical_well)

# ---- manifest update ----

readr::write_csv(
  selected_wells %>%
    select(
      run_id, run_date, assay, lod_merged_well, physical_well,
      sample = lod_sample_label, lod_fraction, accepted_droplets,
      positives, negatives, fractional_abundance, source_layout,
      source_analysis_csv, archive_contents_dir, well_index,
      x_min, x_max, y_min, y_max
    ),
  file.path(lod_out_dir, "selected_wells.csv")
)
readr::write_csv(plot_manifest, file.path(lod_out_dir, "plot_manifest.csv"))

# ---- run summary ----

plot_file_paths <- unlist(plot_manifest %>% select(output_png_path, output_svg_path, output_pdf_path), use.names = FALSE)
known_lod_runs <- c(
  "2020-09-11_SNV_D178N",
  "2020-10-08_SNV_E200K",
  "2021-01-25_SNV_P102L"
)
known_lod_runs_present <- intersect(known_lod_runs, unique(selected_wells$run_id))
known_lod_summary <- if (length(known_lod_runs_present) == 0L) {
  "none"
} else {
  paste(known_lod_runs_present, collapse = ";")
}

summary_lines <- c(
  paste0("LoD raw root: ", lod_raw_root),
  paste0("LoD assay runs selected: ", n_distinct(selected_wells$run_id)),
  paste0("Selected LoD physical wells: ", nrow(selected_wells)),
  paste0("Individual LoD scatterplots written this run: ", sum(plot_manifest$status == "written")),
  paste0("Plot manifest rows: ", nrow(plot_manifest)),
  paste0("Scatterplot files present in manifest: ", sum(!is.na(plot_file_paths) & file.exists(plot_file_paths))),
  paste0("Known LoD runs present: ", known_lod_summary)
)
writeLines(summary_lines, file.path(lod_out_dir, "scatterplot_summary.txt"))
cat(paste(summary_lines, collapse = "\n"))
cat("\n")
