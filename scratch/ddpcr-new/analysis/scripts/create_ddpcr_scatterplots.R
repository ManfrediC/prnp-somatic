library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(jsonlite)
library(ggplot2)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
scratch_root <- file.path(project_root, "scratch", "ddpcr-new")
analysis_root <- file.path(scratch_root, "analysis")
raw_root <- file.path(scratch_root, "raw")
validation_dir <- file.path(analysis_root, "validation")
out_root <- file.path(scratch_root, "figures", "scatterplots")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

source(file.path(analysis_root, "scripts", "ddpcr_raw_import_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
limit_arg <- args[startsWith(args, "--limit=")]
limit <- if (length(limit_arg) == 1L) {
  as.integer(sub("^--limit=", "", limit_arg))
} else {
  NA_integer_
}
start_arg <- args[startsWith(args, "--start=")]
start <- if (length(start_arg) == 1L) {
  as.integer(sub("^--start=", "", start_arg))
} else {
  1L
}

safe_file_component <- function(x) {
  x <- str_squish(as.character(x))
  x <- str_replace_all(x, "[^A-Za-z0-9._-]+", "_")
  x <- str_replace_all(x, "_+", "_")
  str_replace_all(x, "^_+|_+$", "")
}

threshold_entry_value <- function(entry) {
  if (is.null(entry) || length(entry) == 0L) {
    return(NA_real_)
  }
  if (is.list(entry) && !is.null(entry[[1]]) && is.list(entry[[1]]) &&
      !is.null(entry[[1]][["ThresholdValue"]])) {
    return(as.numeric(entry[[1]][["ThresholdValue"]]))
  }
  if (is.list(entry) && !is.null(entry[["ThresholdValue"]])) {
    return(as.numeric(entry[["ThresholdValue"]]))
  }
  if (is.numeric(entry) && length(entry) > 0L) {
    return(as.numeric(entry[[1]]))
  }
  NA_real_
}

threshold_for_channel <- function(metadata, channel) {
  threshold <- threshold_entry_value((metadata$ThresholdValues %||% list())[[channel]])
  threshold_type <- "manual"
  if (is.na(threshold)) {
    threshold <- threshold_entry_value((metadata$AutoThresholdValues %||% list())[[channel]])
    threshold_type <- "auto"
  }
  tibble(
    channel = paste0("Ch", channel),
    threshold = threshold,
    threshold_type = if_else(is.na(threshold), NA_character_, threshold_type)
  )
}

read_well_droplets <- function(extract_path, run_id, well, assay) {
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
  amplitudes <- peak$PeakInfo$Amplitudes
  if (is.null(amplitudes) || length(amplitudes) < 2L) {
    stop("Peak data does not contain two amplitude channels: ", peak_path)
  }

  ch1 <- as.numeric(amplitudes[[1]])
  ch2 <- as.numeric(amplitudes[[2]])
  droplet_count <- min(length(ch1), length(ch2))
  droplets <- tibble(
    droplet_index = seq.int(0L, droplet_count - 1L),
    ch1_amplitude = ch1[seq_len(droplet_count)],
    ch2_amplitude = ch2[seq_len(droplet_count)],
    gate = "Rejected / unassigned"
  )

  targets <- metadata_targets(metadata)
  selected <- selected_target_indices(targets, assay)
  if (is.null(targets) || is.null(selected)) {
    stop("Could not select targets for scatterplot: ", run_id, " ", well)
  }

  channels <- vapply(targets, target_channel, integer(1))
  selected_indices <- c(selected$ref, selected$mut)
  ch1_idx <- selected_indices[channels[selected_indices] == 1L][1]
  ch2_idx <- selected_indices[channels[selected_indices] == 2L][1]
  if (is.na(ch1_idx) || is.na(ch2_idx)) {
    stop("Could not map channels for scatterplot: ", run_id, " ", well)
  }

  # Gate colours come from the cluster assignments saved in the .ddpcr file.
  # Thresholds are plotted below for visual QC only; they do not reclassify
  # droplets in this script.
  for (cluster in metadata$Clusters %||% list()) {
    droplet_indices <- as.integer(cluster$Droplets %||% integer(0))
    if (length(droplet_indices) == 0L) {
      next
    }
    row_indices <- droplet_indices + 1L
    row_indices <- row_indices[row_indices >= 1L & row_indices <= nrow(droplets)]
    if (length(row_indices) == 0L) {
      next
    }

    results <- as.character(unlist(cluster$Results, use.names = FALSE))
    if (length(results) < max(selected_indices) || is_true(cluster$Unassigned)) {
      droplets$gate[row_indices] <- "Gated / unassigned"
      next
    }

    ch1_result <- results[[ch1_idx]]
    ch2_result <- results[[ch2_idx]]
    if (!all(c(ch1_result, ch2_result) %in% c("Negative", "Positive"))) {
      droplets$gate[row_indices] <- "Gated / unassigned"
      next
    }

    droplets$gate[row_indices] <- paste0(
      ifelse(ch1_result == "Positive", "Ch1+Ch2", "Ch1-Ch2"),
      ifelse(ch2_result == "Positive", "+", "-")
    )
  }

  thresholds <- bind_rows(
    threshold_for_channel(metadata, 1L),
    threshold_for_channel(metadata, 2L)
  )

  list(
    droplets = droplets,
    thresholds = thresholds,
    peak_path = peak_path,
    metadata_path = metadata_path,
    rejected_droplet_count = as.integer(peak$RejectedInfo$RejectedDropletCount %||% NA_integer_)
  )
}

plot_well <- function(well_row) {
  run_dir <- file.path(out_root, as.character(well_row$run_date), safe_file_component(well_row$run_id))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  output_path <- file.path(
    run_dir,
    paste0(
      well_row$well,
      "_",
      safe_file_component(well_row$sample),
      "_",
      safe_file_component(well_row$assay),
      ".png"
    )
  )
  extract_path <- file.path(raw_root, well_row$archive_contents_relative_dir)
  peak_path <- file.path(extract_path, "PeakData", paste0(well_row$well, ".ddpeakjson"))
  metadata_path <- file.path(extract_path, "PeakMetaData", paste0(well_row$well, ".ddmetajson"))

  missing_status <- case_when(
    !file.exists(peak_path) & !file.exists(metadata_path) ~ "missing_peak_data_and_metadata",
    !file.exists(peak_path) ~ "missing_peak_data",
    !file.exists(metadata_path) ~ "missing_peak_metadata",
    TRUE ~ NA_character_
  )

  if (!is.na(missing_status)) {
    return(tibble(
      well_index = well_row$well_index,
      run_id = well_row$run_id,
      run_date = as.Date(well_row$run_date),
      assay = well_row$assay,
      well = well_row$well,
      sample = well_row$sample,
      output_path = NA_character_,
      status = missing_status,
      droplet_count = NA_integer_,
      accepted_clustered_droplets = NA_integer_,
      gated_or_unassigned_droplets = NA_integer_,
      rejected_droplet_count = NA_integer_,
      ch1_threshold = NA_real_,
      ch2_threshold = NA_real_,
      peak_path = peak_path,
      metadata_path = metadata_path
    ))
  }

  if (file.exists(output_path) && !force) {
    return(tibble(
      well_index = well_row$well_index,
      run_id = well_row$run_id,
      run_date = as.Date(well_row$run_date),
      assay = well_row$assay,
      well = well_row$well,
      sample = well_row$sample,
      output_path = output_path,
      status = "already_present",
      droplet_count = NA_integer_,
      accepted_clustered_droplets = NA_integer_,
      gated_or_unassigned_droplets = NA_integer_,
      rejected_droplet_count = NA_integer_,
      ch1_threshold = NA_real_,
      ch2_threshold = NA_real_,
      peak_path = peak_path,
      metadata_path = metadata_path
    ))
  }

  parsed <- read_well_droplets(extract_path, well_row$run_id, well_row$well, well_row$assay)
  droplets <- parsed$droplets %>%
    mutate(gate = factor(
      gate,
      levels = c(
        "Ch1+Ch2+",
        "Ch1+Ch2-",
        "Ch1-Ch2+",
        "Ch1-Ch2-",
        "Gated / unassigned",
        "Rejected / unassigned"
      )
    ))

  gate_colours <- c(
    `Ch1+Ch2+` = "#332288",
    `Ch1+Ch2-` = "#117733",
    `Ch1-Ch2+` = "#CC6677",
    `Ch1-Ch2-` = "#6B7280",
    `Gated / unassigned` = "#DDCC77",
    `Rejected / unassigned` = "#D1D5DB"
  )

  ch1_threshold <- parsed$thresholds$threshold[parsed$thresholds$channel == "Ch1"][[1]]
  ch2_threshold <- parsed$thresholds$threshold[parsed$thresholds$channel == "Ch2"][[1]]

  plot <- ggplot(droplets, aes(x = ch1_amplitude, y = ch2_amplitude, colour = gate)) +
    geom_point(size = 0.25, alpha = 0.45, stroke = 0, na.rm = TRUE) +
    scale_colour_manual(values = gate_colours, drop = FALSE, name = "Gate") +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(
      title = paste(well_row$run_id, well_row$well, well_row$sample),
      x = "Ch1 amplitude",
      y = "Ch2 amplitude"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 9)
    )

  if (!is.na(ch1_threshold)) {
    plot <- plot + geom_vline(xintercept = ch1_threshold, linewidth = 0.35, linetype = "dashed")
  }
  if (!is.na(ch2_threshold)) {
    plot <- plot + geom_hline(yintercept = ch2_threshold, linewidth = 0.35, linetype = "dashed")
  }

  ggsave(
    filename = output_path,
    plot = plot,
    width = 6,
    height = 5,
    dpi = 180
  )

  tibble(
    well_index = well_row$well_index,
    run_id = well_row$run_id,
    run_date = as.Date(well_row$run_date),
    assay = well_row$assay,
    well = well_row$well,
    sample = well_row$sample,
    output_path = output_path,
    status = "written",
    droplet_count = nrow(droplets),
    accepted_clustered_droplets = sum(droplets$gate %in% c("Ch1+Ch2+", "Ch1+Ch2-", "Ch1-Ch2+", "Ch1-Ch2-")),
    gated_or_unassigned_droplets = sum(droplets$gate == "Gated / unassigned"),
    rejected_droplet_count = parsed$rejected_droplet_count,
    ch1_threshold = ch1_threshold,
    ch2_threshold = ch2_threshold,
    peak_path = peak_path,
    metadata_path = metadata_path
  )
}

runs <- readr::read_csv(file.path(raw_root, "manifests", "runs.csv"), show_col_types = FALSE) %>%
  filter(status == "active", experiment == "SNV") %>%
  select(run_id, archive_contents_relative_dir)

active_wells <- readr::read_csv(file.path(raw_root, "manifests", "sample_manifest.csv"), show_col_types = FALSE) %>%
  filter(experiment == "SNV", assay %in% mutation_list_raw_import) %>%
  distinct(run_id, run_date, assay, well, sample) %>%
  left_join(runs, by = "run_id") %>%
  arrange(run_date, assay, well) %>%
  mutate(well_index = row_number())

active_well_total <- nrow(active_wells)

if (is.na(start) || start < 1L) {
  start <- 1L
}

if (!is.na(limit)) {
  stop_index <- min(nrow(active_wells), start + limit - 1L)
  active_wells <- active_wells %>%
    slice(seq.int(start, stop_index))
} else if (start > 1L) {
  active_wells <- active_wells %>%
    slice(seq.int(start, nrow(active_wells)))
}

current_manifest <- purrr::map_dfr(seq_len(nrow(active_wells)), function(i) {
  plot_well(active_wells[i, ])
})

scatterplot_manifest_path <- file.path(validation_dir, "scatterplot_manifest.csv")
if (file.exists(scatterplot_manifest_path)) {
  previous_manifest <- readr::read_csv(scatterplot_manifest_path, show_col_types = FALSE)
} else {
  previous_manifest <- tibble()
}

scatterplot_manifest <- bind_rows(previous_manifest, current_manifest) %>%
  arrange(well_index) %>%
  distinct(run_id, well, .keep_all = TRUE)

readr::write_csv(scatterplot_manifest, scatterplot_manifest_path)

summary_lines <- c(
  paste0("Active wells in full set: ", active_well_total),
  paste0("Active wells considered in this run: ", nrow(active_wells)),
  paste0("Scatterplots written this run: ", sum(current_manifest$status == "written")),
  paste0("Scatterplots already present this run: ", sum(current_manifest$status == "already_present")),
  paste0("Active wells missing peak data this run: ", sum(current_manifest$status == "missing_peak_data")),
  paste0("Scatterplot manifest rows: ", nrow(scatterplot_manifest)),
  paste0("Active wells missing peak data in manifest: ", sum(scatterplot_manifest$status == "missing_peak_data")),
  paste0("Scatterplot files present in manifest: ", sum(!is.na(scatterplot_manifest$output_path) & file.exists(scatterplot_manifest$output_path)))
)
writeLines(summary_lines, file.path(validation_dir, "scatterplot_summary.txt"))
cat(paste(summary_lines, collapse = "\n"))
cat("\n")
