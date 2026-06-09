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
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(positive_individual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(positive_merged_dir, recursive = TRUE, showWarnings = FALSE)

# ---- import helper functions ----
source(file.path(project_root, "src", "ddpcr", "ddpcr_raw_import_helpers.R"))

# ---- command line options ----
args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

# ---- small helpers ----

# Keep generated filenames readable while avoiding path-hostile characters.
safe_file_component <- function(x) {
  x <- str_squish(as.character(x))
  x <- str_replace_all(x, "[^A-Za-z0-9._-]+", "_")
  x <- str_replace_all(x, "_+", "_")
  str_replace_all(x, "^_+|_+$", "")
}

# Interpret workbook flags that may be logical values or text from Excel.
as_true_flag <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

# Match final result sample IDs back to raw well labels, allowing for historical
# spelling and suffix differences in the ddPCR manifests.
normalise_sample_key <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("d1789n", "d178n") %>%
    str_replace_all("(?:[_-](?:d178n|e200k|p102l)(?:[_-]?mut)?)$", "") %>%
    str_replace_all("(?:^|[_-])(?:d178n|e200k|p102l)[_-]?mut$", "") %>%
    str_replace_all("^cjd[-_]?", "") %>%
    str_replace_all("pons", "ps") %>%
    str_replace_all("substantia[_ -]*nigra", "sn") %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

gate_levels <- c(
  "Ch1+Ch2+",
  "Ch1+Ch2-",
  "Ch1-Ch2+",
  "Ch1-Ch2-",
  "Gated / unassigned",
  "Rejected / unassigned"
)

# Palette is chosen for visual separation of true quadrants and QC-only
# unassigned/rejected droplets.
gate_colours <- c(
  `Ch1+Ch2+` = "#332288",
  `Ch1+Ch2-` = "#117733",
  `Ch1-Ch2+` = "#CC6677",
  `Ch1-Ch2-` = "#6B7280",
  `Gated / unassigned` = "#DDCC77",
  `Rejected / unassigned` = "#D1D5DB"
)

# Bio-Rad metadata can store thresholds as lists or scalar values depending on
# how the well was analysed.
threshold_entry_value <- function(entry) {
  # Absent threshold entries are valid for some wells, so return a typed missing
  # value rather than failing the whole plot.
  if (is.null(entry) || length(entry) == 0L) {
    return(NA_real_)
  }

  # Some exports wrap the threshold object in a one-element list, with the
  # numeric value stored under ThresholdValue.
  if (is.list(entry) && !is.null(entry[[1]]) && is.list(entry[[1]]) &&
      !is.null(entry[[1]][["ThresholdValue"]])) {
    return(as.numeric(entry[[1]][["ThresholdValue"]]))
  }

  # Other exports store the threshold object directly at this level.
  if (is.list(entry) && !is.null(entry[["ThresholdValue"]])) {
    return(as.numeric(entry[["ThresholdValue"]]))
  }

  # Already-normalised metadata may arrive as a bare numeric vector.
  if (is.numeric(entry) && length(entry) > 0L) {
    return(as.numeric(entry[[1]]))
  }

  # Any unrecognised shape is treated as missing so the plot can still be made
  # without a dashed threshold line.
  NA_real_
}

# Prefer manual thresholds for visual QC (fall-back: auto thresholds)
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

# ---- droplet parsing ----

# Read one well from the extracted .ddpcr archive and attach the saved
# cluster calls needed for colouring the scatterplot.
read_well_droplets <- function(extract_path, run_id, well, assay) {
  peak_path <- file.path(extract_path, "PeakData", paste0(well, ".ddpeakjson"))
  metadata_path <- file.path(extract_path, "PeakMetaData", paste0(well, ".ddmetajson"))

  # Fail early when the extracted archive is incomplete for this well.
  if (!file.exists(peak_path)) {
    stop("Missing peak data: ", peak_path)
  }
  if (!file.exists(metadata_path)) {
    stop("Missing peak metadata: ", metadata_path)
  }

  peak <- jsonlite::fromJSON(peak_path, simplifyVector = FALSE)
  metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)

  # The amplitude arrays are the raw x/y coordinates for each droplet.
  amplitudes <- peak$PeakInfo$Amplitudes
  if (is.null(amplitudes) || length(amplitudes) < 2L) {
    stop("Peak data does not contain two amplitude channels: ", peak_path)
  }

  # Trim to the shorter channel defensively so malformed wells cannot misalign
  # the droplet index and amplitude columns.
  ch1 <- as.numeric(amplitudes[[1]])
  ch2 <- as.numeric(amplitudes[[2]])
  droplet_count <- min(length(ch1), length(ch2))
  droplets <- tibble(
    droplet_index = seq.int(0L, droplet_count - 1L),
    ch1_amplitude = ch1[seq_len(droplet_count)],
    ch2_amplitude = ch2[seq_len(droplet_count)],
    gate = "Rejected / unassigned"
  )

  # Select the reference and mutant targets for this assay, then map them to
  # the physical fluorescence channels used in the scatterplot axes.
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
    # Cluster droplet indices are zero-based in the archive; tibble rows are
    # one-based in R.
    droplet_indices <- as.integer(cluster$Droplets %||% integer(0))
    if (length(droplet_indices) == 0L) {
      next
    }
    row_indices <- droplet_indices + 1L
    row_indices <- row_indices[row_indices >= 1L & row_indices <= nrow(droplets)]
    if (length(row_indices) == 0L) {
      next
    }

    # Missing or unassigned target calls remain separate from rejected droplets
    # because Bio-Rad has at least put them into a cluster.
    results <- as.character(unlist(cluster$Results, use.names = FALSE))
    if (length(results) < max(selected_indices) || is_true(cluster$Unassigned)) {
      droplets$gate[row_indices] <- "Gated / unassigned"
      next
    }

    ch1_result <- results[[ch1_idx]]
    ch2_result <- results[[ch2_idx]]

    # Only the two selected assay targets are used for quadrant labels.
    if (!all(c(ch1_result, ch2_result) %in% c("Negative", "Positive"))) {
      droplets$gate[row_indices] <- "Gated / unassigned"
      next
    }

    droplets$gate[row_indices] <- paste0(
      ifelse(ch1_result == "Positive", "Ch1+Ch2", "Ch1-Ch2"),
      ifelse(ch2_result == "Positive", "+", "-")
    )
  }

  # Return droplets and QC metadata together so the plotting function writes
  # one manifest row per attempted well.
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

# ---- well plotter ----

# Render one well and return the manifest row describing either the output file
# or the reason this well could not be plotted.
plot_well <- function(well_row) {
  contributing_well <- paste(well_row$run_id, well_row$well, sep = "_")

  # Keep run, well, sample, and assay in the filename for manual review.
  output_path <- file.path(
    positive_individual_dir,
    paste0(
      safe_file_component(well_row$run_id),
      "_",
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

  # Record missing archive components in the manifest rather than aborting the
  # whole batch.
  missing_status <- case_when(
    !file.exists(peak_path) & !file.exists(metadata_path) ~ "missing_peak_data_and_metadata",
    !file.exists(peak_path) ~ "missing_peak_data",
    !file.exists(metadata_path) ~ "missing_peak_metadata",
    TRUE ~ NA_character_
  )

  if (!is.na(missing_status)) {
    return(tibble(
      plot_kind = "individual_well",
      positive_id = well_row$positive_id,
      participant = well_row$participant,
      code = well_row$code,
      brain_region = well_row$brain_region,
      mutation = well_row$mutation,
      sample_id = well_row$sample_id,
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
      metadata_path = metadata_path,
      n_wells = 1L,
      contributing_wells = contributing_well
    ))
  }

  # Existing PNGs are skipped by default so interrupted batches can resume
  # without redrawing every well.
  if (file.exists(output_path) && !force) {
    return(tibble(
      plot_kind = "individual_well",
      positive_id = well_row$positive_id,
      participant = well_row$participant,
      code = well_row$code,
      brain_region = well_row$brain_region,
      mutation = well_row$mutation,
      sample_id = well_row$sample_id,
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
      metadata_path = metadata_path,
      n_wells = 1L,
      contributing_wells = contributing_well
    ))
  }

  parsed <- read_well_droplets(extract_path, well_row$run_id, well_row$well, well_row$assay)

  # Fix legend order even when a well lacks one or more quadrants.
  droplets <- parsed$droplets %>%
    mutate(gate = factor(
      gate,
      levels = gate_levels
    ))

  ch1_threshold <- parsed$thresholds$threshold[parsed$thresholds$channel == "Ch1"][[1]]
  ch2_threshold <- parsed$thresholds$threshold[parsed$thresholds$channel == "Ch2"][[1]]

  # Plot raw droplet amplitudes, coloured by the saved cluster assignment.
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

  # Dashed thresholds are overlays for review, not recomputed gates.
  if (!is.na(ch1_threshold)) {
    plot <- plot + geom_vline(xintercept = ch1_threshold, linewidth = 0.35, linetype = "dashed")
  }
  if (!is.na(ch2_threshold)) {
    plot <- plot + geom_hline(yintercept = ch2_threshold, linewidth = 0.35, linetype = "dashed")
  }

  # Save a modest-resolution PNG for fast manual QC browsing.
  ggsave(
    filename = output_path,
    plot = plot,
    width = 6,
    height = 5,
    dpi = 180
  )

  # The manifest doubles as an audit table for droplet counts, thresholds, and
  # source files used for each scatterplot.
  tibble(
    plot_kind = "individual_well",
    positive_id = well_row$positive_id,
    participant = well_row$participant,
    code = well_row$code,
    brain_region = well_row$brain_region,
    mutation = well_row$mutation,
    sample_id = well_row$sample_id,
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
    metadata_path = metadata_path,
    n_wells = 1L,
    contributing_wells = contributing_well
  )
}

# ---- merged sample plotter ----

# Read all wells for one positive sample and pool the droplets into a single
# scatterplot.
plot_merged_sample <- function(sample_rows) {
  row <- sample_rows[1, ]
  output_path <- file.path(
    positive_merged_dir,
    paste0("merged_", row$positive_id, ".png")
  )

  parsed_wells <- purrr::map(seq_len(nrow(sample_rows)), function(i) {
    well_row <- sample_rows[i, ]
    extract_path <- file.path(raw_root, well_row$archive_contents_relative_dir)
    parsed <- read_well_droplets(extract_path, well_row$run_id, well_row$well, well_row$assay)
    list(row = well_row, parsed = parsed)
  })

  droplets <- purrr::map_dfr(parsed_wells, function(entry) {
    entry$parsed$droplets %>%
      mutate(
        gate = factor(gate, levels = gate_levels),
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

  plot <- ggplot(droplets, aes(x = ch1_amplitude, y = ch2_amplitude, colour = gate)) +
    geom_point(size = 0.18, alpha = 0.38, stroke = 0, na.rm = TRUE) +
    scale_colour_manual(values = gate_colours, drop = FALSE, name = "Gate") +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(
      title = paste(row$participant, row$brain_region, row$mutation),
      subtitle = paste0(n_distinct(droplets$well_plot_id), " well(s); ", nrow(droplets), " droplets"),
      x = "Ch1 amplitude",
      y = "Ch2 amplitude"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 9),
      plot.subtitle = element_text(size = 8)
    )

  # Draw each contributing well's exported thresholds; repeated wells can
  # legitimately have slightly different review gates.
  ch1_thresholds <- thresholds %>% filter(channel == "Ch1", !is.na(threshold))
  ch2_thresholds <- thresholds %>% filter(channel == "Ch2", !is.na(threshold))
  if (nrow(ch1_thresholds) > 0L) {
    plot <- plot +
      geom_vline(
        data = ch1_thresholds,
        aes(xintercept = threshold),
        linewidth = 0.25,
        linetype = "dashed",
        alpha = 0.55,
        inherit.aes = FALSE
      )
  }
  if (nrow(ch2_thresholds) > 0L) {
    plot <- plot +
      geom_hline(
        data = ch2_thresholds,
        aes(yintercept = threshold),
        linewidth = 0.25,
        linetype = "dashed",
        alpha = 0.55,
        inherit.aes = FALSE
      )
  }

  status <- "already_present"
  if (!file.exists(output_path) || force) {
    ggsave(
      filename = output_path,
      plot = plot,
      width = 6,
      height = 5,
      dpi = 180
    )
    status <- "written"
  }

  tibble(
    plot_kind = "merged_sample",
    positive_id = row$positive_id,
    participant = row$participant,
    code = row$code,
    brain_region = row$brain_region,
    mutation = row$mutation,
    sample_id = row$sample_id,
    output_path = output_path,
    status = status,
    n_wells = n_distinct(droplets$well_plot_id),
    droplet_count = nrow(droplets),
    contributing_wells = paste(sort(unique(droplets$well_plot_id)), collapse = ";")
  )
}

# ---- active well set ----

# Only active SNV runs are eligible for positive-sample scatterplot generation.
runs <- readr::read_csv(file.path(raw_root, "manifests", "runs.csv"), show_col_types = FALSE) %>%
  filter(status == "active", experiment == "SNV") %>%
  select(run_id, archive_contents_relative_dir)

# Build a stable well index before selecting positive samples.
active_wells <- readr::read_csv(file.path(raw_root, "manifests", "sample_manifest.csv"), show_col_types = FALSE) %>%
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

# SNV_data_final.xlsx is the source of truth for LoB/LoD positivity.
positive_samples_all <- readxl::read_excel(file.path(project_root, "results", "ddPCR", "SNV_data_final.xlsx")) %>%
  as_tibble() %>%
  mutate(
    detected_above_LoB = as_true_flag(detected_above_LoB),
    detected_above_LoD = as_true_flag(detected_above_LoD),
    is_pooled = as_true_flag(is_pooled),
    sample_id = paste(coalesce(as.character(code), as.character(participant)), brain_region, sep = "_"),
    sample_key = normalise_sample_key(sample_id),
    positive_id = safe_file_component(paste(participant, brain_region, mutation, sep = "_"))
  ) %>%
  filter(
    !is_pooled,
    detected_above_LoB,
    detected_above_LoD
  )

positive_sample_total_before_exclusion <- nrow(positive_samples_all)

# CJD30 E200K is heterozygous and is excluded from these positive-sample QC
# plots to match the manuscript gating figure workflow.
positive_samples <- positive_samples_all %>%
  filter(!(participant == "CJD30" & mutation == "E200K")) %>%
  arrange(mutation, participant, brain_region)

# Map each positive sample-region result back to its raw physical well(s).
selected_wells <- positive_samples %>%
  select(
    positive_id, participant, code, brain_region, mutation, sample_id,
    sample_key, fractional_abundance, ci_low, ci_high,
    detected_above_LoB, detected_above_LoD
  ) %>%
  inner_join(
    active_wells,
    by = c("mutation" = "assay", "sample_key" = "sample_key"),
    relationship = "many-to-many"
  ) %>%
  mutate(assay = mutation) %>%
  arrange(mutation, participant, brain_region, run_date, well)

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

# ---- render batch ----

empty_plot_manifest <- tibble(
  plot_kind = character(),
  positive_id = character(),
  participant = character(),
  code = character(),
  brain_region = character(),
  mutation = character(),
  sample_id = character(),
  well_index = integer(),
  run_id = character(),
  run_date = as.Date(character()),
  assay = character(),
  well = character(),
  sample = character(),
  output_path = character(),
  status = character(),
  droplet_count = integer(),
  accepted_clustered_droplets = integer(),
  gated_or_unassigned_droplets = integer(),
  rejected_droplet_count = integer(),
  ch1_threshold = numeric(),
  ch2_threshold = numeric(),
  peak_path = character(),
  metadata_path = character(),
  n_wells = integer(),
  contributing_wells = character()
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
      positive_id, participant, code, brain_region, mutation, sample_id,
      run_id, run_date, assay, well, sample, archive_contents_relative_dir,
      well_index, fractional_abundance, ci_low, ci_high,
      detected_above_LoB, detected_above_LoD
    ),
  file.path(positive_out_dir, "selected_wells.csv")
)
readr::write_csv(plot_manifest, file.path(positive_out_dir, "plot_manifest.csv"))

# ---- run summary ----

# Write a short text summary for quick terminal and file-based review.
summary_lines <- c(
  paste0("Active wells in full set: ", active_well_total),
  paste0("LoB+/LoD+ sample rows before CJD30 E200K exclusion: ", positive_sample_total_before_exclusion),
  paste0("LoB+/LoD+ sample rows after CJD30 E200K exclusion: ", nrow(positive_samples)),
  paste0("Selected positive wells: ", nrow(selected_wells)),
  paste0("Individual well scatterplots written this run: ", sum(plot_manifest$plot_kind == "individual_well" & plot_manifest$status == "written")),
  paste0("Merged sample scatterplots written this run: ", sum(plot_manifest$plot_kind == "merged_sample" & plot_manifest$status == "written")),
  paste0("Plot manifest rows: ", nrow(plot_manifest)),
  paste0("Scatterplot files present in manifest: ", sum(!is.na(plot_manifest$output_path) & file.exists(plot_manifest$output_path)))
)
writeLines(summary_lines, file.path(positive_out_dir, "scatterplot_summary.txt"))
cat(paste(summary_lines, collapse = "\n"))
cat("\n")
