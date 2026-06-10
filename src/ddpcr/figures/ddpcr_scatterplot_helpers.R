library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(jsonlite)
library(ggplot2)
library(grid)

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

called_gate_levels <- c(
  "wt_pos_mut_pos",
  "wt_pos_mut_neg",
  "wt_neg_mut_pos",
  "wt_neg_mut_neg"
)
plot_gate_levels <- called_gate_levels
rejected_gate_level <- "rejected_unassigned"

# Okabe-Ito colours keep the called populations readable for colourblind readers.
gate_colours <- c(
  wt_pos_mut_pos = "#E69F00",
  wt_pos_mut_neg = "#009E73",
  wt_neg_mut_pos = "#CC79A7",
  wt_neg_mut_neg = "#5684E9"
)

# Build manuscript-facing legend labels for the active mutation.
gate_labels_for_mutation <- function(mutation) {
  c(
    wt_pos_mut_pos = paste0("Wildtype+ / ", mutation, "+"),
    wt_pos_mut_neg = paste0("Wildtype+ / ", mutation, "-"),
    wt_neg_mut_pos = paste0("Wildtype- / ", mutation, "+"),
    wt_neg_mut_neg = paste0("Wildtype- / ", mutation, "-")
  )
}

# Convert final-workbook or sample-detail names into display labels.
participant_display_label <- function(x) {
  as.character(x) %>%
    str_replace("^Control", "Ctrl")
}

# Expand common region abbreviations used in raw ddPCR labels.
brain_region_display_label <- function(x) {
  x <- str_to_lower(as.character(x))
  recode(
    x,
    ps = "pons",
    pons = "pons",
    th = "thalamus",
    cb = "cerebellum",
    fr = "frontal_cortex",
    sn = "substantia_nigra",
    .default = x
  )
}

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

# Read only source/manual thresholds from the .ddpcr metadata. Auto thresholds
# can be display artefacts in clustered wells, so they are not used for these
# manuscript-facing gate overlays.
manual_threshold_for_channel <- function(metadata, channel) {
  threshold <- threshold_entry_value((metadata$ThresholdValues %||% list())[[channel]])
  tibble(
    channel = paste0("Ch", channel),
    threshold = threshold,
    threshold_source = if_else(is.na(threshold), NA_character_, "manual_exported")
  )
}

# Estimate a visual threshold from QuantaSoft-labelled clusters when no manual
# source threshold is exported for that axis.
infer_cluster_boundary <- function(droplets, amplitude_col, call_col) {
  values <- droplets[[amplitude_col]]
  calls <- droplets[[call_col]]
  usable <- is.finite(values) & calls %in% c("Negative", "Positive")
  values <- values[usable]
  calls <- calls[usable]

  negative <- values[calls == "Negative"]
  positive <- values[calls == "Positive"]
  if (length(negative) == 0L || length(positive) == 0L) {
    return(NA_real_)
  }

  positive_high <- stats::median(positive, na.rm = TRUE) > stats::median(negative, na.rm = TRUE)
  ordered <- tibble(value = values, call = calls) %>%
    arrange(value) %>%
    mutate(
      next_value = lead(value),
      next_call = lead(call),
      gap = next_value - value
    ) %>%
    filter(is.finite(gap), gap >= 0)

  # Prefer the visible transition between saved QuantaSoft classes. This avoids
  # drawing a line through a sparse population's tail when positives are rare.
  transition <- if (positive_high) {
    ordered %>% filter(call == "Negative", next_call == "Positive")
  } else {
    ordered %>% filter(call == "Positive", next_call == "Negative")
  }
  if (nrow(transition) > 0L) {
    row <- transition[which.max(transition$gap), ]
    return((row$value + row$next_value) / 2)
  }

  if (positive_high) {
    low_edge <- as.numeric(stats::quantile(negative, 0.99, names = FALSE, na.rm = TRUE))
    high_edge <- as.numeric(stats::quantile(positive, 0.01, names = FALSE, na.rm = TRUE))
  } else {
    low_edge <- as.numeric(stats::quantile(positive, 0.99, names = FALSE, na.rm = TRUE))
    high_edge <- as.numeric(stats::quantile(negative, 0.01, names = FALSE, na.rm = TRUE))
  }
  if (is.finite(low_edge) && is.finite(high_edge) && high_edge > low_edge) {
    return((low_edge + high_edge) / 2)
  }

  (stats::median(negative, na.rm = TRUE) + stats::median(positive, na.rm = TRUE)) / 2
}

# Choose the gate line for one plotted axis: source/manual threshold first,
# otherwise an inferred boundary from the saved cluster labels.
threshold_for_axis <- function(metadata, channel, axis, droplets, amplitude_col, call_col) {
  manual <- manual_threshold_for_channel(metadata, channel)
  threshold <- manual$threshold[[1]]
  threshold_source <- manual$threshold_source[[1]]
  if (is.na(threshold)) {
    threshold <- infer_cluster_boundary(droplets, amplitude_col, call_col)
    threshold_source <- if_else(is.na(threshold), NA_character_, "cluster_boundary_inferred")
  }

  tibble(
    axis = axis,
    channel = paste0("Ch", channel),
    threshold = threshold,
    threshold_source = threshold_source
  )
}

# Keep audit strings compact and stable in plot_manifest.csv.
threshold_details <- function(thresholds) {
  details <- thresholds %>%
    filter(!is.na(threshold)) %>%
    mutate(detail = paste0(
      if ("well_plot_id" %in% names(.)) paste0(well_plot_id, " ") else "",
      axis, "=", round(threshold, 1), " ",
      channel, " ", threshold_source
    )) %>%
    pull(detail)
  if (length(details) == 0L) {
    return(NA_character_)
  }
  paste(details, collapse = ";")
}

# Scale each sample's plots together, with a little headroom for dots and gate
# lines while keeping zero anchored when all amplitudes are positive.
axis_limits_for_values <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(c(NA_real_, NA_real_))
  }
  lower <- min(0, floor(min(values) / 500) * 500)
  upper <- ceiling((max(values) * 1.03) / 500) * 500
  if (!is.finite(upper) || upper <= lower) {
    upper <- lower + 500
  }
  c(lower, upper)
}

# Build the three output paths that every scatterplot writer records in its
# manifest row.
plot_output_paths <- function(output_dir, file_stem) {
  c(
    png = file.path(output_dir, paste0(file_stem, ".png")),
    svg = file.path(output_dir, paste0(file_stem, ".svg")),
    pdf = file.path(output_dir, paste0(file_stem, ".pdf"))
  )
}

# Save each plot in the raster/vector formats used for review and manuscripts.
save_plot_outputs <- function(plot, paths, width = 6, height = 5, dpi = 180) {
  ggsave(filename = paths[["png"]], plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(filename = paths[["svg"]], plot = plot, width = width, height = height, device = grDevices::svg, bg = "transparent")
  ggsave(filename = paths[["pdf"]], plot = plot, width = width, height = height, device = grDevices::pdf, bg = "transparent")
}

# Mirror the generated paths into the manifest column layout.
plot_path_columns <- function(paths) {
  tibble(
    output_path = paths[["png"]],
    output_png_path = paths[["png"]],
    output_svg_path = paths[["svg"]],
    output_pdf_path = paths[["pdf"]]
  )
}

# Remove prior rendered plot files before a rerun writes a fresh manifest.
clean_plot_output_dir <- function(output_dir) {
  files <- list.files(output_dir, pattern = "\\.(png|svg|pdf)$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) > 0L) {
    unlink(files)
  }
}

# Replace guide grobs with blanks while preserving their layout columns. This
# keeps plot geometry stable for the downstream SVG panel assembly.
without_legend_keep_layout <- function(plot) {
  grob <- ggplotGrob(plot)
  guide_indexes <- grep("^guide-box", grob$layout$name)
  for (guide_index in guide_indexes) {
    grob$grobs[[guide_index]] <- grid::nullGrob()
  }
  grob
}

# Keep only the ggplot guide in its original right-side layout slot so the
# panel builder can place a single legend without extracting it from a plot.
legend_only_keep_layout <- function(plot) {
  grob <- ggplotGrob(plot)
  guide_indexes <- grep("^guide-box", grob$layout$name)
  if (length(guide_indexes) == 0L) {
    stop("Could not find a ggplot guide box for the scatterplot legend")
  }

  for (index in seq_along(grob$grobs)) {
    if (!index %in% guide_indexes) {
      grob$grobs[[index]] <- grid::nullGrob()
    }
  }
  grob
}

# Prefer consolidated manifests when present, but support reviewer workspaces
# that only have the immutable per-run manifest exports.
read_ddpcr_manifest_table <- function(raw_root, filename) {
  manifest_path <- file.path(raw_root, "manifests", filename)
  if (file.exists(manifest_path)) {
    return(readr::read_csv(manifest_path, show_col_types = FALSE))
  }

  by_run_dir <- file.path(raw_root, "manifests", "_by_run")
  by_run_paths <- list.files(by_run_dir, recursive = TRUE, full.names = TRUE)
  by_run_paths <- by_run_paths[basename(by_run_paths) == filename]
  by_run_paths <- sort(by_run_paths)
  if (length(by_run_paths) == 0L) {
    stop("Could not find ddPCR manifest ", filename, " under ", file.path(raw_root, "manifests"))
  }

  purrr::map_dfr(by_run_paths, readr::read_csv, show_col_types = FALSE) %>%
    distinct()
}

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
    ref_call = NA_character_,
    mut_call = NA_character_,
    gate = rejected_gate_level
  )

  # Select the reference and mutant targets for this assay, then map them to
  # the physical fluorescence channels used in the scatterplot axes.
  targets <- metadata_targets(metadata)
  selected <- selected_target_indices(targets, assay)
  if (is.null(targets) || is.null(selected)) {
    stop("Could not select targets for scatterplot: ", run_id, " ", well)
  }

  channels <- vapply(targets, target_channel, integer(1))
  ref_channel <- as.integer(channels[[selected$ref]])
  mut_channel <- as.integer(channels[[selected$mut]])
  selected_indices <- c(selected$ref, selected$mut)
  if (!all(c(ref_channel, mut_channel) %in% c(1L, 2L)) || ref_channel == mut_channel) {
    stop("Could not map channels for scatterplot: ", run_id, " ", well)
  }

  droplets <- droplets %>%
    mutate(
      x_amplitude = if (mut_channel == 1L) ch1_amplitude else ch2_amplitude,
      y_amplitude = if (ref_channel == 1L) ch1_amplitude else ch2_amplitude
    )

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
      droplets$gate[row_indices] <- "gated_unassigned"
      next
    }

    ref_result <- results[[selected$ref]]
    mut_result <- results[[selected$mut]]

    # Only the two selected assay targets are used for quadrant labels.
    if (!all(c(ref_result, mut_result) %in% c("Negative", "Positive"))) {
      droplets$gate[row_indices] <- "gated_unassigned"
      next
    }

    droplets$ref_call[row_indices] <- ref_result
    droplets$mut_call[row_indices] <- mut_result
    droplets$gate[row_indices] <- paste0(
      ifelse(ref_result == "Positive", "wt_pos", "wt_neg"),
      "_",
      ifelse(mut_result == "Positive", "mut_pos", "mut_neg")
    )
  }

  # Return droplets and QC metadata together so the plotting function writes
  # one manifest row per attempted well.
  thresholds <- bind_rows(
    threshold_for_axis(metadata, mut_channel, "x", droplets, "x_amplitude", "mut_call"),
    threshold_for_axis(metadata, ref_channel, "y", droplets, "y_amplitude", "ref_call")
  )

  list(
    droplets = droplets,
    thresholds = thresholds,
    ref_channel = paste0("Ch", ref_channel),
    mut_channel = paste0("Ch", mut_channel),
    peak_path = peak_path,
    metadata_path = metadata_path,
    rejected_droplet_count = as.integer(peak$RejectedInfo$RejectedDropletCount %||% NA_integer_)
  )
}

# Build one manuscript-facing scatterplot from already-parsed droplets.
build_scatterplot <- function(
  droplets,
  thresholds,
  mutation,
  title,
  subtitle = NULL,
  axis_limits = NULL,
  called_point_size = 1,
  called_alpha = 0.9,
  unassigned_point_size = 0.5,
  unassigned_alpha = 0.3
) {
  plot_droplets <- droplets %>%
    filter(gate %in% called_gate_levels) %>%
    mutate(gate = factor(gate, levels = plot_gate_levels))

  plot <- ggplot() +
    geom_point(
      data = plot_droplets,
      aes(x = x_amplitude, y = y_amplitude, colour = gate),
      size = called_point_size,
      alpha = called_alpha,
      stroke = 0,
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = gate_colours,
      breaks = plot_gate_levels,
      labels = gate_labels_for_mutation(mutation),
      drop = FALSE,
      name = "Droplet class"
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2.7, alpha = 1))) +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste0(mutation, " amplitude"),
      y = "Wildtype amplitude"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8)
    )

  x_thresholds <- thresholds %>% filter(axis == "x", !is.na(threshold))
  y_thresholds <- thresholds %>% filter(axis == "y", !is.na(threshold))
  if (nrow(x_thresholds) > 0L) {
    plot <- plot +
      geom_vline(
        data = x_thresholds,
        aes(xintercept = threshold),
        linewidth = 0.35,
        linetype = "dotted",
        alpha = 0.85
      )
  }
  if (nrow(y_thresholds) > 0L) {
    plot <- plot +
      geom_hline(
        data = y_thresholds,
        aes(yintercept = threshold),
        linewidth = 0.35,
        linetype = "dotted",
        alpha = 0.85
      )
  }

  if (!is.null(axis_limits) &&
      all(!is.na(c(axis_limits$x_min, axis_limits$x_max, axis_limits$y_min, axis_limits$y_max)))) {
    plot <- plot +
      coord_cartesian(
        xlim = c(axis_limits$x_min, axis_limits$x_max),
        ylim = c(axis_limits$y_min, axis_limits$y_max),
        expand = FALSE
      )
  }

  plot
}
