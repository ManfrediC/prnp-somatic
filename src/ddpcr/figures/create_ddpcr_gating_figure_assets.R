library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)
library(jsonlite)
library(ggplot2)
library(openxlsx)

# ---- paths and output directories ----

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
raw_root <- file.path(project_root, "raw", "ddpcr")
positive_out_dir <- file.path(project_root, "manuscript", "figures", "ddpcr_gating_lob_lod_positive")
strategy_out_dir <- file.path(project_root, "manuscript", "figures", "ddpcr_gating_strategy")
positive_individual_dir <- file.path(positive_out_dir, "individual")
strategy_individual_dir <- file.path(strategy_out_dir, "individual")
dir.create(positive_individual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(strategy_individual_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "src", "ddpcr", "ddpcr_raw_import_helpers.R"))

# ---- plotting constants ----

# Keep mutation, region, and control-stage order fixed across manifests and
# rendered panels.
mutation_order <- c("D178N", "E200K", "P102L")
region_labels <- c(
  bg = "basal ganglia",
  bs = "brainstem",
  cb = "cerebellum",
  fr = "frontal cortex",
  hc = "hippocampus",
  th = "thalamus"
)
control_stage_order <- c("NTC", "WT", "Positive control", "Final adjusted gate")
gate_threshold_tolerance <- 1e-6
relaxed_gate_threshold_tolerance <- 100

class_levels <- c(
  "Reference-only",
  "Mutant-only",
  "Double-positive",
  "Double-negative",
  "Gated/unassigned",
  "Rejected/unassigned"
)
plot_class_levels <- setdiff(class_levels, c("Gated/unassigned", "Rejected/unassigned"))
class_colours <- c(
  `Reference-only` = "#009E73",
  `Mutant-only` = "#CC79A7",
  `Double-positive` = "#E69F00",
  `Double-negative` = "#5684E9"
)
draw_order <- c(
  "Rejected/unassigned" = 1L,
  "Gated/unassigned" = 2L,
  "Double-negative" = 3L,
  "Reference-only" = 4L,
  "Mutant-only" = 5L,
  "Double-positive" = 6L
)

# ---- small normalisation helpers ----

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
    str_replace_all("pons", "bs") %>%
    str_replace_all("substantia[_ -]*nigra", "bs") %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

# Format fractional-abundance values for compact plot subtitles.
format_pct <- function(x) {
  format(round(as.numeric(x), 3), nsmall = 3, trim = TRUE, scientific = FALSE)
}

# ---- threshold and class helpers ----

# Estimate the visible boundary between saved QuantaSoft negative and positive
# calls for one fluorescence channel.
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

cluster_boundary_threshold_for_channel <- function(droplets, channel, amplitude_col, call_col) {
  threshold <- infer_cluster_boundary(droplets, amplitude_col, call_col)
  tibble(
    channel = channel,
    threshold_type = "cluster_boundary_inferred",
    threshold = threshold
  ) %>%
    filter(!is.na(threshold), is.finite(threshold), threshold > 0)
}

cluster_boundary_thresholds <- function(droplets) {
  bind_rows(
    cluster_boundary_threshold_for_channel(droplets, "Ch1", "ch1_amplitude", "ch1_call"),
    cluster_boundary_threshold_for_channel(droplets, "Ch2", "ch2_amplitude", "ch2_call")
  )
}

# User-reviewed gates for edited QuantaSoft archives. These are deliberately
# scoped to the exact runs whose .ddpcr files were corrected, so all other wells
# retain the historical cluster-boundary fallback.
provided_run_thresholds <- tribble(
  ~run_id, ~assay, ~physical_ch1_threshold, ~physical_ch2_threshold,
  "2020-11-24_SNV_D178N", "D178N", 1000, 2250,
  "2020-11-05_SNV_E200K", "E200K", 1600, 2450,
  "2021-02-17_SNV_P102L", "P102L", 2240, 2528
)

provided_final_gate_wells <- tribble(
  ~run_id, ~assay, ~well,
  "2020-11-24_SNV_D178N", "D178N", "B07",
  "2020-11-05_SNV_E200K", "E200K", "D06",
  "2021-02-17_SNV_P102L", "P102L", "E06"
)

strategy_axis_overrides <- tribble(
  ~assay, ~x_max, ~y_max,
  "D178N", 6000, 7000,
  "E200K", 6000, 10000,
  "P102L", NA_real_, 7500
)

provided_thresholds_for_well <- function(well_row, mut_channel, ref_channel) {
  provided <- provided_run_thresholds %>%
    filter(run_id == well_row$run_id, assay == well_row$assay)

  if (nrow(provided) != 1L) {
    return(NULL)
  }

  physical_thresholds <- c(
    Ch1 = provided$physical_ch1_threshold,
    Ch2 = provided$physical_ch2_threshold
  )

  # The provided gate coordinates are physical archive Ch1/Ch2 thresholds. The
  # plot axes are semantic: mutant on x and WT/reference on y. Keep the legacy
  # channel labels because downstream plot code treats Ch1 as x and Ch2 as y.
  tibble(
    channel = c("Ch1", "Ch2"),
    threshold_type = "manual",
    threshold = c(
      physical_thresholds[[paste0("Ch", mut_channel)]],
      physical_thresholds[[paste0("Ch", ref_channel)]]
    )
  )
}

# Convert QuantaSoft reference/mutant calls into manuscript droplet classes.
class_from_calls <- function(ref_result, mut_result) {
  if (!all(c(ref_result, mut_result) %in% c("Negative", "Positive"))) {
    return("Gated/unassigned")
  }
  if (ref_result == "Positive" && mut_result == "Positive") {
    return("Double-positive")
  }
  if (ref_result == "Positive" && mut_result == "Negative") {
    return("Reference-only")
  }
  if (ref_result == "Negative" && mut_result == "Positive") {
    return("Mutant-only")
  }
  "Double-negative"
}

# ---- raw well parsing ----

# Read one extracted .ddpcr well and return plotted droplets, thresholds, and
# provenance paths for manifests.
read_well_droplets_for_plot <- function(well_row) {
  extract_path <- file.path(raw_root, well_row$archive_contents_relative_dir)
  peak_path <- file.path(extract_path, "PeakData", paste0(well_row$well, ".ddpeakjson"))
  metadata_path <- file.path(extract_path, "PeakMetaData", paste0(well_row$well, ".ddmetajson"))
  if (!file.exists(peak_path)) {
    stop("Missing peak data: ", peak_path)
  }
  if (!file.exists(metadata_path)) {
    stop("Missing peak metadata: ", metadata_path)
  }

  peak <- jsonlite::fromJSON(peak_path, simplifyVector = FALSE)
  metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)

  # The amplitude arrays are the raw physical-channel coordinates for each
  # droplet. They are remapped below onto semantic plot axes: mutant on x and
  # WT/reference on y.
  amplitudes <- peak$PeakInfo$Amplitudes
  if (is.null(amplitudes) || length(amplitudes) < 2L) {
    stop("Peak data does not contain two amplitude channels: ", peak_path)
  }

  # Select the reference and mutant targets for this assay, then map them to
  # physical fluorescence channels.
  targets <- metadata_targets(metadata)
  selected <- selected_target_indices(targets, well_row$assay)
  if (is.null(targets) || is.null(selected)) {
    stop("Could not select WT and mutant targets for ", well_row$run_id, " ", well_row$well)
  }

  target_names <- vapply(targets, target_name, character(1))
  channels <- vapply(targets, target_channel, integer(1))
  selected_indices <- c(selected$ref, selected$mut)
  ref_channel <- channels[[selected$ref]]
  mut_channel <- channels[[selected$mut]]
  if (
    is.na(ref_channel) || is.na(mut_channel) ||
      !ref_channel %in% seq_along(amplitudes) ||
      !mut_channel %in% seq_along(amplitudes)
  ) {
    stop("Could not map selected WT and mutant targets to amplitude channels for ", well_row$run_id, " ", well_row$well)
  }

  # Start all droplets as rejected/unassigned, then fill class labels from the
  # saved cluster assignments.
  mut_amplitude <- as.numeric(amplitudes[[mut_channel]])
  ref_amplitude <- as.numeric(amplitudes[[ref_channel]])
  droplet_count <- min(length(mut_amplitude), length(ref_amplitude))
  droplets <- tibble(
    droplet_index = seq.int(0L, droplet_count - 1L),
    ch1_amplitude = mut_amplitude[seq_len(droplet_count)],
    ch2_amplitude = ref_amplitude[seq_len(droplet_count)],
    ch1_call = NA_character_,
    ch2_call = NA_character_,
    droplet_class = "Rejected/unassigned"
  )

  # Cluster droplet indices are zero-based in the archive; tibble rows are
  # one-based in R.
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
      droplets$droplet_class[row_indices] <- "Gated/unassigned"
      droplets$ch1_call[row_indices] <- "Unassigned"
      droplets$ch2_call[row_indices] <- "Unassigned"
      next
    }

    ref_result <- results[[selected$ref]]
    mut_result <- results[[selected$mut]]
    droplets$ch1_call[row_indices] <- mut_result
    droplets$ch2_call[row_indices] <- ref_result
    droplets$droplet_class[row_indices] <- class_from_calls(ref_result, mut_result)
  }

  # Return droplets and threshold metadata together so plot and manifest writers
  # use the same parsed source.
  thresholds <- provided_thresholds_for_well(well_row, mut_channel = mut_channel, ref_channel = ref_channel)
  if (is.null(thresholds)) {
    thresholds <- cluster_boundary_thresholds(droplets)
  }
  thresholds <- thresholds %>%
    mutate(
      run_id = well_row$run_id,
      run_date = as.Date(well_row$run_date),
      assay = well_row$assay,
      well = well_row$well,
      sample = well_row$sample
    )

  list(
    droplets = droplets %>%
      mutate(
        droplet_class = factor(droplet_class, levels = class_levels),
        draw_order = unname(draw_order[as.character(droplet_class)]),
        run_id = well_row$run_id,
        run_date = as.Date(well_row$run_date),
        assay = well_row$assay,
        well = well_row$well,
        sample = well_row$sample,
        well_label = paste0(as.Date(well_row$run_date), " ", well_row$well)
      ),
    thresholds = thresholds,
    peak_path = peak_path,
    metadata_path = metadata_path,
    rejected_droplet_count = as.integer(peak$RejectedInfo$RejectedDropletCount %||% NA_integer_),
    ref_target = target_names[[selected$ref]],
    mut_target = target_names[[selected$mut]]
  )
}

# ---- plot assembly helpers ----

# Choose common axes with enough headroom for dense droplet clouds and gate lines.
axis_limits <- function(droplets, thresholds, droplet_quantile = 0.997) {
  x_candidate <- max(
    quantile(droplets$ch1_amplitude, droplet_quantile, na.rm = TRUE, names = FALSE),
    thresholds$threshold[thresholds$channel == "Ch1"],
    na.rm = TRUE
  )
  y_candidate <- max(
    quantile(droplets$ch2_amplitude, droplet_quantile, na.rm = TRUE, names = FALSE),
    thresholds$threshold[thresholds$channel == "Ch2"],
    na.rm = TRUE
  )
  list(
    x = c(0, x_candidate * 1.08),
    y = c(0, y_candidate * 1.08)
  )
}

# Add cluster-boundary gate lines inferred from saved QuantaSoft calls.
add_gate_layers <- function(
  plot,
  thresholds,
  limits,
  show_gates = TRUE,
  show_auto = TRUE,
  show_manual = TRUE
) {
  threshold_data <- thresholds %>%
    filter(!is.na(threshold), is.finite(threshold))
  if (!show_gates || nrow(threshold_data) == 0L) {
    return(plot)
  }

  # Channel 1 is the x-axis, and Channel 2 is the y-axis in these gating plots.
  vline_data <- threshold_data %>% filter(channel == "Ch1")
  hline_data <- threshold_data %>% filter(channel == "Ch2")
  if (nrow(vline_data) > 0L) {
    plot <- plot +
      geom_vline(
        data = vline_data,
        aes(xintercept = threshold),
        linewidth = 0.35,
        colour = "black",
        linetype = "dotted",
        alpha = 0.85,
        show.legend = FALSE
      )
  }
  if (nrow(hline_data) > 0L) {
    plot <- plot +
      geom_hline(
        data = hline_data,
        aes(yintercept = threshold),
        linewidth = 0.35,
        colour = "black",
        linetype = "dotted",
        alpha = 0.85,
        show.legend = FALSE
      )
  }
  plot
}

# Build one droplet scatterplot with optional faceting and gate overlays.
base_scatter_plot <- function(
  droplets,
  thresholds,
  title,
  subtitle,
  mutation = NULL,
  limits,
  facet_by_well = FALSE,
  show_auto = TRUE,
  show_manual = TRUE,
  show_gates = TRUE,
  show_legend = TRUE,
  plot_title_size = 9
) {
  x_label <- if (!is.null(mutation) && isTRUE(nzchar(mutation))) paste0(mutation, " amplitude") else "Channel 1 amplitude"
  y_label <- if (!is.null(mutation) && isTRUE(nzchar(mutation))) "Wildtype amplitude" else "Channel 2 amplitude"

  # Plot only called droplet clusters while preserving unassigned/rejected
  # droplets in parsed metadata for audit outputs.
  plot <- droplets %>%
    filter(droplet_class %in% plot_class_levels) %>%
    arrange(draw_order, droplet_index) %>%
    ggplot(aes(x = ch1_amplitude, y = ch2_amplitude, colour = droplet_class)) +
    geom_point(size = 0.18, alpha = 0.42, stroke = 0, na.rm = TRUE) +
    scale_colour_manual(
      values = class_colours,
      breaks = plot_class_levels,
      drop = FALSE,
      name = "Droplet class"
    ) +
    coord_cartesian(xlim = limits$x, ylim = limits$y, expand = FALSE) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    ) +
    guides(
      colour = if (show_legend) guide_legend(override.aes = list(size = 2.0, alpha = 1)) else "none",
      linetype = "none"
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.box = "vertical",
      plot.title = element_text(size = plot_title_size, face = "bold"),
      plot.subtitle = element_text(size = 7),
      strip.background = element_rect(fill = "#F3F4F6", colour = "#D1D5DB"),
      strip.text = element_text(size = 7)
    )

  # Faceted plots keep contributing wells visible for merged positive samples.
  if (facet_by_well) {
    plot <- plot + facet_wrap(~well_label, ncol = 1)
  }

  plot <- add_gate_layers(
    plot = plot,
    thresholds = thresholds,
    limits = limits,
    show_gates = show_gates,
    show_auto = show_auto,
    show_manual = show_manual
  )
  plot
}

is_real_control_sample <- function(sample_upper) {
  str_detect(sample_upper, "^(?:CJD|CTRL|CONTROL)")
}

select_final_gate_control <- function(
  assay_name,
  run_id_value,
  positive_thresholds,
  candidate_wells,
  tolerance = gate_threshold_tolerance,
  relaxed_tolerance = relaxed_gate_threshold_tolerance
) {
  candidate_wells <- candidate_wells %>%
    filter(assay == assay_name, run_id == run_id_value)

  if (nrow(candidate_wells) == 0L) {
    stop("No CJD/Ctrl candidate wells found for ", assay_name, " on plate ", run_id_value)
  }

  candidate_scores <- vector("list", 0L)
  for (i in seq_len(nrow(candidate_wells))) {
    candidate <- candidate_wells[i, ]
    parsed <- tryCatch(
      read_well_droplets_for_plot(candidate),
      error = function(err) {
        cat(
          "Skipping candidate final-adjusted-gate well ",
          candidate$well,
          " for ",
          assay_name,
          " on plate ",
          run_id_value,
          " while checking gate compatibility: ",
          conditionMessage(err),
          "\n",
          sep = ""
        )
        NULL
      }
    )
    if (is.null(parsed)) {
      next
    }

    ref_thresholds <- positive_thresholds %>%
      filter(!is.na(threshold), is.finite(threshold), channel %in% c("Ch1", "Ch2")) %>%
      arrange(channel, threshold)
    cand_thresholds <- parsed$thresholds %>%
      filter(!is.na(threshold), is.finite(threshold), channel %in% c("Ch1", "Ch2")) %>%
      arrange(channel, threshold)

    if (nrow(ref_thresholds) == 0L || nrow(cand_thresholds) == 0L) {
      next
    }
    if (nrow(ref_thresholds) != nrow(cand_thresholds) || !all(ref_thresholds$channel == cand_thresholds$channel)) {
      next
    }

    deltas <- abs(ref_thresholds$threshold - cand_thresholds$threshold)
    compatible_ok <- all(deltas <= tolerance)
    candidate_scores[[length(candidate_scores) + 1L]] <- candidate %>%
      mutate(
        final_gate_match_distance = max(deltas),
        final_gate_match_mean_distance = mean(deltas),
        final_gate_match_type = if_else(compatible_ok, "strict", "relaxed")
      )
    if (!compatible_ok) {
      next
    }
  }

  candidate_scores <- bind_rows(candidate_scores)
  if (nrow(candidate_scores) == 0L) {
    stop(
      "No same-plate real CJD/Ctrl candidate with comparable Ch1/Ch2 thresholds for ",
      assay_name,
      " on plate ",
      run_id_value
    )
  }

  intended_final_gate <- provided_final_gate_wells %>%
    filter(run_id == run_id_value, assay == assay_name)

  selected <- if (nrow(intended_final_gate) == 1L) {
    intended <- candidate_scores %>%
      filter(
        final_gate_match_type == "strict",
        well == intended_final_gate$well[[1]]
      )
    if (nrow(intended) != 1L) {
      stop(
        "Intended final-adjusted-gate well ",
        intended_final_gate$well[[1]],
        " was not a strict compatible real CJD/Ctrl candidate for ",
        assay_name,
        " on plate ",
        run_id_value
      )
    }
    intended
  } else if (any(candidate_scores$final_gate_match_type == "strict")) {
    candidate_scores %>%
      filter(final_gate_match_type == "strict") %>%
      arrange(
        desc(accepted_droplets),
        sample,
        well
      ) %>%
      slice_head(n = 1L)
  } else {
    relaxed <- candidate_scores %>%
      filter(final_gate_match_type == "relaxed", final_gate_match_distance <= relaxed_tolerance) %>%
      arrange(
        final_gate_match_distance,
        final_gate_match_mean_distance,
        desc(accepted_droplets),
        sample,
        well
      )
    if (nrow(relaxed) == 0L) {
      stop(
        "No same-plate real CJD/Ctrl candidate was within relaxed gate matching for ",
        assay_name,
        " on plate ",
        run_id_value
      )
    }
    relaxed %>%
      slice_head(n = 1L)
  }

  selected <- selected %>%
    mutate(stage = "Final adjusted gate")

  if (selected$final_gate_match_type %>% first() == "strict") {
    selected_count <- nrow(candidate_scores %>% filter(final_gate_match_type == "strict"))
    compatible_label <- paste0(
      selected_count,
      " strict threshold-matching candidate(s)"
    )
  } else {
    compatible_label <- "0 strict threshold-matching candidate(s)"
  }

  cat(
    "Final adjusted gate selected for ",
    assay_name,
    " on ",
    run_id_value,
    ": ",
    compatible_label,
    ", selected ",
    selected$well,
    " (sample ",
    selected$sample,
    "; match=",
    selected$final_gate_match_type,
    ", dmax=",
    format(signif(selected$final_gate_match_distance, 3)),
    ", dmean=",
    format(signif(selected$final_gate_match_mean_distance, 3)),
    ")\n",
    sep = ""
  )

  selected
}

# Resolve complete control wells plus a compatible final-adjusted gate per assay.
resolve_assay_control_set <- function(
  assay_name,
  control_runs_for_assay,
  available_control_wells,
  real_sample_control_wells,
  tolerance = gate_threshold_tolerance
) {
  if (nrow(control_runs_for_assay) == 0L) {
    stop("Could not identify complete control runs for assay ", assay_name)
  }

  for (i in seq_len(nrow(control_runs_for_assay))) {
    run_id_value <- control_runs_for_assay$run_id[[i]]
    controls <- available_control_wells %>%
      filter(assay == assay_name, run_id == run_id_value, stage %in% control_stage_order[1:3]) %>%
      group_by(stage) %>%
      arrange(desc(accepted_droplets), sample, well, .by_group = TRUE) %>%
      slice_head(n = 1L) %>%
      ungroup()

    if (nrow(controls) != length(control_stage_order[1:3])) {
      next
    }

    positive_thresholds <- tryCatch(
      {
        positive_row <- controls %>% filter(stage == "Positive control")
        if (nrow(positive_row) != 1L) {
          stop("Expected one positive-control row for assay ", assay_name, " on run ", run_id_value)
        }
        read_well_droplets_for_plot(positive_row)$thresholds
      },
      error = function(err) {
        cat(
          "Skipping control run ",
          run_id_value,
          " for assay ",
          assay_name,
          ": ",
          conditionMessage(err),
          "\n",
          sep = ""
        )
        NULL
      }
    )

    if (is.null(positive_thresholds)) {
      next
    }

    final_gate <- tryCatch(
      select_final_gate_control(
        assay_name = assay_name,
        run_id_value = run_id_value,
        positive_thresholds = positive_thresholds,
        candidate_wells = real_sample_control_wells,
        tolerance = tolerance
      ),
      error = function(err) {
        cat(
          "Control run ",
          run_id_value,
          " for assay ",
          assay_name,
          " has no same-plate compatible final-adjusted gate: ",
          conditionMessage(err),
          "\n",
          sep = ""
        )
        NULL
      }
    )

    if (is.null(final_gate)) {
      next
    }

    return(bind_rows(controls, final_gate))
  }

  stop(
    "Could not find a complete control run with a compatible same-plate final-adjusted gate for assay ",
    assay_name
  )
}

# Strategy plots use the same cluster-boundary threshold source for every stage.
thresholds_for_control_stage <- function(thresholds, stage) {
  thresholds
}

# Write paired SVG/PDF outputs for each individual asset.
save_plot_pair <- function(plot, output_dir, basename, width, height) {
  svg_path <- file.path(output_dir, paste0(basename, ".svg"))
  pdf_path <- file.path(output_dir, paste0(basename, ".pdf"))

  grDevices::svg(svg_path, width = width, height = height, onefile = FALSE)
  print(plot)
  grDevices::dev.off()

  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()

  tibble(svg_path = svg_path, pdf_path = pdf_path)
}

# ---- source manifests ----

# Active SNV runs define the raw well universe available for gating assets.
runs <- readr::read_csv(file.path(raw_root, "manifests", "runs.csv"), show_col_types = FALSE) %>%
  filter(status == "active", experiment == "SNV") %>%
  select(run_id, archive_contents_relative_dir)

# Collapse the raw sample manifest to one row per physical well.
well_manifest <- readr::read_csv(file.path(raw_root, "manifests", "sample_manifest.csv"), show_col_types = FALSE) %>%
  filter(experiment == "SNV", assay %in% mutation_order) %>%
  group_by(run_id, run_date, assay, well, sample) %>%
  summarise(accepted_droplets = max(accepted_droplets, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    run_date = as.Date(run_date),
    sample_key = normalise_sample_key(sample)
  ) %>%
  left_join(runs, by = "run_id")

# SNV_data_final.xlsx defines the LoB+LoD+ sample-region positives to render.
positive_rows <- openxlsx::read.xlsx(file.path(project_root, "results", "ddPCR", "SNV_data_final.xlsx")) %>%
  as_tibble() %>%
  mutate(
    detected_above_LoB = as_true_flag(detected_above_LoB),
    detected_above_LoD = as_true_flag(detected_above_LoD),
    is_pooled = as_true_flag(is_pooled),
    region_label = recode(brain_region, !!!region_labels),
    sample_id = paste(coalesce(as.character(code), as.character(participant)), brain_region, sep = "_"),
    sample_key = normalise_sample_key(sample_id),
    manuscript_label = paste(participant, region_label, mutation),
    positive_id = safe_file_component(paste(participant, region_label, mutation, sep = "_")),
    row_id = row_number()
  ) %>%
  filter(
    detected_above_LoB,
    detected_above_LoD,
    !is_pooled,
    !(participant == "CJD30" & mutation == "E200K")
  ) %>%
  arrange(match(mutation, mutation_order), participant, brain_region)

if (nrow(positive_rows) == 0L) {
  # Keep outputs deterministic when no positive rows are available.
  unlink(file.path(positive_individual_dir, list.files(positive_individual_dir)), force = TRUE)
  unlink(
    file.path(
      positive_out_dir,
      c(
        "ddpcr_lob_lod_positive_merged_panel.svg",
        "ddpcr_lob_lod_positive_merged_panel.pdf",
        "ddpcr_lob_lod_positive_faceted_panel.svg",
        "ddpcr_lob_lod_positive_faceted_panel.pdf"
      )
    ),
    force = TRUE
  )
  positive_manifest <- tibble(
    figure = character(),
    plot_kind = character(),
    plot_order = integer(),
    assay = character(),
    sample_label = character(),
    stage = character(),
    n_wells = integer(),
    n_droplets = integer(),
    basename = character(),
    svg_path = character(),
    pdf_path = character()
  )
  positive_wells <- tibble(
    row_id = integer(),
    manuscript_label = character(),
    mutation = character(),
    run_id = character(),
    run_date = as.Date(character()),
    well = character(),
    sample = character(),
    accepted_droplets = numeric(),
    archive_contents_relative_dir = character()
  )
  positive_thresholds <- tibble(
    channel = character(),
    threshold_type = character(),
    threshold = numeric(),
    run_id = character(),
    run_date = as.Date(character()),
    assay = character(),
    well = character(),
    sample = character(),
    row_id = integer(),
    positive_id = character(),
    manuscript_label = character(),
    well_plot_id = character()
  )
} else {
  # Map each positive sample-region result back to its contributing raw wells.
  positive_wells <- positive_rows %>%
    select(
      row_id, positive_id, manuscript_label, participant, brain_region,
      region_label, mutation, sample_id, sample_key,
      fractional_abundance, ci_low, ci_high
    ) %>%
    inner_join(
      well_manifest,
      by = c("mutation" = "assay", "sample_key" = "sample_key"),
      relationship = "many-to-many"
    ) %>%
    filter(!is.na(archive_contents_relative_dir)) %>%
    mutate(assay = mutation) %>%
    arrange(row_id, run_date, well)

  if (nrow(positive_wells) == 0L) {
    stop("Could not map LoB+LoD+ rows back to raw wells")
  }

  # Parse every contributing well once, then reuse the parsed droplets for merged
  # and faceted plots.
  positive_parsed <- positive_wells %>%
    split(seq_len(nrow(.))) %>%
    map(function(row) {
      row <- as_tibble(row)
      parsed <- read_well_droplets_for_plot(row)
      list(row = row, parsed = parsed)
    })

  # Attach manuscript identifiers to every parsed droplet.
  positive_droplets <- map2_dfr(positive_parsed, seq_along(positive_parsed), function(entry, i) {
    entry$parsed$droplets %>%
      mutate(
        row_id = entry$row$row_id,
        positive_id = entry$row$positive_id,
        manuscript_label = entry$row$manuscript_label,
        fractional_abundance = entry$row$fractional_abundance,
        ci_low = entry$row$ci_low,
        ci_high = entry$row$ci_high,
        well_plot_id = paste(entry$row$run_id, entry$row$well, sep = "_")
      )
  })

  # Keep thresholds in a separate table so plot manifests can audit gate lines.
  positive_thresholds <- map_dfr(positive_parsed, function(entry) {
    entry$parsed$thresholds %>%
      mutate(
        row_id = entry$row$row_id,
        positive_id = entry$row$positive_id,
        manuscript_label = entry$row$manuscript_label,
        well_plot_id = paste(entry$row$run_id, entry$row$well, sep = "_")
      )
  })

  # Use one axis frame for all positive-sample gating assets.
  positive_limits <- axis_limits(positive_droplets, positive_thresholds)

  # Render merged and faceted views for each LoB+LoD+ sample-region row.
  positive_manifest <- map_dfr(seq_len(nrow(positive_rows)), function(i) {
    row <- positive_rows[i, ]
    droplets <- positive_droplets %>% filter(row_id == row$row_id)
    thresholds <- positive_thresholds %>% filter(row_id == row$row_id)
    n_wells <- n_distinct(droplets$well_plot_id)
    subtitle <- paste0(
      row$mutation,
      "; ", n_wells, " well", ifelse(n_wells == 1L, "", "s"),
      "; FA ", format_pct(row$fractional_abundance),
      "% (95% CI ", format_pct(row$ci_low), "-", format_pct(row$ci_high), "%)"
    )

    merged_plot <- base_scatter_plot(
      droplets = droplets,
      thresholds = thresholds %>% filter(threshold_type == "manual"),
      title = row$manuscript_label,
      subtitle = paste0(subtitle, "; merged droplets"),
      mutation = row$mutation,
      limits = positive_limits,
      facet_by_well = FALSE,
      show_auto = FALSE,
      show_manual = TRUE
    )
    merged_base <- paste0("positive_merged_", row$positive_id)
    merged_files <- save_plot_pair(merged_plot, positive_individual_dir, merged_base, width = 5.2, height = 4.8)

    faceted_plot <- base_scatter_plot(
      droplets = droplets,
      thresholds = thresholds %>% filter(threshold_type == "manual"),
      title = row$manuscript_label,
      subtitle = paste0(subtitle, "; one facet per contributing well"),
      mutation = row$mutation,
      limits = positive_limits,
      facet_by_well = TRUE,
      show_auto = FALSE,
      show_manual = TRUE
    )
    faceted_base <- paste0("positive_faceted_", row$positive_id)
    faceted_files <- save_plot_pair(faceted_plot, positive_individual_dir, faceted_base, width = 5.2, height = max(4.8, 2.4 * n_wells))

    bind_rows(
      tibble(
        figure = "lob_lod_positive",
        plot_kind = "merged",
        plot_order = i,
        assay = row$mutation,
        sample_label = row$manuscript_label,
        stage = NA_character_,
        n_wells = n_wells,
        n_droplets = nrow(droplets),
        basename = merged_base
      ) %>% bind_cols(merged_files),
      tibble(
        figure = "lob_lod_positive",
        plot_kind = "faceted",
        plot_order = i,
        assay = row$mutation,
        sample_label = row$manuscript_label,
        stage = NA_character_,
        n_wells = n_wells,
        n_droplets = nrow(droplets),
        basename = faceted_base
      ) %>% bind_cols(faceted_files)
    )
  })
}

readr::write_csv(positive_manifest, file.path(positive_out_dir, "plot_manifest.csv"))

# Persist positive-panel source wells and thresholds for audit.
readr::write_csv(
  positive_wells %>%
    select(
      row_id, manuscript_label, mutation, run_id, run_date, well, sample,
      accepted_droplets, archive_contents_relative_dir
    ),
  file.path(positive_out_dir, "selected_wells.csv")
)
readr::write_csv(
  positive_thresholds,
  file.path(positive_out_dir, "selected_thresholds.csv")
)

# ---- gating strategy controls ----

# Add manifest-derived paths so we can select both hard-coded control wells and
# same-plate CJD/Ctrl candidate wells from one source table.
well_manifest_with_paths <- well_manifest %>%
  mutate(
    sample_upper = str_to_upper(sample),
    peak_path = file.path(raw_root, archive_contents_relative_dir, "PeakData", paste0(well, ".ddpeakjson")),
    metadata_path = file.path(raw_root, archive_contents_relative_dir, "PeakMetaData", paste0(well, ".ddmetajson")),
    has_peak_data = file.exists(peak_path),
    has_metadata = file.exists(metadata_path)
  )

# Identify NTC, WT, and positive-control wells with usable raw JSON data.
available_control_wells <- well_manifest_with_paths %>%
  mutate(
    stage = case_when(
      str_detect(sample_upper, "NTC") ~ "NTC",
      str_detect(sample_upper, "^WT|[_-]WT|WT[_-]") ~ "WT",
      str_detect(sample_upper, "MUT") ~ "Positive control",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(stage), has_peak_data, has_metadata)

# Candidate final-adjusted-gate wells must be real CJD/Ctrl samples, not controls.
real_sample_control_wells <- well_manifest_with_paths %>%
  filter(
    is_real_control_sample(sample_upper),
    !str_detect(sample_upper, "NTC|[_-]WT|WT[_-]|MUT"),
    has_peak_data,
    has_metadata
  )

# Prefer control runs close to the positive-sample runs when positives exist.
reference_dates <- if (nrow(positive_wells) > 0L) {
  unique(positive_wells$run_date)
} else {
  unique(well_manifest$run_date)
}

# Choose complete control runs for each assay.
control_runs <- available_control_wells %>%
  group_by(assay, run_id, run_date, archive_contents_relative_dir) %>%
  summarise(
    stages_present = paste(sort(unique(stage)), collapse = ";"),
    complete = all(control_stage_order[1:3] %in% stage),
    total_accepted = sum(accepted_droplets, na.rm = TRUE),
    date_distance = min(abs(as.integer(first(run_date) - reference_dates)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(complete) %>%
  arrange(match(assay, mutation_order), date_distance, desc(total_accepted), run_date)

# Resolve control rows and final-adjusted gate rows one assay at a time in
# ranked run order.
strategy_controls <- map_dfr(
  mutation_order,
  function(assay_name) {
    resolve_assay_control_set(
      assay_name = assay_name,
      control_runs_for_assay = control_runs %>% filter(assay == assay_name),
      available_control_wells = available_control_wells,
      real_sample_control_wells = real_sample_control_wells
    )
  }
) %>%
  mutate(
    stage = factor(stage, levels = control_stage_order),
    stage = as.character(stage)
  ) %>%
  arrange(match(assay, mutation_order), match(stage, control_stage_order))

if (nrow(strategy_controls) != length(mutation_order) * length(control_stage_order)) {
  stop("Could not resolve all required control and final-adjusted gate rows")
}

if (
  nrow(
    semi_join(
      filter(strategy_controls, stage == "Final adjusted gate"),
      filter(strategy_controls, stage == "Positive control"),
      by = c("assay", "run_id", "well", "sample")
    )
  ) > 0L
) {
  stop("Final adjusted gate selection duplicated a positive-control sample; expected a real CJD/Ctrl sample")
}

# Parse each selected control well once, including final adjusted gate rows.
strategy_parsed <- strategy_controls %>%
  split(seq_len(nrow(.))) %>%
  map(function(row) {
    row <- as_tibble(row)
    parsed <- read_well_droplets_for_plot(row)
    list(row = row, parsed = parsed)
  })

# Attach stage IDs to parsed control droplets and thresholds.
strategy_droplets <- map_dfr(strategy_parsed, function(entry) {
  entry$parsed$droplets %>%
    mutate(
      stage = entry$row$stage,
      control_plot_id = paste(entry$row$assay, entry$row$stage, sep = "_")
    )
})

strategy_thresholds <- map_dfr(strategy_parsed, function(entry) {
  entry$parsed$thresholds %>%
    mutate(
      stage = entry$row$stage,
      control_plot_id = paste(entry$row$assay, entry$row$stage, sep = "_")
    )
})

# Use one common axis frame per mutation for gating-strategy panels. The
# strategy plots are explanatory gate examples, so keep every accepted droplet
# visible instead of trimming to a high quantile.
strategy_limits_by_assay <- strategy_droplets %>%
  split(.$assay) %>%
  map(function(assay_droplets) {
    assay_name <- unique(assay_droplets$assay)
    assay_thresholds <- strategy_thresholds %>%
      filter(assay == assay_name)
    limits <- axis_limits(assay_droplets, assay_thresholds, droplet_quantile = 1)
    override <- strategy_axis_overrides %>%
      filter(assay == assay_name)
    if (nrow(override) == 1L) {
      if (!is.na(override$x_max)) {
        limits$x[2] <- override$x_max
      }
      if (!is.na(override$y_max)) {
        limits$y[2] <- override$y_max
      }
    }
    limits
  })

# Render one control-stage asset per selected strategy row.
strategy_manifest <- map_dfr(seq_len(nrow(strategy_controls)), function(i) {
  row <- strategy_controls[i, ]
  droplets <- strategy_droplets %>%
    filter(assay == row$assay, stage == row$stage, well == row$well, run_id == row$run_id)
  thresholds <- strategy_thresholds %>%
    filter(assay == row$assay, stage == row$stage, well == row$well, run_id == row$run_id)
  thresholds_to_plot <- thresholds_for_control_stage(thresholds, row$stage)
  title <- paste(row$assay, row$stage)

  plot <- base_scatter_plot(
    droplets = droplets,
    thresholds = thresholds_to_plot,
    title = title,
    subtitle = NULL,
    mutation = row$assay,
    limits = strategy_limits_by_assay[[row$assay]],
    facet_by_well = FALSE,
    show_auto = TRUE,
    show_gates = row$stage != "NTC",
    show_manual = TRUE,
    show_legend = FALSE,
    plot_title_size = 11
  )

  basename <- paste0(
    "strategy_",
    safe_file_component(row$assay),
    "_",
    safe_file_component(row$stage)
  )
  files <- save_plot_pair(plot, strategy_individual_dir, basename, width = 4.2, height = 3.8)
  tibble(
    figure = "gating_strategy",
    plot_kind = "control",
    plot_order = i,
    assay = row$assay,
    sample_label = row$sample,
    stage = row$stage,
    n_wells = 1L,
    n_droplets = nrow(droplets),
    basename = basename,
    run_id = row$run_id,
    run_date = as.Date(row$run_date),
    well = row$well
  ) %>%
    bind_cols(files)
})

readr::write_csv(strategy_manifest, file.path(strategy_out_dir, "plot_manifest.csv"))

# Persist strategy-panel source wells and thresholds for audit.
readr::write_csv(
  strategy_controls %>%
    select(assay, stage, run_id, run_date, well, sample, accepted_droplets, archive_contents_relative_dir),
  file.path(strategy_out_dir, "selected_wells.csv")
)
readr::write_csv(
  strategy_thresholds,
  file.path(strategy_out_dir, "selected_thresholds.csv")
)

# ---- run summary ----

cat("LoB+LoD+ individual plots written: ", nrow(positive_manifest), "\n", sep = "")
cat("Gating-strategy individual plots written: ", nrow(strategy_manifest), "\n", sep = "")
