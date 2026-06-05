library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)
library(jsonlite)
library(ggplot2)
library(openxlsx)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
raw_root <- file.path(project_root, "raw", "ddpcr")
positive_out_dir <- file.path(project_root, "manuscript", "figures", "ddpcr_gating_lob_lod_positive")
strategy_out_dir <- file.path(project_root, "manuscript", "figures", "ddpcr_gating_strategy")
positive_individual_dir <- file.path(positive_out_dir, "individual")
strategy_individual_dir <- file.path(strategy_out_dir, "individual")
dir.create(positive_individual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(strategy_individual_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "src", "ddpcr", "ddpcr_raw_import_helpers.R"))

mutation_order <- c("D178N", "E200K", "P102L")
region_labels <- c(
  bg = "basal ganglia",
  cb = "cerebellum",
  fr = "frontal cortex",
  hc = "hippocampus",
  ps = "pons",
  sn = "substantia nigra",
  th = "thalamus"
)
control_stage_order <- c("NTC", "WT", "Positive control", "Final adjusted gate")

class_levels <- c(
  "Reference-only",
  "Mutant-only",
  "Double-positive",
  "Double-negative",
  "Gated/unassigned",
  "Rejected/unassigned"
)
class_colours <- c(
  `Reference-only` = "#0072B2",
  `Mutant-only` = "#D55E00",
  `Double-positive` = "#CC79A7",
  `Double-negative` = "#9CA3AF",
  `Gated/unassigned` = "#E69F00",
  `Rejected/unassigned` = "#E5E7EB"
)
draw_order <- c(
  "Rejected/unassigned" = 1L,
  "Double-negative" = 2L,
  "Gated/unassigned" = 3L,
  "Reference-only" = 4L,
  "Mutant-only" = 5L,
  "Double-positive" = 6L
)

safe_file_component <- function(x) {
  x <- str_squish(as.character(x))
  x <- str_replace_all(x, "[^A-Za-z0-9._-]+", "_")
  x <- str_replace_all(x, "_+", "_")
  str_replace_all(x, "^_+|_+$", "")
}

as_true_flag <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

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

format_pct <- function(x) {
  format(round(as.numeric(x), 3), nsmall = 3, trim = TRUE, scientific = FALSE)
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

threshold_rows_for_channel <- function(metadata, channel) {
  manual <- threshold_entry_value((metadata$ThresholdValues %||% list())[[channel]])
  auto <- threshold_entry_value((metadata$AutoThresholdValues %||% list())[[channel]])
  tibble(
    channel = paste0("Ch", channel),
    threshold_type = c("manual", "auto"),
    threshold = c(manual, auto)
  ) %>%
    filter(!is.na(threshold), is.finite(threshold), threshold > 0)
}

thresholds_from_metadata <- function(metadata) {
  bind_rows(
    threshold_rows_for_channel(metadata, 1L),
    threshold_rows_for_channel(metadata, 2L)
  )
}

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
  amplitudes <- peak$PeakInfo$Amplitudes
  if (is.null(amplitudes) || length(amplitudes) < 2L) {
    stop("Peak data does not contain two amplitude channels: ", peak_path)
  }

  targets <- metadata_targets(metadata)
  selected <- selected_target_indices(targets, well_row$assay)
  if (is.null(targets) || is.null(selected)) {
    stop("Could not select WT and mutant targets for ", well_row$run_id, " ", well_row$well)
  }

  target_names <- vapply(targets, target_name, character(1))
  channels <- vapply(targets, target_channel, integer(1))
  selected_indices <- c(selected$ref, selected$mut)
  ch1_idx <- selected_indices[channels[selected_indices] == 1L][1]
  ch2_idx <- selected_indices[channels[selected_indices] == 2L][1]
  if (is.na(ch1_idx) || is.na(ch2_idx)) {
    stop("Could not map selected targets to Ch1 and Ch2 for ", well_row$run_id, " ", well_row$well)
  }

  ch1 <- as.numeric(amplitudes[[1]])
  ch2 <- as.numeric(amplitudes[[2]])
  droplet_count <- min(length(ch1), length(ch2))
  droplets <- tibble(
    droplet_index = seq.int(0L, droplet_count - 1L),
    ch1_amplitude = ch1[seq_len(droplet_count)],
    ch2_amplitude = ch2[seq_len(droplet_count)],
    droplet_class = "Rejected/unassigned"
  )

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
      next
    }

    ref_result <- results[[selected$ref]]
    mut_result <- results[[selected$mut]]
    droplets$droplet_class[row_indices] <- class_from_calls(ref_result, mut_result)
  }

  thresholds <- thresholds_from_metadata(metadata) %>%
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

axis_limits <- function(droplets, thresholds) {
  x_candidate <- max(
    quantile(droplets$ch1_amplitude, 0.997, na.rm = TRUE, names = FALSE),
    thresholds$threshold[thresholds$channel == "Ch1"],
    na.rm = TRUE
  )
  y_candidate <- max(
    quantile(droplets$ch2_amplitude, 0.997, na.rm = TRUE, names = FALSE),
    thresholds$threshold[thresholds$channel == "Ch2"],
    na.rm = TRUE
  )
  list(
    x = c(0, x_candidate * 1.08),
    y = c(0, y_candidate * 1.08)
  )
}

add_gate_layers <- function(plot, thresholds, limits, show_auto = TRUE, show_manual = TRUE) {
  threshold_data <- thresholds %>%
    filter(
      (threshold_type == "auto" & show_auto) |
        (threshold_type == "manual" & show_manual)
    )
  if (nrow(threshold_data) == 0L) {
    return(plot +
      annotate(
        "text",
        x = limits$x[[1]] + diff(limits$x) * 0.03,
        y = limits$y[[2]] - diff(limits$y) * 0.04,
        label = "QuantaSoft cluster gate; no threshold line exported",
        hjust = 0,
        vjust = 1,
        size = 2.2,
        colour = "#374151"
      ))
  }

  plot <- plot +
    scale_linetype_manual(
      values = c(manual = "solid", auto = "dotted"),
      guide = "none"
    )

  vline_data <- threshold_data %>% filter(channel == "Ch1")
  hline_data <- threshold_data %>% filter(channel == "Ch2")
  if (nrow(vline_data) > 0L) {
    plot <- plot +
      geom_vline(
        data = vline_data,
        aes(xintercept = threshold, linetype = threshold_type),
        linewidth = 0.35,
        colour = "black",
        alpha = 0.75,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  if (nrow(hline_data) > 0L) {
    plot <- plot +
      geom_hline(
        data = hline_data,
        aes(yintercept = threshold, linetype = threshold_type),
        linewidth = 0.35,
        colour = "black",
        alpha = 0.75,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  plot
}

base_scatter_plot <- function(
  droplets,
  thresholds,
  title,
  subtitle,
  limits,
  facet_by_well = FALSE,
  show_auto = TRUE,
  show_manual = TRUE,
  show_legend = TRUE
) {
  plot <- droplets %>%
    arrange(draw_order, droplet_index) %>%
    ggplot(aes(x = ch1_amplitude, y = ch2_amplitude, colour = droplet_class)) +
    geom_point(size = 0.18, alpha = 0.42, stroke = 0, na.rm = TRUE) +
    scale_colour_manual(values = class_colours, drop = FALSE, name = "Droplet class") +
    coord_cartesian(xlim = limits$x, ylim = limits$y, expand = FALSE) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Channel 1 amplitude",
      y = "Channel 2 amplitude"
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
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7),
      strip.background = element_rect(fill = "#F3F4F6", colour = "#D1D5DB"),
      strip.text = element_text(size = 7)
    )

  if (facet_by_well) {
    plot <- plot + facet_wrap(~well_label, ncol = 1)
  }

  add_gate_layers(
    plot = plot,
    thresholds = thresholds,
    limits = limits,
    show_auto = show_auto,
    show_manual = show_manual
  )
}

thresholds_for_control_stage <- function(thresholds, stage) {
  if (stage != "Final adjusted gate") {
    return(thresholds)
  }
  manual_thresholds <- thresholds %>% filter(threshold_type == "manual")
  if (nrow(manual_thresholds) > 0L) {
    return(manual_thresholds)
  }
  thresholds %>% filter(threshold_type == "auto")
}

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

runs <- readr::read_csv(file.path(raw_root, "manifests", "runs.csv"), show_col_types = FALSE) %>%
  filter(status == "active", experiment == "SNV") %>%
  select(run_id, archive_contents_relative_dir)

well_manifest <- readr::read_csv(file.path(raw_root, "manifests", "sample_manifest.csv"), show_col_types = FALSE) %>%
  filter(experiment == "SNV", assay %in% mutation_order) %>%
  group_by(run_id, run_date, assay, well, sample) %>%
  summarise(accepted_droplets = max(accepted_droplets, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    run_date = as.Date(run_date),
    sample_key = normalise_sample_key(sample)
  ) %>%
  left_join(runs, by = "run_id")

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
  stop("No non-germline sample-region LoB+LoD+ rows found")
}

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

positive_parsed <- positive_wells %>%
  split(seq_len(nrow(.))) %>%
  map(function(row) {
    row <- as_tibble(row)
    parsed <- read_well_droplets_for_plot(row)
    list(row = row, parsed = parsed)
  })

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

positive_thresholds <- map_dfr(positive_parsed, function(entry) {
  entry$parsed$thresholds %>%
    mutate(
      row_id = entry$row$row_id,
      positive_id = entry$row$positive_id,
      manuscript_label = entry$row$manuscript_label,
      well_plot_id = paste(entry$row$run_id, entry$row$well, sep = "_")
    )
})

positive_limits <- axis_limits(positive_droplets, positive_thresholds)

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

readr::write_csv(positive_manifest, file.path(positive_out_dir, "plot_manifest.csv"))
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

available_control_wells <- well_manifest %>%
  mutate(
    sample_upper = str_to_upper(sample),
    stage = case_when(
      str_detect(sample_upper, "NTC") ~ "NTC",
      str_detect(sample_upper, "^WT|[_-]WT|WT[_-]") ~ "WT",
      str_detect(sample_upper, "MUT") ~ "Positive control",
      TRUE ~ NA_character_
    ),
    peak_path = file.path(raw_root, archive_contents_relative_dir, "PeakData", paste0(well, ".ddpeakjson")),
    metadata_path = file.path(raw_root, archive_contents_relative_dir, "PeakMetaData", paste0(well, ".ddmetajson")),
    has_peak_data = file.exists(peak_path),
    has_metadata = file.exists(metadata_path)
  ) %>%
  filter(!is.na(stage), has_peak_data, has_metadata)

reference_dates <- unique(positive_wells$run_date)

control_runs <- available_control_wells %>%
  group_by(assay, run_id, run_date, archive_contents_relative_dir) %>%
  summarise(
    stages_present = paste(sort(unique(stage)), collapse = ";"),
    complete = all(control_stage_order[1:3] %in% stage),
    total_accepted = sum(accepted_droplets, na.rm = TRUE),
    date_distance = min(abs(as.integer(run_date - reference_dates))),
    .groups = "drop"
  ) %>%
  filter(complete) %>%
  arrange(match(assay, mutation_order), date_distance, desc(total_accepted), run_date)

selected_control_runs <- control_runs %>%
  group_by(assay) %>%
  slice_head(n = 1) %>%
  ungroup()

if (nrow(selected_control_runs) != length(mutation_order)) {
  stop("Could not identify complete NTC/WT/positive-control sets for every assay")
}

selected_controls <- selected_control_runs %>%
  select(assay, run_id) %>%
  left_join(
    available_control_wells,
    by = c("assay", "run_id")
  ) %>%
  filter(stage %in% control_stage_order[1:3]) %>%
  group_by(assay, stage) %>%
  arrange(desc(accepted_droplets), well, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

final_gate_controls <- selected_controls %>%
  filter(stage == "Positive control") %>%
  mutate(stage = "Final adjusted gate")

strategy_controls <- bind_rows(selected_controls, final_gate_controls) %>%
  mutate(
    stage = factor(stage, levels = control_stage_order),
    stage = as.character(stage)
  ) %>%
  arrange(match(assay, mutation_order), match(stage, control_stage_order))

strategy_parsed <- strategy_controls %>%
  split(seq_len(nrow(.))) %>%
  map(function(row) {
    row <- as_tibble(row)
    parsed <- read_well_droplets_for_plot(row)
    list(row = row, parsed = parsed)
  })

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

strategy_limits <- axis_limits(strategy_droplets, strategy_thresholds)

strategy_manifest <- map_dfr(seq_len(nrow(strategy_controls)), function(i) {
  row <- strategy_controls[i, ]
  droplets <- strategy_droplets %>%
    filter(assay == row$assay, stage == row$stage, well == row$well, run_id == row$run_id)
  thresholds <- strategy_thresholds %>%
    filter(assay == row$assay, stage == row$stage, well == row$well, run_id == row$run_id)
  thresholds_to_plot <- thresholds_for_control_stage(thresholds, row$stage)
  title <- paste(row$assay, row$stage)
  subtitle <- paste0(
    as.Date(row$run_date), " ", row$well, " ", row$sample,
    "; n=", nrow(droplets)
  )

  plot <- base_scatter_plot(
    droplets = droplets,
    thresholds = thresholds_to_plot,
    title = title,
    subtitle = subtitle,
    limits = strategy_limits,
    facet_by_well = FALSE,
    show_auto = TRUE,
    show_manual = TRUE,
    show_legend = FALSE
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
readr::write_csv(
  strategy_controls %>%
    select(assay, stage, run_id, run_date, well, sample, accepted_droplets, archive_contents_relative_dir),
  file.path(strategy_out_dir, "selected_wells.csv")
)
readr::write_csv(
  strategy_thresholds,
  file.path(strategy_out_dir, "selected_thresholds.csv")
)

cat("LoB+LoD+ individual plots written: ", nrow(positive_manifest), "\n", sep = "")
cat("Gating-strategy individual plots written: ", nrow(strategy_manifest), "\n", sep = "")
