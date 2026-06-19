source("scratch/ddpcr_gating_method_comparison_v2/code/lib_v2.R")

setup_v2_dirs()

summary_plot_dir <- path_v2("plots", "individual", "method_summary")
gating_plot_dir <- path_v2("plots", "individual", "method_gating")
panel_dir <- path_v2("plots", "panels")
report_dir <- path_v2("report")
dir.create(summary_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gating_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

json_text <- function(x) {
  as.character(jsonlite::toJSON(x, auto_unbox = TRUE))
}

current_json_counts <- function() {
  shared <- load_shared_droplets()
  shared %>%
    mutate(target_class = factor(target_class, levels = target_class_levels)) %>%
    group_by(run_id, assay, well, sample, sample_type, sample_key, run_date) %>%
    summarise(
      method_id = "current_json_quanta_gating",
      method_variant = "current exported JSON target-class calls",
      n_total_droplets = n(),
      n_accepted_droplets = sum(target_class %in% c("NN", "WT", "MUT", "DP")),
      n_nn = sum(target_class == "NN"),
      n_wt_only = sum(target_class == "WT"),
      n_mut_only = sum(target_class == "MUT"),
      n_double_positive = sum(target_class == "DP"),
      n_rain = sum(target_class == "Rain"),
      n_unclassified = sum(target_class == "Unclassified"),
      n_mut_expected = n_mut_only + n_double_positive,
      n_wt_expected = n_wt_only + n_double_positive,
      classification_status = "ok",
      classification_message = NA_character_,
      method_parameters_json = json_text(list(source = "exported_quantasoft_json")),
      .groups = "drop"
    )
}

read_method_counts <- function() {
  count_files <- c(
    path_v2("tables", "twoddpcr_well_counts.csv"),
    path_v2("tables", "ddPCRclust_well_counts.csv"),
    path_v2("tables", "dpcp_well_counts.csv"),
    path_v2("tables", "one_channel_well_counts.csv"),
    path_v2("tables", "bayesian_well_counts.csv"),
    path_v2("tables", "polygon_well_counts.csv")
  )

  bind_rows(
    current_json_counts(),
    bind_rows(lapply(count_files, read_csv, show_col_types = FALSE))
  ) %>%
    mutate(
      run_date = as.Date(run_date),
      method_family = case_when(
        str_starts(method_id, "twoddpcr") ~ "twoddpcr",
        str_starts(method_id, "ddPCRclust") ~ "ddPCRclust",
        str_starts(method_id, "dPCP") ~ "dPCP",
        str_starts(method_id, "ddpcRquant") ~ "ddpcRquant",
        str_starts(method_id, "definetherain") ~ "definetherain",
        str_starts(method_id, "bayesian") ~ "bayesian",
        str_starts(method_id, "polygon") ~ "polygon_gates",
        str_starts(method_id, "current") ~ "current",
        TRUE ~ "other"
      )
    )
}

calculate_fa_table <- function(counts) {
  ok_counts <- counts %>%
    filter(classification_status == "ok") %>%
    mutate(
      ref_positive = n_wt_expected,
      mut_positive = n_mut_expected,
      total = n_accepted_droplets,
      ref_negative = total - ref_positive,
      mut_negative = total - mut_positive
    )

  fa_rows <- pmap_dfr(
    list(
      ref_positive = ok_counts$ref_positive,
      ref_negative = ok_counts$ref_negative,
      mut_positive = ok_counts$mut_positive,
      mut_negative = ok_counts$mut_negative,
      total = ok_counts$total
    ),
    function(ref_positive, ref_negative, mut_positive, mut_negative, total) {
      fractional_abundance(
        ref_positive = ref_positive,
        ref_negative = ref_negative,
        mut_positive = mut_positive,
        mut_negative = mut_negative,
        total = total
      )
    }
  )

  ok_counts %>%
    bind_cols(fa_rows) %>%
    rename(
      fractional_abundance = fractional_abundance,
      ci_low = ci_low,
      ci_high = ci_high
    )
}

sample_region_results <- function(counts) {
  pooled <- counts %>%
    group_by(method_id, method_variant, method_family, assay, sample, sample_type, sample_key) %>%
    summarise(
      classification_status = ifelse(any(classification_status == "failed"), "failed", "ok"),
      classification_message = paste(unique(na.omit(classification_message)), collapse = "; "),
      n_wells = n(),
      run_dates = paste(sort(unique(as.character(run_date))), collapse = ";"),
      n_total_droplets = sum(n_total_droplets, na.rm = TRUE),
      n_accepted_droplets = sum(n_accepted_droplets, na.rm = TRUE),
      n_mut_expected = sum(n_mut_expected, na.rm = TRUE),
      n_wt_expected = sum(n_wt_expected, na.rm = TRUE),
      n_nn = sum(n_nn, na.rm = TRUE),
      n_wt_only = sum(n_wt_only, na.rm = TRUE),
      n_mut_only = sum(n_mut_only, na.rm = TRUE),
      n_double_positive = sum(n_double_positive, na.rm = TRUE),
      n_rain = sum(n_rain, na.rm = TRUE),
      n_unclassified = sum(n_unclassified, na.rm = TRUE),
      .groups = "drop"
    )

  fa <- calculate_fa_table(pooled)

  blanks <- counts %>%
    filter(classification_status == "ok", sample_type == "WT_control", n_accepted_droplets >= 10000) %>%
    group_by(method_id, assay, plate = run_date) %>%
    summarise(
      x_blank = sum(n_mut_expected, na.rm = TRUE),
      n_blank = sum(n_accepted_droplets, na.rm = TRUE),
      n_wells_blank = n(),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      p0_upper = binom::binom.confint(
        x = ceiling(x_blank),
        n = round(n_blank),
        methods = "exact"
      )$upper
    ) %>%
    ungroup()

  fallback <- blanks %>%
    group_by(method_id, assay) %>%
    summarise(x_blank = sum(x_blank), n_blank = sum(n_blank), .groups = "drop") %>%
    rowwise() %>%
    mutate(
      p0_upper_fallback = binom::binom.confint(
        x = ceiling(x_blank),
        n = round(n_blank),
        methods = "exact"
      )$upper
    ) %>%
    ungroup() %>%
    select(method_id, assay, p0_upper_fallback)

  sample_plates <- counts %>%
    select(method_id, assay, sample, plate = run_date) %>%
    distinct()

  sample_lob <- sample_plates %>%
    left_join(blanks, by = c("method_id", "assay", "plate")) %>%
    left_join(fallback, by = c("method_id", "assay")) %>%
    mutate(p0_per_plate = coalesce(p0_upper, p0_upper_fallback, 0)) %>%
    group_by(method_id, assay, sample) %>%
    slice_max(p0_per_plate, n = 1, with_ties = FALSE) %>%
    summarise(
      p0_use = first(p0_per_plate),
      plate_p0_max = first(plate),
      .groups = "drop"
    )

  fa %>%
    left_join(sample_lob, by = c("method_id", "assay", "sample")) %>%
    mutate(
      p0_use = coalesce(p0_use, 0),
      lob_count = qbinom(0.95, size = round(n_accepted_droplets), prob = p0_use),
      lob_fa = ifelse(n_accepted_droplets > 0, fa_scale * lob_count / n_accepted_droplets, NA_real_),
      detected_above_LoB = n_mut_expected > lob_count,
      detected_above_LoD = ci_low > lod_cut[assay],
      is_germline_e200k = assay == "E200K" & str_detect(sample, "CJD30")
    )
}

plot_pass_heatmap <- function(pass_matrix) {
  ggplot(pass_matrix, aes(x = method_id, y = assay, fill = n_detected_lod)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = paste0(n_detected_lod, "/", n_biological_samples)), size = 2.1) +
    scale_fill_viridis_c(option = "C", name = "LoD+ samples") +
    labs(x = "Method", y = "Assay", title = "LoD-positive biological sample-regions by method") +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.grid = element_blank()
    )
}

plot_failure_rates <- function(method_summary) {
  ggplot(method_summary, aes(x = method_id, y = failed_well_fraction, fill = method_family)) +
    geom_col(width = 0.75) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Method", y = "Failed well rows", title = "Classifier failure rate by method") +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

safe_package_version <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    return(NA_character_)
  }
  as.character(utils::packageVersion(package))
}

package_versions_table <- function() {
  tibble(
    component = c(
      "R",
      "twoddpcr",
      "ddPCRclust",
      "dPCP",
      "dpcR",
      "ggplot2",
      "readr",
      "binom"
    ),
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      safe_package_version("twoddpcr"),
      safe_package_version("ddPCRclust"),
      safe_package_version("dPCP"),
      safe_package_version("dpcR"),
      safe_package_version("ggplot2"),
      safe_package_version("readr"),
      safe_package_version("binom")
    )
  )
}

representative_gating_methods <- c(
  "current_json_quanta_gating",
  "polygon_control_hull_gates",
  "twoddpcr_kmeans_control_mah_rain",
  "twoddpcr_knn_k11_mah_rain",
  "ddPCRclust_control_projected",
  "dPCP_control_projected",
  "definetherain_channel_combined",
  "ddpcRquant_0999",
  "bayesian_control_mixture_weighted"
)

current_plot_droplets <- function() {
  load_shared_droplets() %>%
    filter(sample_type %in% c("NTC", "WT_control", "positive_control")) %>%
    group_by(assay, sample_type, target_class) %>%
    group_modify(function(.x, .y) {
      .x %>% slice_sample(n = min(nrow(.x), 250L))
    }) %>%
    ungroup() %>%
    transmute(
      method_id = "current_json_quanta_gating",
      method_variant = "current exported JSON target-class calls",
      assay,
      run_id,
      well,
      sample,
      sample_type,
      ch1_amplitude,
      ch2_amplitude,
      target_class,
      assigned_target_class = factor(target_class, levels = target_class_levels),
      method_parameters_json = json_text(list(source = "exported_quantasoft_json"))
    )
}

read_plot_droplets <- function() {
  droplet_files <- c(
    path_v2("data", "droplets", "twoddpcr_plot_droplets.rds"),
    path_v2("data", "droplets", "ddPCRclust_plot_droplets.rds"),
    path_v2("data", "droplets", "dPCP_plot_droplets.rds"),
    path_v2("data", "droplets", "one_channel_plot_droplets.rds"),
    path_v2("data", "droplets", "bayesian_mixture_plot_droplets.rds"),
    path_v2("data", "droplets", "polygon_control_hull_gates_plot_droplets.rds")
  )

  bind_rows(
    current_plot_droplets(),
    bind_rows(lapply(droplet_files[file.exists(droplet_files)], readRDS))
  ) %>%
    filter(
      method_id %in% representative_gating_methods,
      sample_type %in% c("NTC", "WT_control", "positive_control")
    ) %>%
    mutate(
      method_id = factor(method_id, levels = representative_gating_methods),
      assigned_target_class = factor(as.character(assigned_target_class), levels = target_class_levels)
    )
}

write_representative_gating_plots <- function(plot_droplets) {
  plot_rows <- list()
  plot_objects <- list()
  sampled <- plot_droplets %>%
    group_by(method_id, assay, sample_type, assigned_target_class) %>%
    group_modify(function(.x, .y) {
      .x %>% slice_sample(n = min(nrow(.x), 40L))
    }) %>%
    ungroup()

  for (assay_name in mutation_order) {
    assay_data <- sampled %>% filter(assay == assay_name)
    axis_data <- plot_droplets %>% filter(assay == assay_name)
    x_max <- quantile(axis_data$ch1_amplitude, 0.999, na.rm = TRUE) * 1.05
    y_max <- quantile(axis_data$ch2_amplitude, 0.999, na.rm = TRUE) * 1.05

    plot <- ggplot(assay_data, aes(ch1_amplitude, ch2_amplitude, colour = assigned_target_class)) +
      geom_point(size = 0.18, alpha = 0.45, stroke = 0) +
      facet_grid(sample_type ~ method_id) +
      coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max), expand = FALSE) +
      scale_colour_viridis_d(option = "D", end = 0.9, drop = FALSE) +
      labs(
        title = paste(assay_name, "representative control droplet classifications"),
        x = "Channel 1 amplitude",
        y = "Channel 2 amplitude",
        colour = "Assigned class"
      ) +
      theme_bw(base_size = 6) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.text.x = element_text(size = 5),
        strip.text.y = element_text(size = 6)
      )

    svg_path <- file.path(gating_plot_dir, paste0(assay_name, "_representative_gating.svg"))
    pdf_path <- file.path(gating_plot_dir, paste0(assay_name, "_representative_gating.pdf"))
    ggsave(svg_path, plot, width = 15, height = 7)
    ggsave(pdf_path, plot, width = 15, height = 7, useDingbats = FALSE)

    plot_rows[[assay_name]] <- tibble(
      plot = paste0(assay_name, "_representative_gating"),
      assay = assay_name,
      svg_path = svg_path,
      pdf_path = pdf_path
    )
    plot_objects[[assay_name]] <- plot
  }

  list(
    manifest = bind_rows(plot_rows),
    plots = plot_objects
  )
}

write_plot_pair <- function(plot, stem, width, height) {
  svg_path <- file.path(summary_plot_dir, paste0(stem, ".svg"))
  pdf_path <- file.path(summary_plot_dir, paste0(stem, ".pdf"))
  ggsave(svg_path, plot, width = width, height = height)
  ggsave(pdf_path, plot, width = width, height = height, useDingbats = FALSE)
  tibble(plot = stem, svg_path = svg_path, pdf_path = pdf_path)
}

write_panel_pdf <- function(heatmap_plot, failure_plot) {
  pdf_path <- file.path(panel_dir, "method_summary_panel.pdf")
  grDevices::pdf(pdf_path, width = 12, height = 9, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  print(heatmap_plot, vp = grid::viewport(x = 0.5, y = 0.72, width = 0.98, height = 0.52))
  print(failure_plot, vp = grid::viewport(x = 0.5, y = 0.24, width = 0.98, height = 0.40))
  pdf_path
}

write_report_pdf <- function(method_summary, pass_matrix, positive_rows, e2e_summary,
                             package_versions, heatmap_plot, failure_plot,
                             gating_plots) {
  pdf_path <- file.path(report_dir, "ddpcr_gating_method_comparison_v2.pdf")
  grDevices::pdf(pdf_path, width = 11, height = 8.5, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)

  text_page <- function(title, lines) {
    grid::grid.newpage()
    grid::grid.text(title, x = 0.05, y = 0.94, just = "left", gp = grid::gpar(fontsize = 16, fontface = "bold"))
    grid::grid.text(paste(lines, collapse = "\n"), x = 0.05, y = 0.86, just = c("left", "top"), gp = grid::gpar(fontsize = 9))
  }

  table_lines <- function(df, n = 60) {
    shown <- min(nrow(df), n)
    body <- capture.output(print(as.data.frame(head(df, shown)), row.names = FALSE))
    if (shown < nrow(df)) {
      body <- c(body, paste("Showing first", shown, "of", nrow(df), "rows."))
    }
    body
  }

  text_page(
    "ddPCR gating method comparison v2",
    c(
      paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      paste("Methods:", nrow(method_summary)),
      paste("Biological LoD+ rows:", nrow(positive_rows)),
      "Counts use all droplets; plots use deterministic method-level samples.",
      "LoB uses WT genomic controls with same-plate preference and assay-wide fallback.",
      "LoD thresholds are D178N 0.056%, E200K 0.067%, and P102L 0.13%."
    )
  )

  text_page(
    "Evidence basis and interpretation",
    c(
      "dMIQE 2020 motivates transparent reporting of thresholds, controls, analysis software, accepted partitions, and rain.",
      "Bio-Rad QuantaSoft Analysis Pro documents 2D line and cluster modes, including square, circle, and free-form cluster shapes.",
      "The ddpcr package overview frames two-channel assays as four amplitude clusters and notes that manual gates are common when automatic gates are inadequate.",
      "The polygon comparator is therefore included as a control-anchored free-form manual-gate analogue.",
      "One-channel methods are sensitivity comparators because they combine independent channel thresholds rather than modelling the full two-dimensional droplet cloud.",
      "Native package failures are retained as failed rows; control-projected adapters are labelled separately and should not be read as native package success."
    )
  )

  text_page(
    "Package versions",
    table_lines(package_versions, n = 40)
  )

  text_page(
    "Method status summary",
    table_lines(method_summary %>%
      select(method_id, method_family, ok_well_rows, failed_well_rows, n_detected_lod) %>%
      arrange(method_family, method_id), n = 40)
  )

  grid::grid.newpage()
  print(heatmap_plot)
  grid::grid.newpage()
  print(failure_plot)
  for (assay_name in names(gating_plots)) {
    grid::grid.newpage()
    print(gating_plots[[assay_name]])
  }

  text_page(
    "E2E check summary",
    table_lines(e2e_summary, n = 60)
  )

  text_page(
    "LoD-positive biological rows",
    table_lines(positive_rows %>%
      select(method_id, assay, sample, fractional_abundance, ci_low, ci_high, lob_count, detected_above_LoD) %>%
      arrange(method_id, assay, sample), n = 80)
  )

  pdf_path
}

counts <- read_method_counts()
sample_results <- sample_region_results(counts)

pass_matrix <- sample_results %>%
  filter(sample_type == "biological_sample", !is_germline_e200k) %>%
  group_by(method_id, method_family, assay) %>%
  summarise(
    n_biological_samples = n(),
    n_detected_lob = sum(detected_above_LoB, na.rm = TRUE),
    n_detected_lod = sum(detected_above_LoD, na.rm = TRUE),
    n_failed = sum(classification_status == "failed"),
    .groups = "drop"
  )

positive_rows <- sample_results %>%
  filter(sample_type == "biological_sample", !is_germline_e200k) %>%
  filter(detected_above_LoB | detected_above_LoD) %>%
  arrange(method_id, assay, sample)

method_summary <- counts %>%
  group_by(method_id, method_variant, method_family) %>%
  summarise(
    well_rows = n(),
    ok_well_rows = sum(classification_status == "ok"),
    failed_well_rows = sum(classification_status == "failed"),
    failed_well_fraction = failed_well_rows / well_rows,
    .groups = "drop"
  ) %>%
  left_join(
    pass_matrix %>%
      group_by(method_id) %>%
      summarise(
        n_detected_lob = sum(n_detected_lob),
        n_detected_lod = sum(n_detected_lod),
        .groups = "drop"
      ),
    by = "method_id"
  ) %>%
  mutate(
    n_detected_lob = replace_na(n_detected_lob, 0),
    n_detected_lod = replace_na(n_detected_lod, 0)
  )

control_false_positive <- sample_results %>%
  filter(sample_type == "WT_control") %>%
  group_by(method_id, method_family, assay) %>%
  summarise(
    wt_control_rows = n(),
    wt_control_lob_positive = sum(detected_above_LoB, na.rm = TRUE),
    wt_control_lod_positive = sum(detected_above_LoD, na.rm = TRUE),
    median_fractional_abundance = median(fractional_abundance, na.rm = TRUE),
    .groups = "drop"
  )

positive_control_recovery <- sample_results %>%
  filter(sample_type == "positive_control") %>%
  group_by(method_id, method_family, assay) %>%
  summarise(
    positive_control_rows = n(),
    positive_control_lob_positive = sum(detected_above_LoB, na.rm = TRUE),
    positive_control_lod_positive = sum(detected_above_LoD, na.rm = TRUE),
    median_fractional_abundance = median(fractional_abundance, na.rm = TRUE),
    .groups = "drop"
  )

e2e_files <- list.files(path_v2("tables"), pattern = "_e2e_checks[.]csv$", full.names = TRUE)
e2e_files <- e2e_files[
  !basename(e2e_files) %in% c(
    "method_report_e2e_checks.csv",
    "method_svg_panel_e2e_checks.csv"
  )
]
e2e_summary <- bind_rows(lapply(e2e_files, function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(source = basename(path), .before = 1)
}))
package_versions <- package_versions_table()

write_csv(counts, path_v2("tables", "method_well_counts.csv"))
write_csv(sample_results, path_v2("tables", "method_sample_region_results.csv"))
write_csv(pass_matrix, path_v2("tables", "method_lob_lod_pass_matrix.csv"))
write_csv(positive_rows, path_v2("tables", "method_positive_rows.csv"))
write_csv(method_summary, path_v2("tables", "method_summary.csv"))
write_csv(control_false_positive, path_v2("tables", "method_control_false_positive_summary.csv"))
write_csv(positive_control_recovery, path_v2("tables", "method_positive_control_recovery_summary.csv"))
write_csv(e2e_summary, path_v2("tables", "method_e2e_check_summary.csv"))
write_csv(package_versions, path_v2("tables", "method_package_versions.csv"))

heatmap_plot <- plot_pass_heatmap(pass_matrix)
failure_plot <- plot_failure_rates(method_summary)
gating_output <- write_representative_gating_plots(read_plot_droplets())
summary_plot_manifest <- bind_rows(
  write_plot_pair(heatmap_plot, "method_lob_lod_heatmap", width = 13, height = 4.2),
  write_plot_pair(failure_plot, "method_failure_rates", width = 13, height = 4.2)
)
write_csv(gating_output$manifest, path_v2("tables", "method_gating_plot_manifest.csv"))
plot_manifest <- bind_rows(summary_plot_manifest, gating_output$manifest)
write_csv(plot_manifest, path_v2("tables", "method_summary_plot_manifest.csv"))
panel_pdf <- write_panel_pdf(heatmap_plot, failure_plot)
report_pdf <- write_report_pdf(
  method_summary,
  pass_matrix,
  positive_rows,
  e2e_summary,
  package_versions,
  heatmap_plot,
  failure_plot,
  gating_output$plots
)

report_checks <- tibble(
  check = c(
    "method_summary_has_all_methods",
    "polygon_method_present",
    "positive_rows_table_written",
    "control_tables_written",
    "package_versions_written",
    "individual_plots_exist",
    "representative_gating_plots_exist",
    "panel_pdf_exists",
    "report_pdf_exists",
    "all_e2e_checks_pass"
  ),
  passed = c(
    nrow(method_summary) == length(unique(counts$method_id)),
    "polygon_control_hull_gates" %in% method_summary$method_id,
    file.exists(path_v2("tables", "method_positive_rows.csv")),
    file.exists(path_v2("tables", "method_control_false_positive_summary.csv")) &&
      file.exists(path_v2("tables", "method_positive_control_recovery_summary.csv")),
    file.exists(path_v2("tables", "method_package_versions.csv")),
    all(file.exists(c(plot_manifest$svg_path, plot_manifest$pdf_path))),
    all(file.exists(c(gating_output$manifest$svg_path, gating_output$manifest$pdf_path))) &&
      all(file.info(c(gating_output$manifest$svg_path, gating_output$manifest$pdf_path))$size > 0),
    file.exists(panel_pdf) && file.info(panel_pdf)$size > 0,
    file.exists(report_pdf) && file.info(report_pdf)$size > 0,
    all(e2e_summary$passed)
  )
)
write_csv(report_checks, path_v2("tables", "method_report_e2e_checks.csv"))
stopifnot(all(report_checks$passed))

cat("method_count_rows=", nrow(counts), "\n", sep = "")
cat("sample_region_rows=", nrow(sample_results), "\n", sep = "")
cat("method_summary_rows=", nrow(method_summary), "\n", sep = "")
cat("positive_rows=", nrow(positive_rows), "\n", sep = "")
cat("report_pdf=", report_pdf, "\n", sep = "")
