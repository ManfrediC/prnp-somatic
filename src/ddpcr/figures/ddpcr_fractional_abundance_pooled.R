library(tidyverse)
library(readxl)
library(cowplot)

# ---- paths and inputs ----

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_path <- file.path(project_root, "results", "ddPCR", "SNV_pooled_participant.xlsx")
out_dir <- file.path(project_root, "manuscript", "figures", "ddpcr_fractional_abundance_pooled")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read the participant-pooled ddPCR table produced by the main workflow.
pooled_plot_data <- readxl::read_excel(data_path)

# ---- plotting constants ----

# Keep mutation order and LoD cut-offs explicit for consistent single-panel and
# combined-panel output.
mutation_list <- c("D178N", "E200K", "P102L")
lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)
target_mutation_panels <- strsplit(
  Sys.getenv("DDPCR_TARGET_MUTATION_PANELS", paste(mutation_list, collapse = ",")),
  ","
)[[1]] %>%
  trimws()
target_mutation_panels <- intersect(target_mutation_panels, mutation_list)
if (length(target_mutation_panels) == 0L) {
  # Fall back to every mutation when the optional environment filter is absent.
  target_mutation_panels <- mutation_list
}

# Preserve natural participant ordering while dropping unused IDs.
participant_levels <- function(x) {
  intersect(c(paste0("CJD", 1:200), paste0("Control", 1:200)), unique(x))
}

# Build one mutation-specific participant-pooled plot.
make_plot <- function(mut) {
  # Coerce workbook flags before mapping them to point shape/fill.
  df <- pooled_plot_data %>%
    filter(mutation == mut) %>%
    drop_na(fractional_abundance) %>%
    mutate(
      participant = factor(participant, levels = participant_levels(participant)),
      detected_above_lob = as.logical(detected_above_lob),
      detected_above_lod = as.logical(detected_above_lod)
    )

  y_max <- max(c(0.25, df$ci_high, lod_cut[[mut]]), na.rm = TRUE) * 1.08

  # Draw participant-level point estimates, confidence intervals, and the
  # assay-specific LoD line on the percent fractional-abundance scale.
  ggplot(df, aes(x = participant, y = fractional_abundance)) +
    geom_point(aes(shape = detected_above_lob, fill = detected_above_lod),
               colour = "black", size = 2.2, stroke = 0.4, na.rm = TRUE) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0, linewidth = 0.3, na.rm = TRUE) +
    geom_hline(yintercept = lod_cut[[mut]], linetype = "dashed", linewidth = 0.5) +
    coord_cartesian(ylim = c(0, y_max)) +
    scale_shape_manual(
      name = "Above LoB",
      values = c(`FALSE` = 1, `TRUE` = 21),
      labels = c("No", "Yes")
    ) +
    scale_fill_manual(
      name = "Above LoD",
      values = c(`FALSE` = "white", `TRUE` = "#0072B2"),
      labels = c("No", "Yes")
    ) +
    labs(
      title = paste(mut, "participant-pooled ddPCR"),
      subtitle = "Pooled fractional abundance, confidence intervals, and LoB are on the same percent scale as LoD",
      x = "Participant",
      y = paste0("Fractional abundance of ", mut, " allele (%)")
    ) +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

# Prepare all mutation plots once so file writing can be filtered separately.
plots <- setNames(lapply(mutation_list, make_plot), mutation_list)

# Write requested single-mutation participant-pooled panels.
for (mut in target_mutation_panels) {
  ggsave(
    filename = file.path(out_dir, paste0("SNV_pooled_", mut, "_panel.pdf")),
    plot = plots[[mut]],
    width = 9,
    height = 4.5,
    dpi = 300
  )
}

# ---- combined mutation panel ----

# Build the combined table with stable mutation and participant ordering.
combined_df <- pooled_plot_data %>%
  filter(mutation %in% mutation_list) %>%
  drop_na(fractional_abundance) %>%
  mutate(
    mutation = factor(mutation, levels = mutation_list),
    participant = factor(participant, levels = participant_levels(participant)),
    detected_above_lob = as.logical(detected_above_lob),
    detected_above_lod = as.logical(detected_above_lod)
  )

# Store LoD cut-offs in a facetable data frame for per-panel horizontal lines.
lod_df <- tibble(
  mutation = factor(names(lod_cut), levels = mutation_list),
  LoD = as.numeric(lod_cut)
)

# Stack all participant-pooled mutation panels into one figure.
combined_plot <- ggplot(combined_df, aes(x = participant, y = fractional_abundance)) +
  geom_point(aes(shape = detected_above_lob, fill = detected_above_lod),
             colour = "black", size = 2.2, stroke = 0.4, na.rm = TRUE) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, linewidth = 0.3, na.rm = TRUE) +
  geom_hline(data = lod_df, aes(yintercept = LoD), inherit.aes = FALSE,
             linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~mutation, ncol = 1, scales = "free_y") +
  scale_shape_manual(
    name = "Above LoB",
    values = c(`FALSE` = 1, `TRUE` = 21),
    labels = c("No", "Yes")
  ) +
  scale_fill_manual(
    name = "Above LoD",
    values = c(`FALSE` = "white", `TRUE` = "#0072B2"),
    labels = c("No", "Yes")
  ) +
  labs(
    title = "Participant-pooled ddPCR",
    subtitle = "Pooled fractional abundance, confidence intervals, and LoB are on the same percent scale as LoD",
    x = "Participant",
    y = "Fractional abundance (%)"
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# Save the combined panel with the default right-side legend.
ggsave(
  filename = file.path(out_dir, "SNV_pooled_all_mutations.pdf"),
  plot = combined_plot,
  width = 9,
  height = 12,
  dpi = 300
)

# Write SVGs when svglite is available, otherwise remove stale vectors so
# missing optional tooling is visible rather than silently preserving old output.
save_svg_if_available <- function(plot, path, width, height) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      filename = path,
      plot = plot,
      width = width,
      height = height,
      dpi = 300
    )
  } else {
    unlink(path, force = TRUE)
    warning("svglite is not installed; removed stale SVG output: ", path, call. = FALSE)
  }
}

# Move the combined legend below the facets for the final manuscript layout.
one_legend_plot <- cowplot::plot_grid(
  combined_plot + theme(legend.position = "none"),
  cowplot::get_legend(
    combined_plot +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal"
      )
  ),
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# Save the final one-legend combined panel as PDF and SVG.
ggsave(
  filename = file.path(out_dir, "SNV_pooled_all_mutations_one_legend.pdf"),
  plot = one_legend_plot,
  width = 9,
  height = 12,
  dpi = 300
)

save_svg_if_available(
  plot = one_legend_plot,
  path = file.path(out_dir, "SNV_pooled_all_mutations_one_legend.svg"),
  width = 9,
  height = 12
)

# Write one plot per page for easier review of the individual mutation panels.
write_multipage_pdf <- function(path, plots, width, height) {
  pdf(path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  for (plot in plots) {
    print(plot)
  }
  invisible(TRUE)
}

write_multipage_pdf(
  file.path(out_dir, "SNV_pooled_all_mutations_multipage.pdf"),
  plots,
  width = 9,
  height = 4.5
)

# Record LoB+LoD passing participant-pooled rows for audit beside the figure.
pass_rows <- pooled_plot_data %>%
  mutate(
    detected_above_lob = as.logical(detected_above_lob),
    detected_above_lod = as.logical(detected_above_lod)
  ) %>%
  filter(detected_above_lob, detected_above_lod) %>%
  arrange(mutation, participant)

readr::write_csv(
  pooled_plot_data,
  file.path(out_dir, "participant_pooled_rows.csv")
)

readr::write_csv(
  pass_rows,
  file.path(out_dir, "participant_pooled_lob_lod_pass_rows.csv")
)
