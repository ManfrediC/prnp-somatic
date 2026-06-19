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
               colour = "black", size = 2.4, stroke = 0.45, na.rm = TRUE) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0, linewidth = 0.3, na.rm = TRUE) +
    geom_hline(yintercept = lod_cut[[mut]], linetype = "dashed", linewidth = 0.5) +
    coord_cartesian(ylim = c(0, y_max)) +
    scale_shape_manual(
      name = "Above LoB",
      values = c(`FALSE` = 1, `TRUE` = 21),
      labels = c("No", "Yes"),
      drop = FALSE
    ) +
    scale_fill_manual(
      name = "Above LoD",
      values = c(`FALSE` = "white", `TRUE` = "#009E73"),
      labels = c("No", "Yes"),
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(fill = "white", size = 3)),
      fill = guide_legend(order = 2, override.aes = list(shape = 21, size = 3))
    ) +
    labs(
      title = paste(mut, "participant-pooled ddPCR"),
      x = "Participant",
      y = paste0("Fractional abundance of ", mut, " allele (%)")
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 9.5, angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 2)),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12, margin = margin(t = 8)),
      axis.title.y = element_text(size = 12, margin = margin(r = 6)),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.spacing.y = grid::unit(5, "pt")
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

mutation_label_df <- combined_df %>%
  group_by(mutation) %>%
  summarise(
    participant = factor(
      tail(levels(combined_df$participant), 1),
      levels = levels(combined_df$participant)
    ),
    label_y = max(
      c(0.25, ci_high, lod_cut[as.character(first(mutation))]),
      na.rm = TRUE
    ),
    label = as.character(first(mutation)),
    .groups = "drop"
  )

# Stack all participant-pooled mutation panels into one figure.
combined_plot <- ggplot(combined_df, aes(x = participant, y = fractional_abundance)) +
  geom_point(aes(shape = detected_above_lob, fill = detected_above_lod),
             colour = "black", size = 2.4, stroke = 0.45, na.rm = TRUE) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, linewidth = 0.3, na.rm = TRUE) +
  geom_hline(data = lod_df, aes(yintercept = LoD),
             linetype = "dashed", linewidth = 0.5) +
  geom_text(
    data = mutation_label_df,
    aes(x = participant, y = label_y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.15,
    family = "Helvetica",
    fontface = "bold",
    size = 10.5,
    colour = "black",
    show.legend = FALSE
  ) +
  facet_wrap(~mutation, ncol = 1, scales = "free_y") +
  scale_shape_manual(
    name = "Above LoB",
    values = c(`FALSE` = 1, `TRUE` = 21),
    labels = c("No", "Yes"),
    drop = FALSE
  ) +
  scale_fill_manual(
    name = "Above LoD",
    values = c(`FALSE` = "white", `TRUE` = "#009E73"),
    labels = c("No", "Yes"),
    drop = FALSE
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(shape = c(1, 16), fill = c("white", "black"), size = 3)
    ),
    fill = guide_legend(
      order = 2,
      override.aes = list(shape = 21, fill = c("white", "#009E73"), size = 3)
    )
  ) +
  labs(
    title = "Participant-pooled ddPCR",
    x = "Participant",
    y = "Fractional abundance (%)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(size = 16),
    axis.text.x = element_text(size = 9.5, angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 2)),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, margin = margin(r = 6)),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.spacing.y = grid::unit(5, "pt"),
    legend.spacing.x = grid::unit(28, "pt")
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

fix_svg_mutation_label_font <- function(path, labels) {
  if (!file.exists(path)) {
    return(invisible(FALSE))
  }

  svg_lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  for (label in labels) {
    label_line <- grepl(paste0(">", label, "</text>"), svg_lines, fixed = TRUE)
    svg_lines[label_line] <- sub(
      'font-family: "Arial";',
      'font-family: "Helvetica";',
      svg_lines[label_line],
      fixed = TRUE
    )
  }
  writeLines(svg_lines, path, useBytes = TRUE)
  invisible(TRUE)
}

# Move the combined legend below the facets for the final manuscript layout.
one_legend_plot <- cowplot::plot_grid(
  combined_plot + theme(legend.position = "none"),
  cowplot::get_legend(
    combined_plot +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.spacing.x = grid::unit(34, "pt"),
        legend.box.spacing = grid::unit(8, "pt")
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

one_legend_svg_path <- file.path(out_dir, "SNV_pooled_all_mutations_one_legend.svg")

save_svg_if_available(
  plot = one_legend_plot,
  path = one_legend_svg_path,
  width = 9,
  height = 12
)
fix_svg_mutation_label_font(one_legend_svg_path, mutation_list)

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
