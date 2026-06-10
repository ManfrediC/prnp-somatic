library(tidyverse)
library(readxl)
library(cowplot)

# ---- paths and inputs ----

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_path <- file.path(project_root, "results", "ddPCR", "SNV_data_final.xlsx")
out_dir <- file.path(project_root, "manuscript", "figures", "ddpcr_fractional_abundance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read the sample-region ddPCR table produced by the main ddPCR workflow.
snv_data <- readxl::read_excel(data_path)

# ---- plotting constants ----

# Keep mutation ordering, LoD cut-offs, and requested panels explicit so single-panel
# reruns remain reproducible from the environment alone.
lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)
mutation_list <- c("D178N", "E200K", "P102L")
max_y <- 0.4
target_mutation_panels <- strsplit(
  Sys.getenv("DDPCR_TARGET_MUTATION_PANELS", paste(mutation_list, collapse = ",")),
  ","
)[[1]] %>%
  trimws()
target_mutation_panels <- intersect(target_mutation_panels, mutation_list)
if (length(target_mutation_panels) == 0L) {
  # Fall back to the full figure set when the environment variable is empty or invalid.
  target_mutation_panels <- mutation_list
}

# Display labels and colours are centralised so individual and combined panels
# use the same region vocabulary.
region_labels <- c(
  bg = "basal ganglia",
  cb = "cerebellum",
  fr = "frontal cortex",
  hc = "hippocampus",
  ps = "pons",
  sn = "midbrain",
  th = "thalamus"
)
region_colours <- c(
  bg = "#E69F00",
  cb = "#56B4E9",
  fr = "#009E73",
  hc = "#0072B2",
  ps = "#999999",
  sn = "#D55E00",
  th = "#CC79A7"
)

# Shared dodge used by all point / CI layers so points and their intervals line up.
point_dodge <- position_dodge(width = 0.6)

# Build one mutation-specific sample-region plot.
make_plot <- function(mut) {
  # Restrict to plotted rows and keep the historical CJD30 E200K exclusion
  # aligned with the rest of the ddPCR figure workflow.
  df <- snv_data %>%
    filter(mutation == mut) %>%
    drop_na(fractional_abundance) %>%
    filter(mutation != "E200K" | participant != "CJD30") %>%
    mutate(
      participant = factor(
        participant,
        levels = intersect(c(paste0("CJD", 1:200), paste0("Control", 1:200)), unique(participant))
      ),
      brain_region = factor(brain_region, levels = names(region_labels)),
      detected_above_LoB = as.logical(detected_above_LoB),
      is_pooled = as.logical(is_pooled)
    )

  y_max <- max(c(0.25, df$ci_high, lod_cut[[mut]]), na.rm = TRUE) * 1.08

  # Draw point estimates, confidence intervals, and the assay-specific LoD line
  # on the common percent fractional-abundance scale.
  ggplot(df, aes(x = participant, y = fractional_abundance, colour = brain_region)) +
    geom_point(
      aes(shape = detected_above_LoB, size = is_pooled, group = brain_region),
      position = point_dodge, na.rm = TRUE
    ) +
    geom_errorbar(
      aes(ymin = ci_low, ymax = ci_high, group = brain_region),
      position = point_dodge, width = 0, linewidth = 0.3, na.rm = TRUE
    ) +
    geom_hline(yintercept = lod_cut[[mut]], linetype = "dashed", linewidth = 0.5) +
    coord_cartesian(ylim = c(0, y_max)) +
    scale_colour_manual(name = "Brain region", values = region_colours, labels = region_labels, drop = FALSE) +
    scale_shape_manual(name = "Above LoB", values = c(`FALSE` = 1, `TRUE` = 16), labels = c("No", "Yes")) +
    scale_size_manual(name = "Pooled row", values = c(`FALSE` = 1.8, `TRUE` = 3.0), labels = c("No", "Yes")) +
    labs(
      title = paste(mut, "sample-region ddPCR"),
      subtitle = "Pooled rows are marked larger; y-axis expands if pooled values exceed the old range",
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

# Prepare all mutation plots once so panel selection only controls which files
# are written.
plots <- setNames(lapply(mutation_list, make_plot), mutation_list)

# Write requested single-mutation panels.
for (mut in target_mutation_panels) {
  ggsave(
    filename = file.path(out_dir, paste0("SNV_", mut, "_panel.pdf")),
    plot = plots[[mut]],
    width = 11,
    height = 4.5,
    dpi = 300
  )
}

# ---- combined mutation panel ----

# Build the combined long table with fixed mutation and region ordering.
combined_df <- snv_data %>%
  filter(mutation %in% mutation_list) %>%
  drop_na(fractional_abundance) %>%
  filter(mutation != "E200K" | participant != "CJD30") %>%
  mutate(
    mutation = factor(mutation, levels = mutation_list),
    participant = factor(
      participant,
      levels = intersect(c(paste0("CJD", 1:200), paste0("Control", 1:200)), unique(participant))
    ),
    brain_region = factor(brain_region, levels = names(region_labels)),
    detected_above_LoB = as.logical(detected_above_LoB),
    is_pooled = as.logical(is_pooled)
  )

# Store LoD cut-offs in a facetable data frame for per-panel horizontal lines.
lod_df <- tibble(
  mutation = factor(names(lod_cut), levels = mutation_list),
  LoD = as.numeric(lod_cut)
)

# Stack all mutations into one manuscript-facing figure.
combined_plot <- ggplot(
  combined_df,
  aes(x = participant, y = fractional_abundance, colour = brain_region)
) +
  geom_point(
    aes(shape = detected_above_LoB, size = is_pooled, group = brain_region),
    position = point_dodge,
    na.rm = TRUE
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high, group = brain_region),
    position = point_dodge,
    width = 0,
    linewidth = 0.3,
    na.rm = TRUE
  ) +
  geom_hline(data = lod_df, aes(yintercept = LoD), inherit.aes = FALSE,
             linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~mutation, ncol = 1) +
  coord_cartesian(ylim = c(0, max_y)) +
  scale_colour_manual(name = "Brain region", values = region_colours, labels = region_labels, drop = FALSE) +
  scale_shape_manual(name = "Above LoB", values = c(`FALSE` = 1, `TRUE` = 16), labels = c("No", "Yes")) +
  scale_size_manual(name = "Pooled row", values = c(`FALSE` = 1.8, `TRUE` = 3.0), labels = c("No", "Yes")) +
  labs(
    title = "Sample-region ddPCR",
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
  filename = file.path(out_dir, "SNV_all_mutations.pdf"),
  plot = combined_plot,
  width = 11,
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
legend_bottom_plot <- cowplot::plot_grid(
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

# Save the final combined panel as PDF and SVG.
ggsave(
  filename = file.path(out_dir, "SNV_all_mutations_legend_bottom_final.pdf"),
  plot = legend_bottom_plot,
  width = 11,
  height = 12,
  dpi = 300
)

save_svg_if_available(
  plot = legend_bottom_plot,
  path = file.path(out_dir, "SNV_all_mutations_legend_bottom_final.svg"),
  width = 11,
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
  file.path(out_dir, "SNV_all_mutations_multipage.pdf"),
  plots,
  width = 11,
  height = 4.5
)

# Record pooled rows that affect the sample-region figure so reviewers can audit
# the larger point markers separately from the plot.
affected_rows <- snv_data %>%
  mutate(is_pooled = as.logical(is_pooled)) %>%
  filter(is_pooled) %>%
  arrange(mutation, participant, brain_region)

readr::write_csv(
  affected_rows,
  file.path(out_dir, "pooled_rows_in_main_ddpcr.csv")
)
