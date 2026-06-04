library(tidyverse)
library(readxl)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
scratch_root <- file.path(project_root, "scratch", "ddpcr-new")
analysis_root <- file.path(scratch_root, "analysis")
data_path <- file.path(analysis_root, "results", "ddPCR", "SNV_data_final.xlsx")
out_dir <- file.path(scratch_root, "figures", "ddpcr", "main_ddpcr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

snv_data <- readxl::read_excel(data_path)

lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)
mutation_list <- c("D178N", "E200K", "P102L")
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

make_plot <- function(mut) {
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

  ggplot(df, aes(x = participant, y = fractional_abundance, colour = brain_region)) +
    geom_point(aes(shape = detected_above_LoB, size = is_pooled), position = position_dodge(width = 0.6), na.rm = TRUE) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high, group = brain_region), position = position_dodge(width = 0.6), width = 0, linewidth = 0.3, na.rm = TRUE) +
    geom_hline(yintercept = lod_cut[[mut]], linetype = "dashed", linewidth = 0.5) +
    coord_cartesian(ylim = c(0, y_max)) +
    scale_colour_manual(name = "Brain region", values = region_colours, labels = region_labels, drop = FALSE) +
    scale_shape_manual(name = "Above LoB", values = c(`FALSE` = 1, `TRUE` = 16), labels = c("No", "Yes")) +
    scale_size_manual(name = "Pooled row", values = c(`FALSE` = 1.8, `TRUE` = 3.0), labels = c("No", "Yes")) +
    labs(
      title = paste(mut, "corrected sample-region ddPCR"),
      subtitle = "Pooled rows are marked larger; y-axis expands if corrected pooled values exceed the old range",
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

plots <- setNames(lapply(mutation_list, make_plot), mutation_list)

for (mut in mutation_list) {
  ggsave(
    filename = file.path(out_dir, paste0("SNV_", mut, "_panel_corrected.pdf")),
    plot = plots[[mut]],
    width = 11,
    height = 4.5,
    dpi = 300
  )
}

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

lod_df <- tibble(
  mutation = factor(names(lod_cut), levels = mutation_list),
  LoD = as.numeric(lod_cut)
)

combined_plot <- ggplot(
  combined_df,
  aes(x = participant, y = fractional_abundance, colour = brain_region)
) +
  geom_point(aes(shape = detected_above_LoB, size = is_pooled),
             position = position_dodge(width = 0.6), na.rm = TRUE) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, group = brain_region),
                position = position_dodge(width = 0.6), width = 0, linewidth = 0.3, na.rm = TRUE) +
  geom_hline(data = lod_df, aes(yintercept = LoD), inherit.aes = FALSE,
             linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~mutation, ncol = 1, scales = "free_y") +
  scale_colour_manual(name = "Brain region", values = region_colours, labels = region_labels, drop = FALSE) +
  scale_shape_manual(name = "Above LoB", values = c(`FALSE` = 1, `TRUE` = 16), labels = c("No", "Yes")) +
  scale_size_manual(name = "Pooled row", values = c(`FALSE` = 1.8, `TRUE` = 3.0), labels = c("No", "Yes")) +
  labs(
    title = "Corrected sample-region ddPCR",
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

ggsave(
  filename = file.path(out_dir, "SNV_all_mutations_corrected_single_page.pdf"),
  plot = combined_plot,
  width = 11,
  height = 12,
  dpi = 300
)

write_multipage_pdf <- function(path, plots, width, height) {
  pdf(path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  for (plot in plots) {
    print(plot)
  }
  invisible(TRUE)
}

write_multipage_pdf(
  file.path(out_dir, "SNV_all_mutations_corrected_multipage.pdf"),
  plots,
  width = 11,
  height = 4.5
)

affected_rows <- snv_data %>%
  mutate(is_pooled = as.logical(is_pooled)) %>%
  filter(is_pooled) %>%
  arrange(mutation, participant, brain_region)

readr::write_csv(
  affected_rows,
  file.path(out_dir, "pooled_rows_in_main_ddpcr_corrected.csv")
)
