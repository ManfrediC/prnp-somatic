library(tidyverse)
library(readxl)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_path <- file.path(project_root, "results", "ddPCR", "SNV_pooled_participant.xlsx")
out_dir <- file.path(project_root, "manuscript", "figures", "ddpcr_fractional_abundance_pooled")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pooled_plot_data <- readxl::read_excel(data_path)

mutation_list <- c("D178N", "E200K", "P102L")
lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)

participant_levels <- function(x) {
  intersect(c(paste0("CJD", 1:200), paste0("Control", 1:200)), unique(x))
}

make_plot <- function(mut) {
  df <- pooled_plot_data %>%
    filter(mutation == mut) %>%
    drop_na(fractional_abundance) %>%
    mutate(
      participant = factor(participant, levels = participant_levels(participant)),
      detected_above_lob = as.logical(detected_above_lob),
      detected_above_lod = as.logical(detected_above_lod)
    )

  y_max <- max(c(0.25, df$ci_high, lod_cut[[mut]]), na.rm = TRUE) * 1.08

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

plots <- setNames(lapply(mutation_list, make_plot), mutation_list)

for (mut in mutation_list) {
  ggsave(
    filename = file.path(out_dir, paste0("SNV_pooled_", mut, "_panel.pdf")),
    plot = plots[[mut]],
    width = 9,
    height = 4.5,
    dpi = 300
  )
}

combined_df <- pooled_plot_data %>%
  filter(mutation %in% mutation_list) %>%
  drop_na(fractional_abundance) %>%
  mutate(
    mutation = factor(mutation, levels = mutation_list),
    participant = factor(participant, levels = participant_levels(participant)),
    detected_above_lob = as.logical(detected_above_lob),
    detected_above_lod = as.logical(detected_above_lod)
  )

lod_df <- tibble(
  mutation = factor(names(lod_cut), levels = mutation_list),
  LoD = as.numeric(lod_cut)
)

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

ggsave(
  filename = file.path(out_dir, "SNV_pooled_all_mutations.pdf"),
  plot = combined_plot,
  width = 9,
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
  file.path(out_dir, "SNV_pooled_all_mutations_multipage.pdf"),
  plots,
  width = 9,
  height = 4.5
)

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
