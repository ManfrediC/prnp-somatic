# Well-level E200K ddPCR results for the two positive samples
# (CJD4 pons, CJD21 thalamus): figure + table data.
#
# Inputs (existing repository artefacts, no recomputation from raw archives):
#   results/ddPCR/scatterplots/lob_lod_positive/selected_wells.csv
#     per-well FA/CI/LoB/LoD calls written by create_ddpcr_scatterplots.R
#   raw/ddpcr/manifests/sample_manifest.csv
#     per-well accepted/positive droplet counts
#   results/ddPCR/SNV_data_final.xlsx
#     pooled sample-level results written by create_snv_dataframe.R
#
# Outputs:
#   manuscript/figures/ddpcr_e200k_positive_well_fa/ddpcr_e200k_positive_well_fa.{pdf,svg}
#   manuscript/tables/ddpcr_e200k_positive_well_results/ddpcr_e200k_positive_well_results.csv
#   manuscript/tables/ddpcr_e200k_positive_well_results/ddpcr_e200k_positive_well_results_rows.tex

library(tidyverse)
library(readxl)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
fig_dir <- file.path(project_root, "manuscript", "figures", "ddpcr_e200k_positive_well_fa")
tab_dir <- file.path(project_root, "manuscript", "tables", "ddpcr_e200k_positive_well_results")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

lod_e200k <- 0.067  # E200K assay limit of detection (% FA), as in create_snv_dataframe.R

# Region display vocabulary shared with ddpcr_fractional_abundance.R.
region_labels <- c(ps = "pons", th = "thalamus")
region_colours <- c(ps = "#999999", th = "#CC79A7")

# ---- per-well data ---------------------------------------------------------

# Well-level FA/CI/LoB/LoD calls for the five E200K-positive wells.
# Rename the pooled sample-level columns first: they are kept for the
# cross-check against SNV_data_final.xlsx and must not shadow the well-level
# ci_low/ci_high created in the transmute below.
wells <- read_csv(
  file.path(project_root, "results", "ddPCR", "scatterplots", "lob_lod_positive", "selected_wells.csv"),
  show_col_types = FALSE
) %>%
  rename(pooled_fa = fractional_abundance, pooled_ci_low = ci_low, pooled_ci_high = ci_high) %>%
  transmute(
    participant = participant_display,
    brain_region,
    sample = paste(participant_display, brain_region_display),
    replicate = as.character(replicate_index),
    run_id, run_date, well,
    fa = individual_fractional_abundance,
    ci_low = individual_ci_low,
    ci_high = individual_ci_high,
    lob_fa = individual_lob_fa,
    above_lob = detected_above_LoB,
    above_lod = detected_above_LoD,
    pooled_fa, pooled_ci_low, pooled_ci_high
  )

# Accepted/mutant-positive droplet counts per well from the raw manifest
# (the mutant-target row of each well carries the totals for both channels).
well_counts <- read_csv(
  file.path(project_root, "raw", "ddpcr", "manifests", "sample_manifest.csv"),
  show_col_types = FALSE
) %>%
  filter(target_role == "mutant") %>%
  select(run_id, well, accepted = accepted_droplets, positives)

wells <- wells %>% left_join(well_counts, by = c("run_id", "well"))

# ---- pooled sample-level data ------------------------------------------------

pooled <- read_excel(file.path(project_root, "results", "ddPCR", "SNV_data_final.xlsx")) %>%
  filter(mutation == "E200K", paste(participant, brain_region) %in% c("CJD4 ps", "CJD21 th")) %>%
  transmute(
    participant,
    brain_region,
    sample = paste(participant, recode(brain_region, !!!region_labels)),
    replicate = "Pooled",
    run_date = as.Date(NA),
    well = NA_character_,
    accepted = n_total_droplets,
    positives = n_mut_droplets,
    fa = fractional_abundance,
    ci_low, ci_high, lob_fa,
    above_lob = detected_above_LoB,
    above_lod = detected_above_LoD
  )

# ---- cross-checks ------------------------------------------------------------

# Row counts and complete droplet data.
stopifnot(nrow(wells) == 5L, nrow(pooled) == 2L, !any(is.na(wells$accepted)))

# Per-well droplet counts must sum to the pooled sample totals.
sums_ok <- wells %>%
  group_by(sample) %>%
  summarise(well_accepted = sum(accepted), well_positives = sum(positives), .groups = "drop") %>%
  inner_join(pooled, by = "sample") %>%
  summarise(ok = all(well_accepted == accepted & well_positives == positives)) %>%
  pull(ok)
stopifnot(sums_ok)

# The pooled estimates carried in selected_wells.csv must match SNV_data_final.xlsx.
pooled_ok <- wells %>%
  distinct(sample, pooled_fa, pooled_ci_low, pooled_ci_high) %>%
  inner_join(pooled, by = "sample") %>%
  summarise(ok = all(near(pooled_fa, fa) & near(pooled_ci_low, ci_low) & near(pooled_ci_high, ci_high))) %>%
  pull(ok)
stopifnot(pooled_ok)

wells <- wells %>% select(-pooled_fa, -pooled_ci_low, -pooled_ci_high)

# ---- figure -------------------------------------------------------------------

# One row per measurement: five wells plus the two pooled sample summaries.
measurements <- bind_rows(
  wells %>% mutate(measurement = if_else(above_lob, "Well (above LoB)", "Well (below LoB)")),
  pooled %>% mutate(measurement = "Pooled")
) %>%
  mutate(
    sample = factor(sample, levels = c("CJD4 pons", "CJD21 thalamus")),
    replicate = factor(replicate, levels = c("1", "2", "3", "Pooled")),
    measurement = factor(measurement, levels = c("Well (below LoB)", "Well (above LoB)", "Pooled"))
  ) %>%
  arrange(sample, replicate) %>%
  # Numeric x position per panel for the LoB tick segments
  # (free_x facets drop unused replicate levels, so index per sample).
  group_by(sample) %>%
  mutate(x = match(as.character(replicate), unique(as.character(replicate)))) %>%
  ungroup()

y_max <- max(measurements$ci_high, measurements$lob_fa, lod_e200k) * 1.08

# Point estimates, confidence intervals, per-measurement LoB ticks and the
# assay LoD line, following ddpcr_fractional_abundance.R styling.
fa_plot <- ggplot(measurements, aes(x = replicate, y = fa)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, linewidth = 0.3) +
  geom_segment(
    aes(x = x - 0.18, xend = x + 0.18, y = lob_fa, yend = lob_fa, linetype = "LoB"),
    linewidth = 0.4, colour = "grey30"
  ) +
  geom_hline(aes(yintercept = lod_e200k, linetype = "LoD (E200K)"), linewidth = 0.5) +
  geom_point(aes(colour = brain_region, shape = measurement), size = 2.0) +
  facet_grid(. ~ sample, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(0, y_max)) +
  scale_colour_manual(name = "Brain region", values = region_colours, labels = region_labels) +
  scale_shape_manual(name = "Measurement", values = c("Well (below LoB)" = 1, "Well (above LoB)" = 16, "Pooled" = 18)) +
  scale_linetype_manual(
    name = NULL,
    values = c("LoB" = "solid", "LoD (E200K)" = "dashed"),
    guide = guide_legend(override.aes = list(colour = c("grey30", "black")))
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(shape = 16, size = 3)),
    shape = guide_legend(order = 2, override.aes = list(size = 3)),
    linetype = guide_legend(order = 3)
  ) +
  labs(x = "Replicate", y = "Fractional abundance of E200K allele (%)") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 9.5),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, margin = margin(r = 6)),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

ggsave(
  filename = file.path(fig_dir, "ddpcr_e200k_positive_well_fa.pdf"),
  plot = fa_plot, device = grDevices::pdf, width = 7.2, height = 4.0, dpi = 300
)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = file.path(fig_dir, "ddpcr_e200k_positive_well_fa.svg"),
    plot = fa_plot, width = 7.2, height = 4.0, device = svglite::svglite, dpi = 300
  )
} else {
  ggsave(
    filename = file.path(fig_dir, "ddpcr_e200k_positive_well_fa.svg"),
    plot = fa_plot, width = 7.2, height = 4.0, device = grDevices::svg, dpi = 300
  )
}

# ---- table data: CSV + lean TeX rows ------------------------------------------

# Machine-readable table (5 well rows + 2 pooled rows).
table_data <- measurements %>%
  transmute(
    participant,
    brain_region = recode(brain_region, !!!region_labels),
    measurement = if_else(replicate == "Pooled", "pooled", "well"),
    replicate_index = if_else(replicate == "Pooled", NA_integer_, as.integer(as.character(replicate))),
    accepted_droplets = accepted,
    mut_positive_droplets = positives,
    fa_pct = round(fa, 3),
    ci_low_pct = round(ci_low, 3),
    ci_high_pct = round(ci_high, 3),
    lob_fa_pct = round(lob_fa, 3),
    above_lob, above_lod
  )
write_csv(table_data, file.path(tab_dir, "ddpcr_e200k_positive_well_results.csv"))

# Lean TeX rows: one `a & b & ... \\` line per row, with a section marker
# before the pooled block. Formatting is added by the styling step.
format_row <- function(row, replicate_label) {
  paste0(
    paste(
      row$participant, row$brain_region, replicate_label,
      row$accepted_droplets, row$mut_positive_droplets,
      sprintf("%.3f", row$fa_pct),
      sprintf("%.3f--%.3f", row$ci_low_pct, row$ci_high_pct),
      sprintf("%.3f", row$lob_fa_pct),
      if_else(row$above_lob, "Yes", "No"),
      if_else(row$above_lod, "Yes", "No"),
      sep = " & "
    ),
    " \\\\"
  )
}

well_data <- table_data %>% filter(measurement == "well")
pooled_data <- table_data %>% filter(measurement == "pooled")

well_rows <- vapply(
  seq_len(nrow(well_data)),
  function(i) format_row(well_data[i, ], as.character(well_data$replicate_index[i])),
  character(1)
)
pooled_rows <- vapply(
  seq_len(nrow(pooled_data)),
  function(i) format_row(pooled_data[i, ], "pooled"),
  character(1)
)

writeLines(
  c(well_rows, "%%SECTION:pooled%%", pooled_rows),
  file.path(tab_dir, "ddpcr_e200k_positive_well_results_rows.tex"),
  useBytes = TRUE
)

cat("Wrote figure and table outputs to manuscript/figures/ddpcr_e200k_positive_well_fa and manuscript/tables/ddpcr_e200k_positive_well_results\n")
