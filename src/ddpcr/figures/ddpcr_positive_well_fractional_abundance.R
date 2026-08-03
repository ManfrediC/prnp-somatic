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
#   manuscript/tables/supplement/table_ddpcr_e200k_positive_well_results.tex

library(tidyverse)
library(readxl)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
fig_dir <- file.path(project_root, "manuscript", "figures", "ddpcr_e200k_positive_well_fa")
tab_dir <- file.path(project_root, "manuscript", "tables", "ddpcr_e200k_positive_well_results")
supplement_tab_dir <- file.path(project_root, "manuscript", "tables", "supplement")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplement_tab_dir, recursive = TRUE, showWarnings = FALSE)

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
  filter(
    mutation == "E200K",
    paste(participant_display, brain_region) %in% c("CJD4 ps", "CJD21 th")
  ) %>%
  rename(pooled_fa = fractional_abundance, pooled_ci_low = ci_low, pooled_ci_high = ci_high) %>%
  transmute(
    participant = participant_display,
    brain_region,
    sample = paste(participant_display, brain_region_display),
    replicate = as.character(replicate_index),
    run_id, run_date = as.Date(run_date), well,
    fa = individual_fractional_abundance,
    ci_low = individual_ci_low,
    ci_high = individual_ci_high,
    lob_count = individual_lob_count,
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

wells <- wells %>%
  left_join(
    well_counts,
    by = c("run_id", "well"),
    relationship = "many-to-one"
  )

# ---- pooled sample-level data ------------------------------------------------

pooled <- read_excel(file.path(project_root, "results", "ddPCR", "SNV_data_final.xlsx")) %>%
  filter(
    mutation == "E200K",
    is_pooled,
    paste(participant, brain_region) %in% c("CJD4 ps", "CJD21 th")
  ) %>%
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
    ci_low, ci_high, lob_count, lob_fa,
    above_lob = detected_above_LoB,
    above_lod = detected_above_LoD
  )

# ---- cross-checks ------------------------------------------------------------

expected_preparation_counts <- c("CJD4 pons" = 3L, "CJD21 thalamus" = 2L)
preparation_counts <- table(wells$sample)

stopifnot(
  nrow(wells) == 5L,
  nrow(pooled) == 2L,
  setequal(names(preparation_counts), names(expected_preparation_counts)),
  setequal(pooled$sample, names(expected_preparation_counts)),
  all(as.integer(preparation_counts[names(expected_preparation_counts)]) == expected_preparation_counts),
  !any(c(wells$participant, pooled$participant) == "CJD30"),
  n_distinct(wells$run_id, wells$well) == nrow(wells),
  !anyNA(wells %>% select(run_id, run_date, well, accepted, positives, fa, ci_low, ci_high, lob_count, lob_fa, above_lob, above_lod)),
  !anyNA(pooled %>% select(accepted, positives, fa, ci_low, ci_high, lob_count, lob_fa, above_lob, above_lod)),
  all(near(wells$lob_fa, 100 * wells$lob_count / wells$accepted)),
  all(wells$above_lob == (wells$positives > wells$lob_count)),
  all(wells$above_lod == (wells$ci_low > lod_e200k)),
  all(near(pooled$lob_fa, 100 * pooled$lob_count / pooled$accepted)),
  all(pooled$above_lob == (pooled$positives > pooled$lob_count)),
  all(pooled$above_lod == (pooled$ci_low > lod_e200k))
)

sums_ok <- wells %>%
  group_by(sample) %>%
  summarise(accepted = sum(accepted), positives = sum(positives), .groups = "drop") %>%
  inner_join(pooled, by = c("sample", "accepted", "positives")) %>%
  nrow() == 2L

pooled_ok <- wells %>%
  distinct(sample, pooled_fa, pooled_ci_low, pooled_ci_high) %>%
  inner_join(pooled, by = "sample") %>%
  summarise(ok = all(near(pooled_fa, fa) & near(pooled_ci_low, ci_low) & near(pooled_ci_high, ci_high))) %>%
  pull(ok)

stopifnot(sums_ok, pooled_ok)

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

# ---- table data: audit CSV + complete TeX ------------------------------------

# Machine-readable table: all five independent preparations followed by the
# two pooled sample-region summaries. Run date, well ID and LoB count are kept
# in the CSV for direct auditing but are omitted from the publication table.
table_data <- bind_rows(
  wells %>%
    mutate(
      measurement = "preparation",
      preparation_index = as.integer(replicate)
    ),
  pooled %>%
    mutate(
      measurement = "pooled",
      preparation_index = NA_integer_
    )
) %>%
  mutate(
    measurement = factor(measurement, levels = c("preparation", "pooled")),
    sample = factor(sample, levels = names(expected_preparation_counts))
  ) %>%
  arrange(measurement, sample, preparation_index) %>%
  transmute(
    participant,
    brain_region = recode(brain_region, !!!region_labels),
    measurement = as.character(measurement),
    preparation_index,
    run_date,
    well,
    accepted_droplets = as.integer(accepted),
    mut_positive_droplets = as.integer(positives),
    fa_pct = round(fa, 3),
    ci_low_pct = round(ci_low, 3),
    ci_high_pct = round(ci_high, 3),
    lob_count = as.integer(lob_count),
    lob_fa_pct = round(lob_fa, 3),
    above_lob,
    above_lod
  )

stopifnot(
  nrow(table_data) == 7L,
  sum(table_data$measurement == "preparation") == 5L,
  sum(table_data$measurement == "pooled") == 2L
)

write_csv(table_data, file.path(tab_dir, "ddpcr_e200k_positive_well_results.csv"))

format_count <- function(value) {
  format(as.integer(value), big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_tex_row <- function(row, preparation_label) {
  paste0(
    "\\rule{0pt}{1.2em}",
    paste(
      row$participant,
      row$brain_region,
      preparation_label,
      format_count(row$accepted_droplets),
      format_count(row$mut_positive_droplets),
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

preparation_data <- table_data %>% filter(measurement == "preparation")
pooled_data <- table_data %>% filter(measurement == "pooled")

preparation_rows <- vapply(
  seq_len(nrow(preparation_data)),
  function(i) {
    format_tex_row(
      preparation_data[i, ],
      as.character(preparation_data$preparation_index[i])
    )
  },
  character(1)
)
pooled_rows <- vapply(
  seq_len(nrow(pooled_data)),
  function(i) format_tex_row(pooled_data[i, ], "Pooled"),
  character(1)
)

column_spec <- "cccccccccc"

caption_text <- paste0(
  "\\textbf{Detailed ddPCR results for the E200K-positive samples CJD4 pons ",
  "and CJD21 thalamus.} Preparation rows report individual ddPCR wells from ",
  "independent DNA preparations; pooled rows combine all preparations from ",
  "one sample-region. FA: fractional abundance of the E200K allele; 95\\% CI: ",
  "Poisson-based confidence interval; LoB FA: measurement-specific limit of ",
  "blank on the FA scale. The LoB criterion is met when the mutant-positive ",
  "droplet count exceeds the 95th percentile of expected blank noise; the LoD ",
  "criterion is met when the lower 95\\% CI bound exceeds the E200K limit of ",
  "detection (0.067\\%). CJD30, the germline E200K carrier, is not included. ",
  "Abbreviations: FA, fractional abundance; CI, confidence interval; LoB, ",
  "limit of blank; LoD, limit of detection."
)

table_tex <- c(
  "% Auto-generated by src/ddpcr/figures/ddpcr_positive_well_fractional_abundance.R",
  "\\begin{table}[htbp]",
  "\t\\ra{1.0}",
  "\t\\centering",
  "\t{\\normalsize",
  "\t\t\\begin{adjustbox}{max width=\\textwidth}",
  paste0("\t\t\t\\begin{tabular}{", column_spec, "}"),
  "\t\t\t\t\\toprule[2pt]",
  paste0(
    "\t\t\t\t\\multirow{2}{*}[-0.6ex]{\\textbf{Sample}} & ",
    "\\multirow{2}{*}[-0.6ex]{\\makecell{\\textbf{Brain}\\\\\\textbf{region}}} & ",
    "\\multirow{2}{*}[-0.6ex]{\\textbf{Preparation}} & ",
    "\\multicolumn{2}{c}{\\textbf{Droplets}} & ",
    "\\multicolumn{3}{c}{\\textbf{Fractional abundance (\\%)}} & ",
    "\\multicolumn{2}{c}{\\textbf{Criteria passed}} \\\\"
  ),
  "\t\t\t\t\\cmidrule(lr){4-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-10}",
  paste0(
    "\t\t\t\t& & & \\textbf{Accepted} & ",
    "\\makecell{\\textbf{E200K}\\\\\\textbf{positive}} & ",
    "\\textbf{Estimate} & \\textbf{95\\% CI} & \\textbf{LoB} & ",
    "\\textbf{LoB} & \\textbf{LoD} \\\\"
  ),
  "\t\t\t\t\\midrule",
  paste0("\t\t\t\t", preparation_rows),
  "\t\t\t\t\\addlinespace[2pt]",
  paste0("\t\t\t\t", pooled_rows),
  "\t\t\t\t\\bottomrule[2pt]",
  "\t\t\t\\end{tabular}",
  "\t\t\\end{adjustbox}",
  "\t}",
  "",
  "\t\\vspace{0.5em}",
  "",
  "\t\\begin{adjustbox}{max width=\\textwidth}",
  "\t\t\\begin{minipage}{\\linewidth}",
  paste0("\t\t\t\\caption{", caption_text, "}"),
  "\t\t\t\\label{tab:ddpcr_e200k_positive_well_results}",
  "\t\t\\end{minipage}",
  "\t\\end{adjustbox}",
  "",
  "\\end{table}",
  "",
  "\\clearpage"
)

writeLines(
  table_tex,
  file.path(supplement_tab_dir, "table_ddpcr_e200k_positive_well_results.tex"),
  useBytes = TRUE
)

cat(
  "Wrote E200K positive-well figure, audit CSV and complete supplementary Table S5 TeX\n"
)
