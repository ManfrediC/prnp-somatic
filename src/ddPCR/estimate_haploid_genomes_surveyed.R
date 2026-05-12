library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)

# -------------------------------------
# reproducible, repo-relative paths
# -------------------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
if (length(script_path) == 0) {
  stop("Could not determine script path. Please run with: Rscript src/ddPCR/estimate_haploid_genomes_surveyed.R")
}

script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

partition_counts_path <- file.path(project_root, "results", "ddPCR", "ddpcr_partition_counts_by_sample_assay.csv")
curated_ddpcr_path <- file.path(project_root, "results", "ddPCR", "SNV_data_final.xlsx")
output_dir <- file.path(project_root, "results", "ddPCR", "haploid_genomes_surveyed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# The partition-count CSV is produced by create_snv_dataframe.R after all
# existing ddPCR sample filtering and harmonisation. This script deliberately
# consumes that intermediate rather than touching the raw exports again.
if (!file.exists(partition_counts_path)) {
  stop("Missing partition-count input: ", partition_counts_path,
       "\nRun bash src/ddPCR/run_ddpcr.sh first.")
}

if (!file.exists(curated_ddpcr_path)) {
  stop("Missing curated ddPCR workbook: ", curated_ddpcr_path)
}

# -------------------------------------
# helpers
# -------------------------------------

# Poisson occupancy estimate from the fraction of droplets with no signal for
# the target class being counted. This returns genome-equivalent molecules, not
# observed positive droplets.
poisson_from_negative_fraction <- function(negative_fraction, n_total) {
  case_when(
    is.na(negative_fraction) | is.na(n_total) ~ NA_real_,
    n_total <= 0 ~ NA_real_,
    negative_fraction < 0 | negative_fraction > 1 ~ NA_real_,
    negative_fraction == 0 ~ Inf,
    TRUE ~ -n_total * log(negative_fraction)
  )
}

# Convenience wrapper for rows where the negative count is available directly.
poisson_from_negative_count <- function(n_negative, n_total) {
  poisson_from_negative_fraction(n_negative / n_total, n_total)
}

# Exact binomial interval for the negative-droplet fraction. The Poisson
# transform is monotone decreasing, so upper negative-fraction bounds become
# lower genome-count bounds, and vice versa.
exact_negative_fraction_ci <- function(n_negative, n_total, conf.level = 0.95) {
  alpha <- 1 - conf.level
  n_negative <- as.numeric(n_negative)
  n_total <- as.numeric(n_total)

  lower <- rep(NA_real_, length(n_negative))
  upper <- rep(NA_real_, length(n_negative))
  ok <- !is.na(n_negative) & !is.na(n_total) & n_total > 0

  lower[ok & n_negative == 0] <- 0
  upper[ok & n_negative == n_total] <- 1

  idx_lower <- ok & n_negative > 0
  lower[idx_lower] <- qbeta(
    alpha / 2,
    n_negative[idx_lower],
    n_total[idx_lower] - n_negative[idx_lower] + 1
  )

  idx_upper <- ok & n_negative < n_total
  upper[idx_upper] <- qbeta(
    1 - alpha / 2,
    n_negative[idx_upper] + 1,
    n_total[idx_upper] - n_negative[idx_upper]
  )

  tibble(lower = lower, upper = upper)
}

# Add REF, MUT, and total haploid genome-equivalent estimates to a table that
# already has accepted droplets and the three corresponding negative counts.
add_genome_estimates <- function(.data) {
  ref_ci <- exact_negative_fraction_ci(.data$n_ref_negative_droplets, .data$n_accepted_droplets)
  mut_ci <- exact_negative_fraction_ci(.data$n_mut_negative_droplets, .data$n_accepted_droplets)
  total_ci <- exact_negative_fraction_ci(.data$n_signal_negative_droplets, .data$n_accepted_droplets)

  .data %>%
    mutate(
      n_ref_genomes_poisson = poisson_from_negative_count(
        n_ref_negative_droplets,
        n_accepted_droplets
      ),
      n_ref_genomes_poisson_low = poisson_from_negative_fraction(
        # Upper negative-fraction CI bound maps to the lower genome estimate.
        ref_ci$upper,
        n_accepted_droplets
      ),
      n_ref_genomes_poisson_high = poisson_from_negative_fraction(
        # Lower negative-fraction CI bound maps to the upper genome estimate.
        ref_ci$lower,
        n_accepted_droplets
      ),
      n_mut_genomes_poisson = poisson_from_negative_count(
        n_mut_negative_droplets,
        n_accepted_droplets
      ),
      n_mut_genomes_poisson_low = poisson_from_negative_fraction(
        # Same reversed-bound logic for the MUT-negative fraction.
        mut_ci$upper,
        n_accepted_droplets
      ),
      n_mut_genomes_poisson_high = poisson_from_negative_fraction(
        mut_ci$lower,
        n_accepted_droplets
      ),
      n_haploid_genomes_poisson = poisson_from_negative_count(
        # Total genome equivalents are estimated from signal-negative droplets
        # directly, not by summing rounded REF and MUT estimates.
        n_signal_negative_droplets,
        n_accepted_droplets
      ),
      n_haploid_genomes_poisson_low = poisson_from_negative_fraction(
        total_ci$upper,
        n_accepted_droplets
      ),
      n_haploid_genomes_poisson_high = poisson_from_negative_fraction(
        total_ci$lower,
        n_accepted_droplets
      ),
      haploid_genomes_per_accepted_droplet =
        n_haploid_genomes_poisson / n_accepted_droplets
    )
}

# Sum droplets first, then estimate genomes from the aggregate negative
# fraction. Averaging or summing row-level estimates would give different
# confidence intervals and would be harder to audit against the displayed row.
summarise_counts <- function(.data, group_cols) {
  .data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_wells = sum(n_wells, na.rm = TRUE),
      n_accepted_droplets = sum(n_accepted_droplets, na.rm = TRUE),
      n_double_positive_droplets = sum(n_double_positive_droplets, na.rm = TRUE),
      n_ref_only_droplets = sum(n_ref_only_droplets, na.rm = TRUE),
      n_mut_only_droplets = sum(n_mut_only_droplets, na.rm = TRUE),
      n_signal_positive_droplets = sum(n_signal_positive_droplets, na.rm = TRUE),
      n_signal_negative_droplets = sum(n_signal_negative_droplets, na.rm = TRUE),
      n_ref_positive_droplets = sum(n_ref_positive_droplets, na.rm = TRUE),
      n_mut_positive_droplets = sum(n_mut_positive_droplets, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      n_ref_negative_droplets = n_accepted_droplets - n_ref_positive_droplets,
      n_mut_negative_droplets = n_accepted_droplets - n_mut_positive_droplets
    ) %>%
    add_genome_estimates()
}

# Formatting helpers are kept here so the CSV outputs remain numeric while the
# TeX table gets manuscript-friendly strings.
format_int <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_num <- function(x, digits = 3) {
  format(round(x, digits), nsmall = digits, scientific = FALSE, trim = TRUE)
}

format_ci <- function(estimate, low, high) {
  paste0(format_int(estimate), " (", format_int(low), "-", format_int(high), ")")
}

escape_latex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("_", "\\\\_", x)
  x
}

# -------------------------------------
# read and validate inputs
# -------------------------------------

partition_counts <- readr::read_csv(partition_counts_path, show_col_types = FALSE)

# Fail early if create_snv_dataframe.R stops emitting any field that the
# downstream supplementary calculation depends on.
required_partition_cols <- c(
  "participant", "group", "histotype", "code", "brain_region", "sample_id", "mutation",
  "n_wells", "n_accepted_droplets",
  "n_double_positive_droplets", "n_ref_only_droplets", "n_mut_only_droplets",
  "n_signal_positive_droplets", "n_signal_negative_droplets",
  "n_ref_positive_droplets", "n_mut_positive_droplets"
)

missing_partition_cols <- setdiff(required_partition_cols, names(partition_counts))
if (length(missing_partition_cols) > 0) {
  stop("Partition-count input is missing columns: ", paste(missing_partition_cols, collapse = ", "))
}

# These checks catch impossible droplet arithmetic before we join to
# SNV_data_final.xlsx. They are deliberately about physical count constraints,
# not manuscript interpretation.
impossible_counts <- partition_counts %>%
  filter(
    n_accepted_droplets < 0 |
      n_double_positive_droplets < 0 |
      n_ref_only_droplets < 0 |
      n_mut_only_droplets < 0 |
      n_signal_negative_droplets < 0 |
      n_signal_positive_droplets + n_signal_negative_droplets != n_accepted_droplets |
      n_ref_positive_droplets > n_accepted_droplets |
      n_mut_positive_droplets > n_accepted_droplets
  )

if (nrow(impossible_counts) > 0) {
  stop("Impossible partition counts detected:\n",
       paste(utils::capture.output(print(impossible_counts)), collapse = "\n"))
}

curated_ddpcr <- openxlsx::read.xlsx(curated_ddpcr_path) %>%
  as_tibble()

# SNV_data_final.xlsx is validation-only here. The supplementary script should
# not change its schema or use it as a second source of droplet truth.
required_curated_cols <- c(
  "participant", "group", "histotype", "code", "brain_region", "mutation",
  "n_mut_droplets", "n_total_droplets"
)

missing_curated_cols <- setdiff(required_curated_cols, names(curated_ddpcr))
if (length(missing_curated_cols) > 0) {
  stop("Curated ddPCR workbook is missing columns: ", paste(missing_curated_cols, collapse = ", "))
}

key_cols <- c("code", "brain_region", "mutation")

# Duplicate keys would make the full join ambiguous and could hide accidental
# row multiplication, so stop before any totals are calculated.
partition_duplicate_keys <- partition_counts %>%
  count(across(all_of(key_cols)), name = "n") %>%
  filter(n > 1)

curated_duplicate_keys <- curated_ddpcr %>%
  count(across(all_of(key_cols)), name = "n") %>%
  filter(n > 1)

if (nrow(partition_duplicate_keys) > 0 || nrow(curated_duplicate_keys) > 0) {
  stop("Duplicate validation keys detected in partition counts or SNV_data_final.")
}

# Validate the new denominator intermediate against the existing curated table
# at the exact same grain. This protects the older ddPCR outputs from quiet
# behavioural drift.
validation <- partition_counts %>%
  select(all_of(key_cols), participant, group, histotype,
         n_accepted_droplets, n_mut_positive_droplets) %>%
  full_join(
    curated_ddpcr %>%
      select(all_of(key_cols),
             participant_curated = participant,
             group_curated = group,
             histotype_curated = histotype,
             n_total_droplets,
             n_mut_droplets),
    by = key_cols
  )

validation_errors <- validation %>%
  filter(
    is.na(participant) |
      is.na(participant_curated) |
      participant != participant_curated |
      group != group_curated |
      histotype != histotype_curated |
      n_accepted_droplets != n_total_droplets |
      n_mut_positive_droplets != n_mut_droplets
  )

if (nrow(validation_errors) > 0) {
  stop("Partition counts do not validate against SNV_data_final:\n",
       paste(utils::capture.output(print(validation_errors)), collapse = "\n"))
}

# -------------------------------------
# sample-region and participant-level outputs
# -------------------------------------

# Sample-region output is the most detailed supplementary table: one row per
# existing biological sample-region by mutation assay. The known heterozygous
# E200K sample is retained because denominator depth, not somatic status, is
# the target of this analysis.
sample_region <- partition_counts %>%
  mutate(
    known_heterozygous_e200k = code == "14-2" & mutation == "E200K",
    n_ref_negative_droplets = n_accepted_droplets - n_ref_positive_droplets,
    n_mut_negative_droplets = n_accepted_droplets - n_mut_positive_droplets
  ) %>%
  add_genome_estimates()

# Participant-pooled output collapses brain regions within each assay while
# preserving the same Poisson-from-negative-fraction logic.
participant_pooled <- summarise_counts(
  sample_region,
  c("participant", "group", "histotype", "code", "mutation", "known_heterozygous_e200k")
)

readr::write_csv(sample_region, file.path(output_dir, "sample_region_haploid_genomes.csv"))
readr::write_csv(participant_pooled, file.path(output_dir, "participant_pooled_haploid_genomes.csv"))

# -------------------------------------
# cohort-level TeX supplement table
# -------------------------------------

# The TeX table is the manuscript-facing summary. The CSVs above keep the
# detailed rows, while this table reports CJD/control/all rows by mutation plus
# clearly labelled all-assay rows.
cohort_by_group_mutation <- sample_region %>%
  group_by(group, mutation) %>%
  summarise(
    n_participants = n_distinct(participant),
    n_sample_regions = n_distinct(paste(code, brain_region, sep = "__")),
    n_wells = sum(n_wells, na.rm = TRUE),
    n_accepted_droplets = sum(n_accepted_droplets, na.rm = TRUE),
    n_double_positive_droplets = sum(n_double_positive_droplets, na.rm = TRUE),
    n_ref_only_droplets = sum(n_ref_only_droplets, na.rm = TRUE),
    n_mut_only_droplets = sum(n_mut_only_droplets, na.rm = TRUE),
    n_signal_positive_droplets = sum(n_signal_positive_droplets, na.rm = TRUE),
    n_signal_negative_droplets = sum(n_signal_negative_droplets, na.rm = TRUE),
    n_ref_positive_droplets = sum(n_ref_positive_droplets, na.rm = TRUE),
    n_mut_positive_droplets = sum(n_mut_positive_droplets, na.rm = TRUE),
    .groups = "drop"
  )

# All-assay rows are assay-level genome-equivalent observations. They should
# not be read as unique genomes, because the same DNA sample may be assayed for
# multiple mutations.
cohort_by_group_all <- sample_region %>%
  mutate(mutation = "all_assays") %>%
  group_by(group, mutation) %>%
  summarise(
    n_participants = n_distinct(participant),
    n_sample_regions = n_distinct(paste(code, brain_region, sep = "__")),
    n_wells = sum(n_wells, na.rm = TRUE),
    n_accepted_droplets = sum(n_accepted_droplets, na.rm = TRUE),
    n_double_positive_droplets = sum(n_double_positive_droplets, na.rm = TRUE),
    n_ref_only_droplets = sum(n_ref_only_droplets, na.rm = TRUE),
    n_mut_only_droplets = sum(n_mut_only_droplets, na.rm = TRUE),
    n_signal_positive_droplets = sum(n_signal_positive_droplets, na.rm = TRUE),
    n_signal_negative_droplets = sum(n_signal_negative_droplets, na.rm = TRUE),
    n_ref_positive_droplets = sum(n_ref_positive_droplets, na.rm = TRUE),
    n_mut_positive_droplets = sum(n_mut_positive_droplets, na.rm = TRUE),
    .groups = "drop"
  )

# Add mutation-specific totals across CJD and control groups for readers who
# want the overall survey depth per assay.
cohort_all_groups_mutation <- sample_region %>%
  mutate(group = "all") %>%
  group_by(group, mutation) %>%
  summarise(
    n_participants = n_distinct(participant),
    n_sample_regions = n_distinct(paste(code, brain_region, sep = "__")),
    n_wells = sum(n_wells, na.rm = TRUE),
    n_accepted_droplets = sum(n_accepted_droplets, na.rm = TRUE),
    n_double_positive_droplets = sum(n_double_positive_droplets, na.rm = TRUE),
    n_ref_only_droplets = sum(n_ref_only_droplets, na.rm = TRUE),
    n_mut_only_droplets = sum(n_mut_only_droplets, na.rm = TRUE),
    n_signal_positive_droplets = sum(n_signal_positive_droplets, na.rm = TRUE),
    n_signal_negative_droplets = sum(n_signal_negative_droplets, na.rm = TRUE),
    n_ref_positive_droplets = sum(n_ref_positive_droplets, na.rm = TRUE),
    n_mut_positive_droplets = sum(n_mut_positive_droplets, na.rm = TRUE),
    .groups = "drop"
  )

# Final grand-total row, again at assay-observation level rather than unique
# biological genomes.
cohort_all_groups_all <- sample_region %>%
  mutate(group = "all", mutation = "all_assays") %>%
  group_by(group, mutation) %>%
  summarise(
    n_participants = n_distinct(participant),
    n_sample_regions = n_distinct(paste(code, brain_region, sep = "__")),
    n_wells = sum(n_wells, na.rm = TRUE),
    n_accepted_droplets = sum(n_accepted_droplets, na.rm = TRUE),
    n_double_positive_droplets = sum(n_double_positive_droplets, na.rm = TRUE),
    n_ref_only_droplets = sum(n_ref_only_droplets, na.rm = TRUE),
    n_mut_only_droplets = sum(n_mut_only_droplets, na.rm = TRUE),
    n_signal_positive_droplets = sum(n_signal_positive_droplets, na.rm = TRUE),
    n_signal_negative_droplets = sum(n_signal_negative_droplets, na.rm = TRUE),
    n_ref_positive_droplets = sum(n_ref_positive_droplets, na.rm = TRUE),
    n_mut_positive_droplets = sum(n_mut_positive_droplets, na.rm = TRUE),
    .groups = "drop"
  )

# Bind the four summary blocks, derive negative counts from their own aggregate
# droplet totals, then calculate estimates and CIs from those same rows.
cohort_table <- bind_rows(
  cohort_by_group_mutation,
  cohort_by_group_all,
  cohort_all_groups_mutation,
  cohort_all_groups_all
) %>%
  mutate(
    n_ref_negative_droplets = n_accepted_droplets - n_ref_positive_droplets,
    n_mut_negative_droplets = n_accepted_droplets - n_mut_positive_droplets
  ) %>%
  add_genome_estimates() %>%
  mutate(
    group_label = recode(group, prion = "CJD", control = "Control", all = "All"),
    mutation_label = recode(mutation, all_assays = "All assays"),
    group_label = factor(group_label, levels = c("CJD", "Control", "All")),
    mutation_label = factor(mutation_label, levels = c("D178N", "E200K", "P102L", "All assays"))
  ) %>%
  arrange(group_label, mutation_label)

# Keep calculated numeric columns in the CSVs; convert to formatted strings
# only for the TeX table.
tex_table <- cohort_table %>%
  transmute(
    Group = escape_latex(as.character(group_label)),
    Mutation = escape_latex(as.character(mutation_label)),
    Participants = format_int(n_participants),
    `Sample-regions` = format_int(n_sample_regions),
    Wells = format_int(n_wells),
    `Accepted droplets` = format_int(n_accepted_droplets),
    `REF+ droplets` = format_int(n_ref_positive_droplets),
    `MUT+ droplets` = format_int(n_mut_positive_droplets),
    `Signal+ droplets` = format_int(n_signal_positive_droplets),
    `Signal- droplets` = format_int(n_signal_negative_droplets),
    `REF haploid genomes (95% CI)` = format_ci(
      n_ref_genomes_poisson,
      n_ref_genomes_poisson_low,
      n_ref_genomes_poisson_high
    ),
    `MUT haploid genomes (95% CI)` = format_ci(
      n_mut_genomes_poisson,
      n_mut_genomes_poisson_low,
      n_mut_genomes_poisson_high
    ),
    `Total haploid genomes (95% CI)` = format_ci(
      n_haploid_genomes_poisson,
      n_haploid_genomes_poisson_low,
      n_haploid_genomes_poisson_high
    ),
    `Genomes per accepted droplet` = format_num(haploid_genomes_per_accepted_droplet, 3)
  )

# Build a small booktabs table without extra dependencies.
tex_header <- paste(escape_latex(names(tex_table)), collapse = " & ")
tex_body <- apply(tex_table, 1, function(row) paste0(paste(row, collapse = " & "), " \\\\"))

tex_lines <- c(
  "% Auto-generated by src/ddPCR/estimate_haploid_genomes_surveyed.R",
  "\\begin{table}[t]",
  "\\centering",
  "\\scriptsize",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{llrrrrrrrrlllr}",
  "\\toprule",
  paste0(tex_header, " \\\\"),
  "\\midrule",
  tex_body,
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\caption{\\textbf{ddPCR droplet counts and estimated haploid genome equivalents surveyed.} Counts are shown after applying the existing ddPCR sample inclusion and quality-control workflow. REF+ and MUT+ droplets include double-positive droplets. Haploid genome-equivalent estimates and 95\\% confidence intervals were derived from the corresponding negative-droplet fractions using a Poisson occupancy model. All-assay rows are assay-level genome-equivalent observations, not unique biological genomes, because the same DNA sample can contribute to multiple mutation assays.}",
  "\\label{tab:ddpcr_haploid_genomes}",
  "\\end{table}",
  ""
)

writeLines(tex_lines, file.path(output_dir, "ddpcr_haploid_genomes_supplement_table.tex"))

# -------------------------------------
# run summary
# -------------------------------------

# The summary is a quick audit file: it records which inputs were validated,
# the formula used, and the headline grand total without requiring TeX parsing.
grand_total <- cohort_table %>%
  filter(group_label == "All", mutation_label == "All assays") %>%
  slice(1)

summary_lines <- c(
  "ddPCR droplet counts and haploid genome-equivalent estimates",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Project root: ", project_root),
  paste0("Partition-count rows read: ", nrow(partition_counts)),
  paste0("Curated ddPCR rows validated: ", nrow(curated_ddpcr)),
  paste0("Participants retained: ", n_distinct(partition_counts$participant)),
  paste0("Known heterozygous E200K sample-region rows retained: ",
         sum(sample_region$known_heterozygous_e200k)),
  "",
  "Formula:",
  "  n_genomes = -n_accepted_droplets * log(n_negative_droplets / n_accepted_droplets)",
  "",
  "Grand total across all biological samples and assays:",
  paste0("  participants = ", grand_total$n_participants),
  paste0("  sample-regions = ", grand_total$n_sample_regions),
  paste0("  wells = ", grand_total$n_wells),
  paste0("  accepted droplets = ", grand_total$n_accepted_droplets),
  paste0("  signal-positive droplets = ", grand_total$n_signal_positive_droplets),
  paste0("  signal-negative droplets = ", grand_total$n_signal_negative_droplets),
  paste0("  total haploid genomes = ", round(grand_total$n_haploid_genomes_poisson, 3)),
  paste0("  total haploid genomes 95% CI = ",
         round(grand_total$n_haploid_genomes_poisson_low, 3),
         "-",
         round(grand_total$n_haploid_genomes_poisson_high, 3))
)

writeLines(summary_lines, file.path(output_dir, "run_summary.txt"))
