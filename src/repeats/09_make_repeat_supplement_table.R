#!/usr/bin/env Rscript

# This script builds a compact one-row-per-sample supplement table from the
# PRNP ORR repeat workflow outputs. The table includes both controls and CJD
# samples and writes a CSV plus a LaTeX row output that matches manuscript table
# assembly style.

# This block loads the tidyverse packages used for reading, joining, formatting,
# and writing the supplement table.
suppressPackageStartupMessages(library(tidyverse))

# This helper finds the directory that contains the running script so the input
# and output paths work regardless of the current working directory.
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]),
                                 winslash = "/",
                                 mustWork = TRUE)))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile,
                                 winslash = "/",
                                 mustWork = TRUE)))
  }

  return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
}

# This helper reads a tab-delimited workflow file with quiet column parsing.
read_tsv <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE)
}

# This helper converts internal group labels into cleaner display labels.
format_group <- function(x) {
  dplyr::case_when(
    x == "control" ~ "Control",
    x == "cjd" ~ "CJD",
    TRUE ~ x
  )
}

# This helper expands control sample IDs for readability in the manuscript table.
format_sample <- function(x) {
  dplyr::if_else(stringr::str_starts(x, "Ctrl"), str_replace(x, "^Ctrl", "Control"), x)
}

# This helper converts internal call labels into cleaner display labels.
format_call <- function(x) {
  dplyr::case_when(
    x == "reference" ~ "Reference",
    x == "candidate_OPRI" ~ "Candidate OPRI",
    x == "candidate_OPRD" ~ "Candidate OPRD",
    x == "uncertain" ~ "Uncertain",
    x == "failed_qc" ~ "Failed QC",
    x == "disabled" ~ "Disabled",
    TRUE ~ x
  )
}

# This helper formats the CJD-only background-aware filter for a shared table.
format_background_filter <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    x == "yes" ~ "Pass",
    x == "no" ~ "No candidate",
    TRUE ~ x
  )
}

# This helper formats numeric coverage values with one decimal place.
format_one_decimal <- function(x) {
  dplyr::if_else(is.na(x), NA_character_, sprintf("%.1f", as.numeric(x)))
}

# This helper extracts the numeric suffix from sample IDs for natural sorting.
sample_number <- function(x) {
  suppressWarnings(as.integer(sub("^[A-Za-z]+", "", x)))
}

# Define the input and output paths relative to the repository root.
script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", ".."),
                             winslash = "/",
                             mustWork = TRUE)
results_dir <- file.path(repo_root, "results", "repeats")
output_dir <- file.path(results_dir, "summary")
output_csv <- file.path(output_dir, "repeat_supplement_table.csv")
output_pub_tsv <- file.path(output_dir, "repeat_supplement_table_publication.tsv")
output_tex_rows <- file.path(output_dir, "repeat_supplement_table_rows.tex")

manifest_path <- file.path(results_dir, "sample_manifest.tsv")
sample_calls_path <- file.path(results_dir, "sample_calls.tsv")
gangstr_path <- file.path(results_dir, "gangstr_calls.tsv")
manual_controls_path <- file.path(results_dir,
                                  "manual_cohort",
                                  "controls",
                                  "cohort_summary.tsv")
manual_cjd_path <- file.path(results_dir,
                             "manual_cohort",
                             "cjd",
                             "cohort_summary.tsv")
filtered_cjd_path <- file.path(results_dir,
                               "manual_cohort",
                               "cjd",
                               "filtered",
                               "default.annotated.tsv")

# Read the workflow tables that feed the supplement summary.
manifest <- read_tsv(manifest_path)
sample_calls <- read_tsv(sample_calls_path)
gangstr_calls <- read_tsv(gangstr_path)
manual_controls <- read_tsv(manual_controls_path)
manual_cjd <- read_tsv(manual_cjd_path)
filtered_cjd <- read_tsv(filtered_cjd_path)

# Keep only the columns needed for the compact supplement table.
eh_selected <- sample_calls %>%
  dplyr::select(sample_id,
                interpretation,
                LC,
                ADSP,
                ADFL,
                ADIR)

gangstr_selected <- gangstr_calls %>%
  dplyr::select(sample_id,
                gangstr_interpretation,
                DP,
                enclosing_nonreference_reads)

manual_selected <- dplyr::bind_rows(
  manual_controls %>%
    dplyr::select(sample_id,
                  exact_nonreference_reads,
                  synthetic_high_or_medium_nonreference_reads,
                  one_sided_indel_or_softclip_reads),
  manual_cjd %>%
    dplyr::select(sample_id,
                  exact_nonreference_reads,
                  synthetic_high_or_medium_nonreference_reads,
                  one_sided_indel_or_softclip_reads)
)

filter_selected <- filtered_cjd %>%
  dplyr::select(sample_id,
                filter_pass)

# Join all tool outputs into a single one-row-per-sample table.
supplement_table <- manifest %>%
  dplyr::select(sample_id, group, orr_overlapping_reads) %>%
  dplyr::left_join(eh_selected, by = "sample_id") %>%
  dplyr::left_join(gangstr_selected, by = "sample_id") %>%
  dplyr::left_join(manual_selected, by = "sample_id") %>%
  dplyr::left_join(filter_selected, by = "sample_id")

# This block formats values, applies natural sample ordering, and sets clean
# supplement-friendly column headers.
supplement_table <- supplement_table %>%
  dplyr::mutate(
    sample_id = format_sample(sample_id),
    group = format_group(group),
    interpretation = format_call(interpretation),
    gangstr_interpretation = format_call(gangstr_interpretation),
    filter_pass = format_background_filter(filter_pass),
    LC = format_one_decimal(LC),
    sample_order_group = dplyr::case_when(
      group == "CJD" ~ 1L,
      group == "Control" ~ 2L,
      TRUE ~ 3L
    ),
    sample_order_number = sample_number(sample_id)
  ) %>%
  dplyr::arrange(sample_order_group, sample_order_number, sample_id) %>%
  dplyr::select(
    -sample_order_group,
    -sample_order_number,
    -group
  ) %>%
  dplyr::rename(
    "Sample" = sample_id,
    "ORR overlapping reads" = orr_overlapping_reads,
    "ExpansionHunter call" = interpretation,
    "LC" = LC,
    "ADSP" = ADSP,
    "ADFL" = ADFL,
    "ADIR" = ADIR,
    "GangSTR call" = gangstr_interpretation,
    "GangSTR depth (DP)" = DP,
    "GangSTR enclosing nonreference reads" = enclosing_nonreference_reads,
    "Manual exact nonreference reads" = exact_nonreference_reads,
    "Manual synthetic nonreference reads" = synthetic_high_or_medium_nonreference_reads,
    "One-sided indel/soft-clip reads" = one_sided_indel_or_softclip_reads,
    "Control-aware filter" = filter_pass
  )

# Ensure the summary output directory exists.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Write the supplement table as a CSV for manuscript submission tables.
readr::write_csv(supplement_table,
                 file = output_csv,
                 na = "")

publication_table <- supplement_table %>%
  select(
    Sample,
    `ORR overlapping reads`,
    `ExpansionHunter call`,
    LC,
    ADSP,
    ADFL,
    ADIR,
    `GangSTR call`,
    `GangSTR depth (DP)`,
    `GangSTR enclosing nonreference reads`,
    `Manual exact nonreference reads`,
    `Manual synthetic nonreference reads`,
    `One-sided indel/soft-clip reads`,
    `Control-aware filter`
  )
write.table(
  publication_table,
  file = output_pub_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

# create LaTeX row output
latex_safe <- function(x) {
  x <- ifelse(is.na(x), "NA", as.character(x))
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#\\$])", "\\\\\\1", x)
  x
}
latex_rows <- apply(publication_table, 1, function(row_vals) {
  paste0(paste(latex_safe(row_vals), collapse = " & "), " \\\\")
})
writeLines(latex_rows, con = output_tex_rows, useBytes = TRUE)

# Print a short completion message for command-line use.
message("Wrote supplement table: ", output_csv)
message("Wrote publication TSV: ", output_pub_tsv)
message("Wrote LaTeX rows: ", output_tex_rows)
