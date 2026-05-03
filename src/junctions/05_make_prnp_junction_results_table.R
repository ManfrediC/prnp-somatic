#!/usr/bin/env Rscript

# 05_make_prnp_junction_results_table.R
# Build a tidy per-sample PRNP junction QC table as a TSV.

library(Rsamtools)
library(dplyr)

## ---- resolve script and repository paths ----
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

## ---- configure default inputs, outputs and analysis thresholds ----
default_junc_dir <- file.path(repo_root, "results", "junctions", "junction_align")
default_out_dir <- file.path(repo_root, "results", "junctions", "junction_counts")
default_out_file <- file.path(default_out_dir, "prnp_junction_qc_table.tsv")
default_counts_file <- file.path(default_out_dir, "prnp_junction_counts.tsv")

JUNC_DIR <- Sys.getenv("PRNP_JUNCTION_ALIGN_DIR", unset = default_junc_dir)
OUT_FILE <- Sys.getenv("PRNP_JUNCTION_QC_TABLE", unset = default_out_file)
COUNTS_FILE <- Sys.getenv("PRNP_JUNCTION_COUNTS_FILE", unset = default_counts_file)

MIN_MAPQ <- 20L

if (!dir.exists(JUNC_DIR)) stop("Junction alignment directory not found: ", JUNC_DIR)
if (!file.exists(COUNTS_FILE)) stop("Required count table not found: ", COUNTS_FILE)
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)

## ---- load script-04 output used for final retained-read counts ----
counts_tbl <- read.delim(COUNTS_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample", "rname", "n_reads") %in% names(counts_tbl))) {
  stop("Count table is missing required columns: sample, rname, n_reads")
}

junction_rnames <- unique(counts_tbl$rname)
retained_tbl <- counts_tbl %>%
  group_by(sample) %>%
  summarise(junction_overhang10_reads = sum(n_reads, na.rm = TRUE), .groups = "drop")

## ---- collect all PRNP junction BAMs (one per sample) ----
bam_files <- list.files(JUNC_DIR, pattern = "\\.PRNP\\.toJunc\\.bam$", full.names = TRUE)
if (!length(bam_files)) stop("No *.PRNP.toJunc.bam found in ", JUNC_DIR)

## ---- helper: sort rows as CJD1, CJD2, ... then Ctrl1, Ctrl2, ... ----
sample_sort_index <- function(sample_names) {
  group <- ifelse(grepl("^CJD\\d+$", sample_names), "CJD",
                  ifelse(grepl("^Ctrl\\d+$", sample_names), "Ctrl", "Other"))
  num <- suppressWarnings(as.integer(sub("^[A-Za-z]+", "", sample_names)))
  num[is.na(num)] <- 999999L
  order(factor(group, levels = c("CJD", "Ctrl", "Other")), num, sample_names)
}

## ---- helper: format numeric values with two decimal places ----
fmt_two_dp <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.2f", round(x, 2)))
}

## ---- helper: format internal sample IDs for manuscript display ----
sample_display_name <- function(sample_names) {
  sub("^Ctrl", "Control", sample_names)
}

## ---- per-sample QC counting: flags, MAPQ and junction-overhang support ----
qc_rows <- lapply(bam_files, function(bam_junc) {
  sample <- sub("\\.PRNP\\.toJunc\\.bam$", "", basename(bam_junc))

  x <- scanBam(
    bam_junc,
    param = ScanBamParam(what = c("flag", "rname", "mapq"))
  )[[1]]

  if (length(x$flag) == 0L) {
    return(data.frame(
      sample = sample,
      unmapped_aln = 0L,
      secondary_aln = 0L,
      supplementary_aln = 0L,
      duplicate_aln = 0L,
      primary_nondup_mapped_aln = 0L,
      mapq_ge20_aln = 0L,
      junction_overhang10_reads = 0L,
      stringsAsFactors = FALSE
    ))
  }

  flag <- x$flag
  rname <- as.character(x$rname)
  mapq <- x$mapq

  unmapped_aln <- sum(bitwAnd(flag, 0x4) != 0L)
  secondary_aln <- sum(bitwAnd(flag, 0x100) != 0L)
  supplementary_aln <- sum(bitwAnd(flag, 0x800) != 0L)
  duplicate_aln <- sum(bitwAnd(flag, 0x400) != 0L)

  primary_nondup_mapped <- (bitwAnd(flag, 0x4) == 0L) &
    (bitwAnd(flag, 0x100) == 0L) &
    (bitwAnd(flag, 0x800) == 0L) &
    (bitwAnd(flag, 0x400) == 0L) &
    (rname %in% junction_rnames)

  high_mapq_keep <- primary_nondup_mapped & (mapq >= MIN_MAPQ)
  primary_nondup_mapped_aln <- sum(primary_nondup_mapped)
  mapq_ge20_aln <- sum(high_mapq_keep)

  data.frame(
    sample = sample,
    unmapped_aln = as.integer(unmapped_aln),
    secondary_aln = as.integer(secondary_aln),
    supplementary_aln = as.integer(supplementary_aln),
    duplicate_aln = as.integer(duplicate_aln),
    primary_nondup_mapped_aln = as.integer(primary_nondup_mapped_aln),
    mapq_ge20_aln = as.integer(mapq_ge20_aln),
    stringsAsFactors = FALSE
  )
})

## ---- combine rows, apply requested sort/format, and enforce column order ----
qc_table <- bind_rows(qc_rows)
qc_table <- qc_table %>% left_join(retained_tbl, by = "sample")
qc_table$junction_overhang10_reads <- dplyr::coalesce(qc_table$junction_overhang10_reads, 0)
qc_table <- qc_table[sample_sort_index(qc_table$sample), , drop = FALSE]

qc_table <- qc_table %>%
  select(
    sample,
    unmapped_aln,
    secondary_aln,
    supplementary_aln,
    duplicate_aln,
    primary_nondup_mapped_aln,
    mapq_ge20_aln,
    junction_overhang10_reads
  )

numeric_cols <- setdiff(names(qc_table), "sample")
qc_table[numeric_cols] <- lapply(qc_table[numeric_cols], fmt_two_dp)

## ---- write final tidy table as TSV ----
write.table(
  qc_table,
  file = OUT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("Wrote:\n  ", OUT_FILE, "\n")

# ---- create publication-ready table ----

# rename colnames and remove Duplicates column
publication_table <- qc_table %>%
  rename(
    Unmapped = unmapped_aln,
    Secondary = secondary_aln,
    Supplementary = supplementary_aln,
    Duplicates = duplicate_aln,
    `Primary mapped non-duplicate` = primary_nondup_mapped_aln,
    `MAPQ ≥ 20` = mapq_ge20_aln,
    `Junction-spanning reads (≥10 bp each side)` = junction_overhang10_reads
  ) %>%
  select(-Duplicates) # Duplicates are 0, due to Picard de-duplication

publication_table$sample <- sample_display_name(publication_table$sample)

pub_out_file <- sub("\\.tsv$", "_publication.tsv", OUT_FILE)
if (identical(pub_out_file, OUT_FILE)) {
  pub_out_file <- paste0(OUT_FILE, "_publication.tsv")
}

write.table(
  publication_table,
  file = pub_out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("Wrote:\n  ", pub_out_file, "\n")

# create LaTeX row output for manual booktabs table assembly
latex_integerize <- function(df) {
  out <- df
  for (nm in names(out)) {
    if (nm == "sample") next
    vals <- suppressWarnings(as.numeric(out[[nm]]))
    if (all(!is.na(vals))) {
      out[[nm]] <- as.character(as.integer(round(vals)))
    }
  }
  out
}

latex_safe <- function(x) {
  x <- ifelse(is.na(x), "NA", as.character(x))
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#$])", "\\\\\\1", x)
  x
}

publication_latex_table <- publication_table
publication_latex_table <- latex_integerize(publication_latex_table)

latex_rows <- apply(publication_latex_table, 1, function(row_vals) {
  paste0(paste(latex_safe(row_vals), collapse = " & "), "\\\\")
})

latex_out_file <- sub("\\.tsv$", "_rows.tex", pub_out_file)
if (identical(latex_out_file, pub_out_file)) {
  latex_out_file <- paste0(pub_out_file, "_rows.tex")
}
writeLines(latex_rows, con = latex_out_file, useBytes = TRUE)

cat("Wrote:\n  ", latex_out_file, "\n")
