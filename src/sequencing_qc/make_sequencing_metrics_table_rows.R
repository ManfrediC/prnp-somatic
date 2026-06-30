#!/usr/bin/env Rscript

# Build LaTeX data rows for the sequencing metrics manuscript table.

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

## ---- configure inputs and generated row output ----
metrics_tsv <- Sys.getenv(
  "PRNP_SEQUENCING_METRICS_TSV",
  unset = file.path(repo_root, "results", "sequencing_qc", "sequencing_metrics_per_sample.tsv")
)

rows_tex <- Sys.getenv(
  "PRNP_SEQUENCING_METRICS_ROWS",
  unset = file.path(repo_root, "results", "sequencing_qc", "sequencing_metrics_table_rows.tex")
)

sample_manifest_tsv <- Sys.getenv(
  "PRNP_MANIFEST",
  unset = file.path(repo_root, "config", "sample_manifest.tsv")
)

## ---- define the manuscript-facing spike-in rows ----
manuscript_spike_ids <- c("NA100_undil", "NA995A05_undil", "NA99A1_undil")

load_spike_display_names <- function(manifest_path, sample_ids) {
  if (!file.exists(manifest_path)) {
    stop("Sample manifest does not exist: ", manifest_path)
  }

  manifest <- read.delim(
    manifest_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c("sample_id", "display_label")
  missing_cols <- setdiff(required_cols, names(manifest))
  if (length(missing_cols) > 0L) {
    stop("Sample manifest is missing required display-label column(s): ", paste(missing_cols, collapse = ", "))
  }

  rows <- manifest$sample_id %in% sample_ids
  if (any(duplicated(manifest$sample_id[rows]))) {
    stop("Sample manifest has duplicate rows for manuscript spike-in IDs")
  }

  labels <- trimws(as.character(manifest$display_label[rows]))
  display_names <- setNames(labels, manifest$sample_id[rows])
  display_names <- display_names[sample_ids]

  missing_labels <- names(display_names)[is.na(display_names) | display_names == "" | display_names == "NA"]
  if (length(missing_labels) > 0L) {
    stop("Sample manifest is missing display_label values for: ", paste(missing_labels, collapse = ", "))
  }

  display_names
}

spike_display_names <- load_spike_display_names(sample_manifest_tsv, manuscript_spike_ids)

## ---- load the canonical sequencing metrics table ----
metrics <- read.delim(
  metrics_tsv,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

## ---- validate the columns needed by the publication table ----
required_cols <- c(
  "sample_id",
  "group",
  "reads_primary_mapped",
  "on_target_percent",
  "coding_mean_depth",
  "coding_p20_depth",
  "coding_fold80",
  "coding_pct_bases_ge_x1",
  "coding_pct_bases_ge_x2"
)

missing_cols <- setdiff(required_cols, names(metrics))
if (length(missing_cols) > 0L) {
  stop("Sequencing metrics table is missing required columns: ", paste(missing_cols, collapse = ", "))
}

## ---- keep the same sample set shown in the manuscript table ----
keep_row <- metrics$group %in% c("CJD", "control") | metrics$sample_id %in% names(spike_display_names)
publication_rows <- metrics[keep_row, required_cols, drop = FALSE]

## ---- format internal sample IDs for manuscript display ----
sample_display_name <- function(sample_ids) {
  out <- sample_ids
  spike_idx <- match(sample_ids, names(spike_display_names))
  has_spike_name <- !is.na(spike_idx)
  out[has_spike_name] <- unname(spike_display_names[spike_idx[has_spike_name]])
  sub("^Ctrl", "Control", out)
}

## ---- format group labels for the selected manuscript rows ----
group_display_name <- function(sample_ids, groups) {
  out <- groups
  out[sample_ids %in% names(spike_display_names)] <- "spike-in"
  out
}

## ---- sort as CJD, Control, then spike-in rows ----
sample_number <- function(sample_ids) {
  out <- suppressWarnings(as.integer(sub("^[A-Za-z]+", "", sample_ids)))
  out[is.na(out)] <- 999999L
  out
}

row_group_rank <- ifelse(
  grepl("^CJD[0-9]+$", publication_rows$sample_id),
  1L,
  ifelse(grepl("^Ctrl[0-9]+$", publication_rows$sample_id), 2L, 3L)
)

row_number_rank <- sample_number(publication_rows$sample_id)
spike_rank <- match(publication_rows$sample_id, names(spike_display_names))
row_number_rank[row_group_rank == 3L] <- spike_rank[row_group_rank == 3L]

row_order <- order(row_group_rank, row_number_rank, publication_rows$sample_id)
publication_rows <- publication_rows[row_order, , drop = FALSE]
row_group_rank <- row_group_rank[row_order]

## ---- format numeric values exactly as used in the manuscript table ----
fmt_int <- function(x) {
  vals <- suppressWarnings(as.numeric(x))
  ifelse(is.na(vals), "NA", as.character(as.integer(round(vals))))
}

fmt_one_dp <- function(x) {
  vals <- suppressWarnings(as.numeric(x))
  ifelse(is.na(vals), "NA", sprintf("%.1f", round(vals, 1)))
}

fmt_two_dp <- function(x) {
  vals <- suppressWarnings(as.numeric(x))
  ifelse(is.na(vals), "NA", sprintf("%.2f", round(vals, 2)))
}

## ---- assemble the manuscript table columns in display order ----
publication_table <- data.frame(
  Sample = sample_display_name(publication_rows$sample_id),
  Group = group_display_name(publication_rows$sample_id, publication_rows$group),
  `Primary mapped (n)` = fmt_int(publication_rows$reads_primary_mapped),
  `On-target (%)` = fmt_one_dp(publication_rows$on_target_percent),
  `Mean depth (x)` = fmt_int(publication_rows$coding_mean_depth),
  `P20 depth (x)` = fmt_int(publication_rows$coding_p20_depth),
  `Fold-80` = fmt_two_dp(publication_rows$coding_fold80),
  `Bases >=100x (%)` = fmt_one_dp(publication_rows$coding_pct_bases_ge_x1),
  `Bases >=500x (%)` = fmt_one_dp(publication_rows$coding_pct_bases_ge_x2),
  check.names = FALSE
)

## ---- escape data cells for safe LaTeX row output ----
latex_safe <- function(x) {
  x <- ifelse(is.na(x), "NA", as.character(x))
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#$])", "\\\\\\1", x)
  x
}

## ---- create data rows and keep visual spacing between table sections ----
latex_rows <- character(0)
section_rank <- row_group_rank
for (i in seq_len(nrow(publication_table))) {
  if (i > 1L && section_rank[i] != section_rank[i - 1L]) {
    latex_rows <- c(latex_rows, "\\addlinespace")
  }

  latex_rows <- c(
    latex_rows,
    paste0(paste(latex_safe(unlist(publication_table[i, ], use.names = FALSE)), collapse = " & "), "\\\\")
  )
}

## ---- write the generated rows for manual booktabs table assembly ----
dir.create(dirname(rows_tex), showWarnings = FALSE, recursive = TRUE)
writeLines(latex_rows, con = rows_tex, useBytes = TRUE)

cat("Wrote:\n  ", rows_tex, "\n")
