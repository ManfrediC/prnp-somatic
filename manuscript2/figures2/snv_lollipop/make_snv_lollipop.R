#!/usr/bin/env Rscript

# Build the sequencing2 Figure 8 through the original reviewed SVG workflow.

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Run this script with Rscript.", call. = FALSE)
  }
  dirname(normalizePath(
    sub("^--file=", "", file_arg[[1]]),
    winslash = "/",
    mustWork = TRUE
  ))
}

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Input does not exist: ", path, call. = FALSE)
  }
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character"
  )
}

require_columns <- function(data, columns, source) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop(source, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), winslash = "/", mustWork = TRUE)
cjd_dir <- file.path(
  repo_root,
  "results2", "sequencing2", "results", "mutect2_cjd_dilutions_with_pon",
  "variant_qc", "cjd"
)
calls_path <- file.path(cjd_dir, "filtered_prnp_variants.tsv")
settings_path <- file.path(cjd_dir, "run_settings.tsv")
lod_path <- file.path(repo_root, "results2", "spikein", "read_recovery", "empirical_lod.tsv")

calls <- read_tsv(calls_path)
settings <- read_tsv(settings_path)
lod <- read_tsv(lod_path)
require_columns(
  calls,
  c("sample", "CHROM", "POS", "REF", "ALT", "AAF", "gene", "mutation_type"),
  "CJD final PRNP table"
)
require_columns(settings, c("key", "value"), "CJD run settings")
require_columns(
  lod,
  c("bam_readcount_fraction", "bam_readcount_vaf", "filtermutectcalls_status"),
  "Empirical LoD table"
)

expected_keys <- c(
  "CJD2|chr20|4691920|G|A",
  "CJD23|chr20|4691920|G|A",
  "CJD23|chr20|4694249|T|C"
)
observed_keys <- paste(calls$sample, calls$CHROM, calls$POS, calls$REF, calls$ALT, sep = "|")
if (!identical(sort(observed_keys), sort(expected_keys))) {
  stop("Unexpected sequencing2 CJD membership: ", paste(observed_keys, collapse = ", "), call. = FALSE)
}
calls <- calls[match(expected_keys, observed_keys), , drop = FALSE]

threshold_text <- settings$value[settings$key == "aaf_threshold"]
valid_lod <- c(
  length(threshold_text) == 1L,
  threshold_text == "0.006682086867129272",
  settings$value[settings$key == "enable_aaf_filter"] == "TRUE",
  settings$value[settings$key == "aaf_filter_applied"] == "TRUE",
  nrow(lod) == 1L,
  lod$bam_readcount_fraction == "26/3891",
  lod$bam_readcount_vaf == threshold_text,
  lod$filtermutectcalls_status == "PASS"
)
if (!all(valid_lod)) {
  stop("Sequencing2 LoD provenance does not match 26/3891 PASS", call. = FALSE)
}

# Define SNVs in the same order and units used by the original R script.
positions <- as.integer(calls$POS)
samples <- calls$sample
aafs <- as.numeric(calls$AAF) * 100
refs <- calls$REF
alts <- calls$ALT

data_path <- file.path(script_dir, "SNV_lollipop_data.csv")
write.csv(
  data.frame(
    sample = samples,
    position = positions,
    vaf = aafs,
    ref = refs,
    alt = alts,
    stringsAsFactors = FALSE
  ),
  data_path,
  row.names = FALSE,
  quote = TRUE
)

python_script <- file.path(script_dir, "build_snv_lollipop_svg.py")
baseline_svg <- file.path(script_dir, "SNV_lollipop_baseline.svg")
reviewed_template_svg <- file.path(script_dir, "SNV_lollipop_template.svg")
for (source in c(python_script, baseline_svg, reviewed_template_svg)) {
  if (!file.exists(source) || file.info(source)$size <= 0) {
    stop("Lollipop source is missing or empty: ", source, call. = FALSE)
  }
}

python_exec <- Sys.which("python3")
if (!nzchar(python_exec)) {
  python_exec <- Sys.which("python")
}
if (!nzchar(python_exec)) {
  stop("Python executable not found", call. = FALSE)
}
inkscape_exec <- Sys.getenv(
  "INKSCAPE_PATH",
  unset = "/mnt/c/Program Files/Inkscape/bin/inkscape.com"
)

lod_percent <- as.numeric(threshold_text) * 100
status <- system2(
  python_exec,
  c(
    python_script,
    "--input", data_path,
    "--baseline", baseline_svg,
    "--reviewed-template", reviewed_template_svg,
    "--output-dir", script_dir,
    "--lod-percent", sprintf("%.15f", lod_percent),
    "--inkscape", shQuote(inkscape_exec)
  )
)
if (!identical(status, 0L)) {
  stop("Lollipop baseline/diff renderer failed with exit status ", status, call. = FALSE)
}

required_outputs <- c(
  data_path,
  file.path(script_dir, "SNV_lollipop.svg"),
  file.path(script_dir, "SNV_lollipop.pdf"),
  file.path(script_dir, "SNV_lollipop_edit_diff.json")
)
missing_outputs <- required_outputs[
  !file.exists(required_outputs) |
    is.na(file.info(required_outputs)$size) |
    file.info(required_outputs)$size <= 0
]
if (length(missing_outputs) > 0L) {
  stop("Lollipop outputs are missing or empty: ", paste(missing_outputs, collapse = ", "))
}

caption <- paste0(
  "\\begin{figure}[p]\n",
  "  \\centering\n",
  "  \\includegraphics[width=0.9\\linewidth]{manuscript2/figures2/snv_lollipop/SNV_lollipop.pdf}\n",
  "  \\caption{\\textbf{PRNP intronic calls retained by sequencing2.} ",
  "Three sample-level calls representing two genomic sites passed all configured filters. ",
  "The empirical complete-pipeline LoD was 26/3891 (0.668\\%). ",
  "CJD23 and CJD2 carried chr20:4691920 G\\textgreater{}A at VAFs of 0.95\\% and 0.89\\%, ",
  "respectively; CJD23 also carried chr20:4694249 T\\textgreater{}C at 0.82\\%. ",
  "Horizontal axis: GRCh38 genomic position. Vertical axis: VAF (\\%).}\n",
  "  \\label{fig:sequencing2_snv_lollipop}\n",
  "\\end{figure}\n"
)
writeLines(caption, file.path(script_dir, "figure_snv_lollipop_variants.tex"), useBytes = TRUE)

message("Wrote: ", paste(required_outputs, collapse = ", "))
