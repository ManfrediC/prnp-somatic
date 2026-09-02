#!/usr/bin/env Rscript

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Run this test with Rscript.", call. = FALSE)
  }
  dirname(normalizePath(
    sub("^--file=", "", file_arg[[1]]),
    winslash = "/",
    mustWork = TRUE
  ))
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

assert_number <- function(actual, expected, label, tolerance = 1e-15) {
  assert_true(
    length(actual) == 1L && !is.na(actual) &&
      isTRUE(all.equal(actual, expected, tolerance = tolerance)),
    paste0(label, " changed: expected ", expected, ", found ", actual)
  )
}

read_tsv <- function(path) {
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character"
  )
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
cjd_dir <- file.path(
  repo_root,
  "results2", "sequencing2", "results", "mutect2_cjd_dilutions_with_pon",
  "variant_qc", "cjd"
)
controls_dir <- file.path(
  repo_root,
  "results2", "sequencing2", "results", "mutect2_controls_no_pon", "variant_qc"
)
dilutions_dir <- file.path(
  repo_root,
  "results2", "sequencing2", "results", "mutect2_cjd_dilutions_with_pon",
  "variant_qc", "dilutions"
)

candidates <- read_tsv(file.path(cjd_dir, "filtered_prnp_variants.tsv"))
candidate_keys <- paste(
  candidates$sample,
  candidates$CHROM,
  candidates$POS,
  candidates$REF,
  candidates$ALT,
  sep = "|"
)
expected_keys <- c(
  "CJD2|chr20|4691920|G|A",
  "CJD23|chr20|4691920|G|A",
  "CJD23|chr20|4694249|T|C"
)
assert_true(identical(sort(candidate_keys), sort(expected_keys)), "CJD call membership changed.")
assert_true(nrow(candidates) == 3L, "CJD call count changed from three.")
assert_true(
  nrow(unique(candidates[c("CHROM", "POS", "REF", "ALT")])) == 2L,
  "CJD unique-site count changed from two."
)

expected <- data.frame(
  key = expected_keys,
  AAF = c(0.008864499788940482, 0.009463113171108536, 0.008172362555720654),
  DP = c(4738, 5178, 1346),
  REF_count = c(4696, 5129, 1335),
  ALT_count = c(42, 49, 11),
  SB_refF = c(2336, 2538, 687),
  SB_refR = c(2360, 2591, 648),
  SB_altF = c(14, 17, 5),
  SB_altR = c(28, 32, 6),
  MEAN_BQ = c(22.01, 21.00, 20.62),
  MEAN_MQ = c(37.90, 38.49, 30.01)
)
observed <- candidates[match(expected$key, candidate_keys), , drop = FALSE]
for (column in setdiff(names(expected), "key")) {
  assert_true(
    isTRUE(all.equal(as.numeric(observed[[column]]), expected[[column]], tolerance = 1e-15)),
    paste0("CJD ", column, " values changed.")
  )
}

settings <- read_tsv(file.path(cjd_dir, "run_settings.tsv"))
settings <- setNames(settings$value, settings$key)
assert_true(settings[["enable_aaf_filter"]] == "TRUE", "CJD AAF filter is not enabled.")
assert_true(settings[["aaf_filter_applied"]] == "TRUE", "CJD AAF filter was not applied.")
assert_true(
  settings[["aaf_threshold"]] == "0.006682086867129272",
  "CJD AAF threshold changed."
)
assert_number(as.numeric(settings[["aaf_threshold"]]), 26 / 3891, "CJD exact LoD", 1e-18)

lod <- read_tsv(file.path(repo_root, "results2", "spikein", "read_recovery", "empirical_lod.tsv"))
assert_true(nrow(lod) == 1L, "Empirical LoD table row count changed.")
assert_true(lod$bam_readcount_fraction == "26/3891", "Empirical LoD fraction changed.")
assert_true(lod$bam_readcount_vaf == settings[["aaf_threshold"]], "LoD sources disagree.")
assert_true(lod$filtermutectcalls_status == "PASS", "Empirical LoD call is not PASS.")

controls <- read_tsv(file.path(controls_dir, "final_noPoN_variants.tsv"))
assert_true(nrow(controls) == 0L, "The final controls table is no longer empty.")
dilution_settings <- read_tsv(file.path(dilutions_dir, "run_settings.tsv"))
dilution_settings <- setNames(dilution_settings$value, dilution_settings$key)
assert_true(dilution_settings[["enable_aaf_filter"]] == "FALSE", "Dilution AAF filter is enabled.")
assert_true(dilution_settings[["aaf_filter_applied"]] == "FALSE", "Dilution AAF filter was applied.")
assert_true(
  dilution_settings[["aaf_threshold"]] == settings[["aaf_threshold"]],
  "Dilution recorded threshold differs from CJD."
)

figure_csv <- file.path(
  repo_root, "manuscript2", "figures2", "snv_lollipop", "SNV_lollipop_data.csv"
)
figure_data <- read.csv(figure_csv, stringsAsFactors = FALSE, check.names = FALSE)
figure_keys <- paste(
  figure_data$sample,
  "chr20",
  figure_data$position,
  figure_data$ref,
  figure_data$alt,
  sep = "|"
)
assert_true(identical(sort(figure_keys), sort(expected_keys)), "Figure membership differs from Table 6.")
assert_true(
  isTRUE(all.equal(sort(figure_data$vaf), sort(100 * as.numeric(candidates$AAF)), tolerance = 1e-12)),
  "Figure VAFs differ from the final CJD table."
)
figure_svg <- readLines(
  file.path(repo_root, "manuscript2", "figures2", "snv_lollipop", "SNV_lollipop.svg"),
  warn = FALSE
)
assert_true(!any(grepl("CJD6", figure_svg, fixed = TRUE)), "Figure still contains CJD6.")
assert_true(any(grepl("LoD (0.67%)", figure_svg, fixed = TRUE)), "Figure LoD label changed.")
assert_true(any(grepl("figure8-y-axis-label", figure_svg, fixed = TRUE)), "Figure lacks the reviewed y-axis label.")

table_script <- file.path(
  repo_root, "manuscript2", "tables2", "make_prnp_somatic_snv_summary_tex.py"
)
table_output <- file.path(
  repo_root, "manuscript2", "tables2", "table_prnp_somatic_snv_summary.tex"
)
temporary_table <- tempfile(fileext = ".tex")
temporary_preview <- tempfile(fileext = ".tex")
on.exit(unlink(c(temporary_table, temporary_preview)), add = TRUE)
status <- system2(
  "python3",
  c(
    table_script,
    "--output", temporary_table,
    "--preview", temporary_preview
  )
)
assert_true(status == 0L, "Temporary Table 6 generation failed.")
assert_true(
  identical(
    readBin(temporary_table, "raw", file.info(temporary_table)$size),
    readBin(table_output, "raw", file.info(table_output)$size)
  ),
  "Generated Table 6 differs from the retained Table 6."
)

spikein_tex <- readLines(
  file.path(repo_root, "manuscript2", "tables2", "table_spikein_snp_recovery.tex"),
  warn = FALSE
)
assert_true(any(grepl("26/3891 \\(0\\.668\\\\%\\)", spikein_tex)), "Spike-in LoD caption changed.")
assert_true(
  any(grepl("Mutect2-estimated VAF", spikein_tex, fixed = TRUE)),
  "Table 4 Mutect2-estimated VAF heading changed."
)

message("Sequencing2 downstream provenance regression test passed.")
