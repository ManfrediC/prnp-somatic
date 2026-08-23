#!/usr/bin/env Rscript

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) {
    stop("Run this test with Rscript.", call. = FALSE)
  }
  dirname(normalizePath(
    sub("^--file=", "", file_arg[[1]]),
    winslash = "/",
    mustWork = TRUE
  ))
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

assert_number <- function(actual, expected, label, tolerance = 1e-15) {
  assert_true(
    length(actual) == 1L && !is.na(actual) &&
      isTRUE(all.equal(actual, expected, tolerance = tolerance)),
    paste0(label, " changed: expected ", expected, ", found ", actual)
  )
}

script_dir <- get_script_dir()
repo_root <- normalizePath(
  file.path(script_dir, "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(repo_root, "src", "pipelines", "variant_table_qc_common.R"))

manual_freq <- readr::read_tsv(
  file.path(repo_root, "resources", "annotations", "manual_population_freq.tsv"),
  comment = "#",
  na = c("", "NA"),
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)
curated_fixture <- curate_dbsnp_ids(
  c("rs996098774", "rs1012397146"),
  manual_freq
)
assert_true(is.na(curated_fixture[[1]]), "Invalid raw rsID was not curated to NA.")
assert_true(
  identical(curated_fixture[[2]], "rs1012397146"),
  "A valid raw rsID was changed by curation."
)

candidate_path <- file.path(
  repo_root,
  "results", "mutect2_cjd_dilutions_with_pon", "variant_qc", "cjd",
  "filtered_prnp_variants.tsv"
)
candidates <- readr::read_tsv(
  candidate_path,
  na = c("", "NA"),
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)
assert_true("dbsnp_id_raw" %in% colnames(candidates), "dbsnp_id_raw is missing.")
assert_true(nrow(candidates) == 4L, "Candidate row count changed from four.")
assert_true(
  nrow(dplyr::distinct(candidates, CHROM, POS, REF, ALT)) == 2L,
  "Candidate site count changed from two."
)

target <- candidates %>%
  dplyr::filter(CHROM == "chr20", POS == "4694249", REF == "T", ALT == "C")
assert_true(nrow(target) == 1L, "chr20:4694249 T>C is not present exactly once.")
assert_true(
  identical(target$dbsnp_id_raw[[1]], "rs996098774"),
  "chr20:4694249 T>C raw rsID changed."
)
assert_true(
  identical(
    target$dbsnp_id_raw[[1]],
    stringr::str_extract(target$FUNCOTATION_PRIMARY[[1]], "rs\\d+")
  ),
  "Raw rsID is not the direct FUNCOTATION_PRIMARY extraction."
)
assert_true(is.na(target$dbsnp_id[[1]]), "chr20:4694249 T>C curated rsID is not NA.")
assert_true(
  is.na(target$population_frequency[[1]]),
  "chr20:4694249 T>C population frequency is not NA."
)
assert_true(identical(target$FILTER[[1]], "PASS"), "Target filtering result changed.")
assert_number(as.numeric(target$AAF), 0.00818452380952381, "Target VAF")
assert_number(as.numeric(target$DP), 1344, "Target depth")
assert_number(as.numeric(target$REF_count), 1333, "Target reference count")
assert_number(as.numeric(target$ALT_count), 11, "Target alternate count")
assert_number(as.numeric(target$SB_refF), 686, "Target forward reference count")
assert_number(as.numeric(target$SB_refR), 647, "Target reverse reference count")
assert_number(as.numeric(target$SB_altF), 5, "Target forward alternate count")
assert_number(as.numeric(target$SB_altR), 6, "Target reverse alternate count")

comparison <- candidates %>%
  dplyr::filter(CHROM == "chr20", POS == "4691920", REF == "G", ALT == "A") %>%
  dplyr::arrange(sample)
assert_true(nrow(comparison) == 3L, "chr20:4691920 G>A candidate rows changed.")
assert_true(
  identical(comparison$sample, c("CJD2", "CJD23", "CJD6")),
  "chr20:4691920 G>A sample membership changed."
)
assert_true(all(is.na(comparison$dbsnp_id_raw)), "Comparison-site raw rsID changed.")
assert_true(all(is.na(comparison$dbsnp_id)), "Comparison-site curated rsID changed.")
assert_true(
  all(is.na(comparison$population_frequency)),
  "Comparison-site population frequency changed."
)
assert_true(all(comparison$FILTER == "PASS"), "Comparison-site filtering result changed.")

expected <- data.frame(
  sample = c("CJD2", "CJD23", "CJD6"),
  AAF = c(0.0055248618784530384, 0.006039353204753555, 0.005169442848937392),
  DP = c(4706, 5133, 3482),
  REF_count = c(4680, 5102, 3464),
  ALT_count = c(26, 31, 18),
  SB_refF = c(2327, 2525, 1787),
  SB_refR = c(2353, 2577, 1677),
  SB_altF = c(10, 12, 6),
  SB_altR = c(16, 19, 12)
)
for (column in setdiff(colnames(expected), "sample")) {
  actual <- as.numeric(comparison[[column]])
  wanted <- expected[[column]]
  assert_true(
    isTRUE(all.equal(actual, wanted, tolerance = 1e-15)),
    paste0("chr20:4691920 G>A ", column, " values changed.")
  )
}

table_script <- file.path(
  repo_root, "manuscript", "tables", "main",
  "make_table_prnp_somatic_snv_summary_tex.py"
)
table_canonical <- file.path(
  repo_root, "manuscript", "tables", "main",
  "table_prnp_somatic_snv_summary.tex"
)
table_comparison <- tempfile(fileext = ".tex")
on.exit(unlink(table_comparison), add = TRUE)
status <- system2(
  "python3",
  c(
    shQuote(table_script),
    "--input", shQuote(candidate_path),
    "--output", shQuote(table_comparison)
  )
)
assert_true(status == 0L, "Temporary Table 6 generation failed.")
assert_true(
  identical(readBin(table_comparison, "raw", file.info(table_comparison)$size),
            readBin(table_canonical, "raw", file.info(table_canonical)$size)),
  "Generated Table 6 content differs from the canonical Table 6."
)

message("Variant-table provenance regression test passed.")
