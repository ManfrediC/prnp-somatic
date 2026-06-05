source("src/ddpcr/ddpcr_raw_import_helpers.R")

expect_close <- function(actual, expected, tolerance = 1e-9, label = deparse(substitute(actual))) {
  if (is.na(actual) || abs(actual - expected) > tolerance) {
    stop(label, " was ", actual, ", expected ", expected, call. = FALSE)
  }
}

fa <- fractional_abundance(
  ref_positive = 1000,
  ref_negative = 9000,
  mut_positive = 10,
  mut_negative = 9990,
  total = 10000
)

expected_lambda_ref <- -log1p(-1000 / 10000)
expected_lambda_mut <- -log1p(-10 / 10000)
expected_fa_percent <- 100 * expected_lambda_mut / (expected_lambda_ref + expected_lambda_mut)

expect_close(fa$lambda_ref, expected_lambda_ref, label = "lambda_ref")
expect_close(fa$lambda_mut, expected_lambda_mut, label = "lambda_mut")
expect_close(fa$fractional_abundance, expected_fa_percent, label = "fractional_abundance")

if (fa$fractional_abundance <= 0 || fa$fractional_abundance >= 100) {
  stop("fractional_abundance should be a percentage bounded by 0 and 100", call. = FALSE)
}

helper_text <- readLines("src/ddpcr/ddpcr_raw_import_helpers.R", warn = FALSE)
retired_markers <- c(
  "quantasoft_count_bound_table",
  "quanta_like_fractional_abundance",
  "quanta_like_ci_low",
  "quanta_like_ci_high"
)
for (marker in retired_markers) {
  if (any(grepl(marker, helper_text, fixed = TRUE))) {
    stop("Retired helper marker is still present: ", marker, call. = FALSE)
  }
}

message("ddPCR fractional abundance regression test passed")
