#!/usr/bin/env Rscript
# Controls-only variant-table QC entrypoint (no PoN).

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[1]),
      winslash = "/",
      mustWork = TRUE
    )))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
source(file.path(script_dir, "variant_table_qc_common.R"))

cli_args <- parse_variant_table_qc_cli(commandArgs(trailingOnly = TRUE))
do.call(
  run_variant_table_qc,
  c(cli_args, list(output_prefix = "noPoN"))
)
