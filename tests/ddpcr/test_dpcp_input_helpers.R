source("src/ddpcr/dpcp/dpcp_input_helpers.R")

expect_equal <- function(actual, expected, label = deparse(substitute(actual))) {
  if (!identical(actual, expected)) {
    stop(label, " was ", paste(actual, collapse = ", "), ", expected ", paste(expected, collapse = ", "), call. = FALSE)
  }
}

expect_true <- function(value, label = deparse(substitute(value))) {
  if (!isTRUE(value)) {
    stop(label, " was not TRUE", call. = FALSE)
  }
}

project_root <- get_dpcp_project_root()

# A01 in the 2020-10-12 D178N run is an NTC well with 14,654 accepted droplets.
# The metadata also contains gated droplets, so this fixture verifies that the
# helper exports the accepted cluster only rather than every peak in PeakData.
experimental_archive <- file.path(
  project_root,
  "raw", "ddpcr", "archive_contents", "2020-10-12", "SNV_D178N"
)
experimental_well <- read_dpcp_physical_well(
  archive_dir = experimental_archive,
  well = "A01",
  assay = "D178N"
)

expect_equal(nrow(experimental_well$data), 14654L, "experimental accepted droplets")
expect_equal(names(experimental_well$data), c("Ch1.Amplitude", "Ch2.Amplitude"), "amplitude column names")
expect_true(all(vapply(experimental_well$data, is.numeric, logical(1))), "experimental amplitude columns numeric")
expect_equal(experimental_well$ch1_target, "WT", "experimental Ch1 target")
expect_equal(experimental_well$ch2_target, "D178N", "experimental Ch2 target")
expect_equal(as.integer(experimental_well$partition_counts[["Ch1-Ch2-"]]), 14654L, "experimental empty droplets")

# M02 in the D178N LoD manifest merges C01 and C02. The dPCP compatibility file
# must therefore concatenate the accepted droplets from both physical wells and
# reconcile to the merged accepted count recorded in the LoD manifest.
lod_archive <- file.path(project_root, "raw", "ddpcr_lod", "d178n", "archive_contents")
lod_merged <- combine_dpcp_physical_wells(
  archive_dir = lod_archive,
  wells = c("C01", "C02"),
  assay = "D178N"
)

expect_equal(nrow(lod_merged$data), 22879L, "LoD merged accepted droplets")
expect_equal(lod_merged$accepted_count, 22879L, "LoD accepted_count")
expect_equal(lod_merged$ch1_target, "WT", "LoD Ch1 target")
expect_equal(lod_merged$ch2_target, "D178N", "LoD Ch2 target")
expect_true(sum(lod_merged$partition_counts) == 22879L, "LoD partition counts sum")

sample_table <- make_dpcp_sample_table(
  tibble(
    sample = "D178N_1",
    dpcp_well_id = "20200911_SNV_D178N_M02",
    ch1_target = lod_merged$ch1_target,
    ch2_target = lod_merged$ch2_target
  ),
  reference_filename = "20200911_SNV_D178N_M02_Amplitude.csv"
)

expect_equal(names(sample_table), dpcp_sample_table_columns, "sample table columns")
expect_equal(sample_table$`Chip ID/Well ID`[[1]], "20200911_SNV_D178N_M02", "sample table well ID")

message("dPCP input helper regression test passed")
