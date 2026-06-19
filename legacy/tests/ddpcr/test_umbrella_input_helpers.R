source("src/ddpcr/umbrella/umbrella_input_helpers.R")

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

# A positive D178N LoD well exercises all four cluster labels, which is the
# fragile Umbrella-specific part that dPCP does not need.
lod_archive <- file.path(project_root, "raw", "ddpcr_lod", "d178n", "archive_contents")
lod_well <- read_umbrella_physical_well(
  archive_dir = lod_archive,
  well = "C01",
  assay = "D178N"
)

expect_equal(names(lod_well$data), c("Ch1.Amplitude", "Ch2.Amplitude", "Cluster"), "Umbrella data columns")
expect_true(all(vapply(lod_well$data[1:2], is.numeric, logical(1))), "Umbrella amplitude columns numeric")
expect_equal(sort(unique(lod_well$data$Cluster)), 1:4, "Umbrella cluster levels")
expect_equal(sum(lod_well$partition_counts), nrow(lod_well$data), "partition count sum")

# M02 in the D178N LoD manifest merges C01 and C02. The merged Umbrella object
# must keep row-level clusters while reconciling to the manifest accepted count.
lod_merged <- combine_umbrella_physical_wells(
  archive_dir = lod_archive,
  wells = c("C01", "C02"),
  assay = "D178N"
)

expect_equal(nrow(lod_merged$data), 22879L, "LoD merged Umbrella accepted droplets")
expect_equal(sort(unique(lod_merged$data$Cluster)), 1:4, "LoD merged Umbrella cluster levels")
expect_true(sum(lod_merged$partition_counts) == 22879L, "LoD merged partition counts sum")

# A clean NTC well can have only the double-negative cluster; this is valid for
# an individual partition set even though the full plate list must contain four
# clusters somewhere if Umbrella is to use prior cluster labels.
experimental_archive <- file.path(
  project_root,
  "raw", "ddpcr", "archive_contents", "2020-10-12", "SNV_D178N"
)
experimental_ntc <- read_umbrella_physical_well(
  archive_dir = experimental_archive,
  well = "A01",
  assay = "D178N"
)

expect_equal(nrow(experimental_ntc$data), 14654L, "experimental Umbrella accepted droplets")
expect_equal(unique(experimental_ntc$data$Cluster), 1L, "experimental NTC cluster")

message("Umbrella input helper regression test passed")
