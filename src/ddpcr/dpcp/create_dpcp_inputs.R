#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

# This script creates dPCP-compatible input folders from the reviewed ddPCR raw
# database. dPCP does not consume the repository's QuantaSoft/QX Manager summary
# CSVs; it needs one droplet-amplitude CSV per sample/reference plus a sample
# table. The raw archive JSONs are therefore the source of truth for amplitudes,
# while the manifests remain the source of truth for which wells are active.
command_line <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_line, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  file.path(getwd(), "src", "ddpcr", "dpcp")
}
source(file.path(script_dir, "dpcp_input_helpers.R"))
project_root <- get_dpcp_project_root()

args <- commandArgs(trailingOnly = TRUE)

# The CLI is intentionally small because this converter is a reproducibility
# bridge, not a general import framework. By default it refreshes every dPCP
# input artefact; --dataset narrows that to one source while keeping the same
# validation code path.
arg_or_default <- function(pattern, default) {
  matches <- args[str_detect(args, pattern)]
  if (length(matches) == 0L || is.na(matches[[1]])) {
    return(default)
  }
  matches[[1]]
}

dataset_arg <- arg_or_default("^--dataset=", "--dataset=all")
dataset <- str_remove(dataset_arg, "^--dataset=")
if (!dataset %in% c("all", "experimental", "lod")) {
  stop("--dataset must be one of: all, experimental, lod", call. = FALSE)
}

min_ref_arg <- arg_or_default("^--min-reference-cluster=", "--min-reference-cluster=50")
min_reference_cluster <- as.integer(str_remove(min_ref_arg, "^--min-reference-cluster="))
if (is.na(min_reference_cluster) || min_reference_cluster < 1L) {
  stop("--min-reference-cluster must be a positive integer.", call. = FALSE)
}

experimental_raw_root <- file.path(project_root, "raw", "ddpcr")
lod_raw_root <- file.path(project_root, "raw", "ddpcr_lod")
input_root <- file.path(project_root, "results", "ddPCR", "dpcp", "inputs")
validation_root <- file.path(project_root, "results", "ddPCR", "dpcp", "validation")
analysis_root <- file.path(project_root, "results", "ddPCR", "dpcp", "analysis")

dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
dir.create(validation_root, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)

# Write one amplitude CSV and one manifest row for an already-extracted dPCP
# well record. This small shared writer keeps experimental physical wells and
# LoD merged pseudo-wells on the same validation contract.
write_export_record <- function(
  dataset,
  assay,
  run_id,
  sample,
  source_wells,
  expected_accepted,
  out_dir,
  well_record,
  dpcp_well_id
) {
  amplitude_file <- paste0(dpcp_well_id, "_Amplitude.csv")
  amplitude_path <- file.path(out_dir, amplitude_file)

  write_dpcp_amplitude_csv(well_record$data, amplitude_path)

  tibble(
    dataset = dataset,
    assay = assay,
    run_id = run_id,
    sample = sample,
    source_wells = paste(source_wells, collapse = "|"),
    dpcp_well_id = dpcp_well_id,
    amplitude_file = amplitude_file,
    amplitude_path = amplitude_path,
    accepted_droplets_exported = nrow(well_record$data),
    accepted_droplets_manifest = expected_accepted,
    double_positive_droplets = as.integer(well_record$partition_counts[["Ch1+Ch2+"]]),
    ch1_only_droplets = as.integer(well_record$partition_counts[["Ch1+Ch2-"]]),
    ch2_only_droplets = as.integer(well_record$partition_counts[["Ch1-Ch2+"]]),
    empty_droplets = as.integer(well_record$partition_counts[["Ch1-Ch2-"]]),
    gated_or_unassigned_droplets = well_record$gated_or_unassigned,
    ch1_target = well_record$ch1_target,
    ch2_target = well_record$ch2_target,
    ch1_raw_target = well_record$ch1_raw_target,
    ch2_raw_target = well_record$ch2_raw_target,
    source_peak_files = well_record$source_peak_files %||% well_record$peak_path,
    source_metadata_files = well_record$source_metadata_files %||% well_record$metadata_path
  )
}

# The active experimental dataset already has one physical well per sample row
# in sample_manifest.csv. We use the manifest to decide which wells to export,
# and the archive contents to reconstruct the droplet amplitudes for those wells.
process_experimental_run <- function(run_row, sample_manifest, output_root) {
  run_id <- run_row$run_id[[1]]
  assay <- normalise_dpcp_assay(run_row$assay[[1]])
  archive_dir <- file.path(experimental_raw_root, run_row$archive_contents_relative_dir[[1]])
  out_dir <- file.path(output_root, run_id)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  run_manifest <- sample_manifest %>%
    filter(.data$run_id == !!run_id) %>%
    mutate(
      target_clean = clean_dpcp_target(target_clean),
      channel = as.character(channel)
    )

  well_rows <- run_manifest %>%
    distinct(well, sample, accepted_droplets) %>%
    arrange(well)

  export_manifest <- purrr::map_dfr(seq_len(nrow(well_rows)), function(i) {
    well_row <- well_rows[i, ]
    well <- well_row$well[[1]]

    well_record <- read_dpcp_physical_well(
      archive_dir = archive_dir,
      well = well,
      assay = assay
    )

    # The manifest and metadata should agree on the channel-to-target mapping.
    # A mismatch here would mean the dPCP sample table could silently swap WT
    # and mutant target labels, so we stop before writing misleading inputs.
    manifest_channels <- run_manifest %>%
      filter(.data$well == !!well) %>%
      distinct(channel, target_clean)
    manifest_ch1 <- manifest_channels %>% filter(channel == "Ch1") %>% pull(target_clean)
    manifest_ch2 <- manifest_channels %>% filter(channel == "Ch2") %>% pull(target_clean)
    if (length(unique(manifest_ch1)) != 1L || length(unique(manifest_ch2)) != 1L ||
        unique(manifest_ch1) != well_record$ch1_target ||
        unique(manifest_ch2) != well_record$ch2_target) {
      stop("Manifest channel mapping disagrees with archive metadata for ", run_id, " ", well, call. = FALSE)
    }

    write_export_record(
      dataset = "experimental",
      assay = assay,
      run_id = run_id,
      sample = well_row$sample[[1]],
      source_wells = well,
      expected_accepted = as.integer(well_row$accepted_droplets[[1]]),
      out_dir = out_dir,
      well_record = well_record,
      dpcp_well_id = make_dpcp_id(run_id, well)
    )
  })

  reference_candidates <- score_reference_candidates(
    export_manifest,
    min_partition_count = min_reference_cluster
  )
  reference_filename <- choose_reference_filename(reference_candidates)
  sample_table <- make_dpcp_sample_table(export_manifest, reference_filename)
  validation <- validate_export_manifest(export_manifest)

  # Each run directory is self-contained for dPCP: the sample table references
  # only amplitude files in the same directory, and the companion manifest lets
  # us audit how each dPCP well maps back to the raw archive and reviewed manifest.
  readr::write_csv(sample_table, file.path(out_dir, "sample_table.csv"))
  readr::write_csv(export_manifest, file.path(out_dir, "dpcp_input_manifest.csv"))
  readr::write_csv(reference_candidates, file.path(out_dir, "reference_candidates.csv"))
  readr::write_csv(validation, file.path(out_dir, "dpcp_input_validation.csv"))

  list(
    export_manifest = export_manifest,
    reference_candidates = reference_candidates,
    validation = validation,
    sample_table = sample_table
  )
}

process_experimental_dataset <- function() {
  runs_path <- file.path(experimental_raw_root, "manifests", "runs.csv")
  sample_manifest_path <- file.path(experimental_raw_root, "manifests", "sample_manifest.csv")

  if (!file.exists(runs_path) || !file.exists(sample_manifest_path)) {
    stop("Missing experimental ddPCR manifests under ", experimental_raw_root, call. = FALSE)
  }

  runs <- readr::read_csv(runs_path, show_col_types = FALSE) %>%
    filter(status == "active", experiment == "SNV") %>%
    mutate(
      assay = normalise_dpcp_assay(assay),
      run_date = as.Date(run_date),
      archive_dir = file.path(experimental_raw_root, archive_contents_relative_dir),
      peak_data_dir = file.path(archive_dir, "PeakData"),
      peak_data_files = vapply(
        peak_data_dir,
        function(path) {
          if (!dir.exists(path)) {
            return(0L)
          }
          length(list.files(path, pattern = "\\.ddpeakjson$", full.names = TRUE))
        },
        integer(1)
      )
    ) %>%
    arrange(run_date, assay, run_id)

  skipped <- runs %>%
    filter(peak_data_files == 0L) %>%
    transmute(
      dataset = "experimental",
      assay,
      run_id,
      skip_reason = "missing_peak_data",
      notes
    )
  runs_to_export <- runs %>%
    filter(peak_data_files > 0L)

  sample_manifest <- readr::read_csv(sample_manifest_path, show_col_types = FALSE) %>%
    filter(experiment == "SNV", normalise_dpcp_assay(assay) %in% c("D178N", "E200K", "P102L")) %>%
    mutate(
      assay = normalise_dpcp_assay(assay),
      target_clean = clean_dpcp_target(target_clean)
    )

  output_root <- file.path(input_root, "experimental")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  message(
    "Exporting experimental dPCP inputs for ",
    nrow(runs_to_export),
    " active SNV runs; skipping ",
    nrow(skipped),
    " run(s) without PeakData"
  )
  results <- purrr::map(seq_len(nrow(runs_to_export)), function(i) {
    process_experimental_run(runs_to_export[i, ], sample_manifest, output_root)
  })

  list(
    export_manifest = bind_rows(purrr::map(results, "export_manifest")),
    reference_candidates = bind_rows(purrr::map(results, "reference_candidates")),
    validation = bind_rows(purrr::map(results, "validation")),
    skipped = skipped
  )
}

# The LoD dataset stores merged pseudo-wells such as M05 in its manifest. dPCP
# can analyse those rows only if we physically concatenate the accepted droplets
# from each listed source well into one amplitude CSV.
process_lod_run <- function(assay_dir, assay, run_id, output_root) {
  archive_dir <- file.path(assay_dir, "archive_contents")
  sample_manifest_path <- file.path(assay_dir, "manifests", "sample_manifest.csv")
  if (!file.exists(sample_manifest_path)) {
    stop("Missing LoD sample manifest: ", sample_manifest_path, call. = FALSE)
  }

  sample_manifest <- readr::read_csv(sample_manifest_path, show_col_types = FALSE) %>%
    mutate(
      assay = normalise_dpcp_assay(assay),
      target_normalised = clean_dpcp_target(target_normalised)
    )

  out_dir <- file.path(output_root, tolower(assay), run_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  merged_rows <- sample_manifest %>%
    distinct(well, sample, sample_normalised, accepted_droplets, merged_wells) %>%
    arrange(well)

  export_manifest <- purrr::map_dfr(seq_len(nrow(merged_rows)), function(i) {
    merged_row <- merged_rows[i, ]
    physical_wells <- split_merged_wells(merged_row$merged_wells[[1]])
    if (length(physical_wells) == 0L) {
      stop("LoD merged well has no physical source wells: ", run_id, " ", merged_row$well[[1]], call. = FALSE)
    }

    well_record <- combine_dpcp_physical_wells(
      archive_dir = archive_dir,
      wells = physical_wells,
      assay = assay
    )

    # The LoD manifest does not carry explicit dye/channel columns, so validate
    # that the two targets present in the manifest match the archive-derived
    # Ch1/Ch2 target pair before writing the sample table.
    manifest_targets <- sample_manifest %>%
      filter(.data$well == !!merged_row$well[[1]]) %>%
      pull(target_normalised) %>%
      unique() %>%
      sort()
    archive_targets <- sort(c(well_record$ch1_target, well_record$ch2_target))
    if (!identical(manifest_targets, archive_targets)) {
      stop("LoD manifest targets disagree with archive metadata for ", run_id, " ", merged_row$well[[1]], call. = FALSE)
    }

    write_export_record(
      dataset = "lod",
      assay = assay,
      run_id = run_id,
      sample = merged_row$sample[[1]],
      source_wells = physical_wells,
      expected_accepted = as.integer(merged_row$accepted_droplets[[1]]),
      out_dir = out_dir,
      well_record = well_record,
      dpcp_well_id = make_dpcp_id(run_id, merged_row$well[[1]])
    )
  })

  reference_candidates <- score_reference_candidates(
    export_manifest,
    min_partition_count = min_reference_cluster
  )
  reference_filename <- choose_reference_filename(reference_candidates)
  sample_table <- make_dpcp_sample_table(export_manifest, reference_filename)
  validation <- validate_export_manifest(export_manifest)

  readr::write_csv(sample_table, file.path(out_dir, "sample_table.csv"))
  readr::write_csv(export_manifest, file.path(out_dir, "dpcp_input_manifest.csv"))
  readr::write_csv(reference_candidates, file.path(out_dir, "reference_candidates.csv"))
  readr::write_csv(validation, file.path(out_dir, "dpcp_input_validation.csv"))

  list(
    export_manifest = export_manifest,
    reference_candidates = reference_candidates,
    validation = validation,
    sample_table = sample_table
  )
}

process_lod_dataset <- function() {
  index_path <- file.path(lod_raw_root, "manifest_index.csv")
  if (!file.exists(index_path)) {
    stop("Missing LoD manifest index: ", index_path, call. = FALSE)
  }

  manifest_index <- readr::read_csv(index_path, show_col_types = FALSE) %>%
    mutate(
      assay = normalise_dpcp_assay(assay),
      assay_dir = file.path(lod_raw_root, tolower(assay))
    ) %>%
    arrange(assay, run_id)

  output_root <- file.path(input_root, "lod")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  message("Exporting LoD dPCP inputs for ", nrow(manifest_index), " assay runs")
  results <- purrr::map(seq_len(nrow(manifest_index)), function(i) {
    row <- manifest_index[i, ]
    process_lod_run(
      assay_dir = row$assay_dir[[1]],
      assay = row$assay[[1]],
      run_id = row$run_id[[1]],
      output_root = output_root
    )
  })

  list(
    export_manifest = bind_rows(purrr::map(results, "export_manifest")),
    reference_candidates = bind_rows(purrr::map(results, "reference_candidates")),
    validation = bind_rows(purrr::map(results, "validation")),
    skipped = tibble(
      dataset = character(),
      assay = character(),
      run_id = character(),
      skip_reason = character(),
      notes = character()
    )
  )
}

all_results <- list()
if (dataset %in% c("all", "experimental")) {
  all_results$experimental <- process_experimental_dataset()
}
if (dataset %in% c("all", "lod")) {
  all_results$lod <- process_lod_dataset()
}

combined_manifest <- bind_rows(purrr::map(all_results, "export_manifest"))
combined_reference_candidates <- bind_rows(purrr::map(all_results, "reference_candidates"))
combined_validation <- bind_rows(purrr::map(all_results, "validation"))
combined_skipped <- bind_rows(purrr::map(all_results, "skipped"))

# The combined validation files are not consumed by dPCP directly. They are the
# repository-level audit trail proving that every generated amplitude file exists,
# has the expected accepted-droplet count, and carries a non-degenerate Ch1/Ch2
# target mapping before downstream dPCP analysis begins.
readr::write_csv(combined_manifest, file.path(validation_root, "dpcp_input_manifest.csv"))
readr::write_csv(combined_reference_candidates, file.path(validation_root, "reference_candidates.csv"))
readr::write_csv(combined_validation, file.path(validation_root, "dpcp_input_validation.csv"))
readr::write_csv(combined_skipped, file.path(validation_root, "skipped_runs.csv"))

if (any(!combined_validation$valid)) {
  bad_path <- file.path(validation_root, "dpcp_input_validation.csv")
  stop("dPCP input validation failed. See ", bad_path, call. = FALSE)
}

missing_reference_runs <- combined_reference_candidates %>%
  group_by(dataset, assay, run_id) %>%
  summarise(has_reference = any(reference_candidate), .groups = "drop") %>%
  filter(!has_reference)

if (nrow(missing_reference_runs) > 0L) {
  readr::write_csv(missing_reference_runs, file.path(validation_root, "runs_without_reference_candidates.csv"))
  warning(
    nrow(missing_reference_runs),
    " run(s) have no automatic reference candidate. Sample tables were still written with blank Reference fields."
  )
}

message("Wrote dPCP input manifest: ", file.path(validation_root, "dpcp_input_manifest.csv"))
message("Wrote dPCP input validation: ", file.path(validation_root, "dpcp_input_validation.csv"))
