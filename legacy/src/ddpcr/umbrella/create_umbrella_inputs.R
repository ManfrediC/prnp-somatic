#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

# Umbrella is distributed as R scripts rather than an installed package. The
# converter therefore creates self-contained RData files with the exact data
# objects Umbrella expects, plus CSV mirrors and manifests so the transformation
# remains auditable without loading R workspaces.
command_line <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_line, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  file.path(getwd(), "src", "ddpcr", "umbrella")
}
source(file.path(script_dir, "umbrella_input_helpers.R"))
project_root <- get_dpcp_project_root()

args <- commandArgs(trailingOnly = TRUE)
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

experimental_raw_root <- file.path(project_root, "raw", "ddpcr")
lod_raw_root <- file.path(project_root, "raw", "ddpcr_lod")
input_root <- file.path(project_root, "results", "ddPCR", "umbrella", "inputs")
validation_root <- file.path(project_root, "results", "ddPCR", "umbrella", "validation")
analysis_root <- file.path(project_root, "results", "ddPCR", "umbrella", "analysis")
vendor_root <- file.path(project_root, "src", "ddpcr", "umbrella", "vendor", "statOmics_umbrella")

dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
dir.create(validation_root, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)

# Write one partition-set CSV and return the manifest row that ties Umbrella's
# list entry back to the source run, source well(s), archive JSON files, target
# mapping, and manifest droplet counts.
write_umbrella_export_record <- function(
  dataset,
  assay,
  run_id,
  sample,
  source_wells,
  expected_accepted,
  out_dir,
  well_record,
  umbrella_name
) {
  partition_dir <- file.path(out_dir, "partition_sets")
  dir.create(partition_dir, recursive = TRUE, showWarnings = FALSE)

  partition_set_file <- paste0(umbrella_name, "_Amplitude.csv")
  partition_set_path <- file.path(partition_dir, partition_set_file)
  write_umbrella_partition_csv(well_record$data, partition_set_path)

  tibble(
    dataset = dataset,
    assay = assay,
    run_id = run_id,
    sample = sample,
    is_ntc = is_ntc_sample(sample),
    source_wells = paste(source_wells, collapse = "|"),
    umbrella_name = umbrella_name,
    partition_set_file = partition_set_file,
    partition_set_path = partition_set_path,
    accepted_droplets_exported = nrow(well_record$data),
    accepted_droplets_manifest = expected_accepted,
    double_positive_droplets = as.integer(well_record$partition_counts[["Ch1+Ch2+"]]),
    ch1_only_droplets = as.integer(well_record$partition_counts[["Ch1+Ch2-"]]),
    ch2_only_droplets = as.integer(well_record$partition_counts[["Ch1-Ch2+"]]),
    empty_droplets = as.integer(well_record$partition_counts[["Ch1-Ch2-"]]),
    gated_or_unassigned_droplets = well_record$gated_or_unassigned,
    cluster_levels = paste(sort(unique(well_record$data$Cluster)), collapse = "|"),
    ch1_target = well_record$ch1_target,
    ch2_target = well_record$ch2_target,
    ch1_raw_target = well_record$ch1_raw_target,
    ch2_raw_target = well_record$ch2_raw_target,
    source_peak_files = well_record$source_peak_files %||% well_record$peak_path,
    source_metadata_files = well_record$source_metadata_files %||% well_record$metadata_path
  )
}

save_umbrella_run_objects <- function(out_dir, export_manifest, plate_list) {
  umbrella_plate_2d <- plate_list
  umbrella_ntc_names <- names(umbrella_plate_2d)[export_manifest$is_ntc]
  umbrella_metadata <- export_manifest %>%
    select(
      dataset, assay, run_id, sample, is_ntc, source_wells, umbrella_name,
      accepted_droplets_exported, accepted_droplets_manifest,
      ch1_target, ch2_target, partition_set_file
    )

  save(
    umbrella_plate_2d,
    umbrella_metadata,
    umbrella_ntc_names,
    file = file.path(out_dir, "umbrella_input.RData")
  )

  writeLines(umbrella_ntc_names, file.path(out_dir, "ntc_names.txt"))
  invisible(umbrella_ntc_names)
}

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

  plate_list <- list()
  export_manifest <- purrr::map_dfr(seq_len(nrow(well_rows)), function(i) {
    well_row <- well_rows[i, ]
    well <- well_row$well[[1]]
    umbrella_name <- make_umbrella_id(run_id, well)

    well_record <- read_umbrella_physical_well(
      archive_dir = archive_dir,
      well = well,
      assay = assay
    )

    # Keep the manifest channel assignment and archive target metadata in lock
    # step. If these disagree, the RData object would still run but the channel
    # labels in downstream result tables would be untrustworthy.
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

    plate_list[[umbrella_name]] <<- well_record$data
    write_umbrella_export_record(
      dataset = "experimental",
      assay = assay,
      run_id = run_id,
      sample = well_row$sample[[1]],
      source_wells = well,
      expected_accepted = as.integer(well_row$accepted_droplets[[1]]),
      out_dir = out_dir,
      well_record = well_record,
      umbrella_name = umbrella_name
    )
  })

  ntc_names <- save_umbrella_run_objects(out_dir, export_manifest, plate_list)
  validation <- validate_umbrella_manifest(export_manifest)

  readr::write_csv(export_manifest, file.path(out_dir, "umbrella_input_manifest.csv"))
  readr::write_csv(validation, file.path(out_dir, "umbrella_input_validation.csv"))

  list(
    export_manifest = export_manifest,
    validation = validation,
    ntc_selection = tibble(
      dataset = "experimental",
      assay = assay,
      run_id = run_id,
      ntc_count = length(ntc_names),
      ntc_names = paste(ntc_names, collapse = "|")
    )
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
    "Exporting experimental Umbrella inputs for ",
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
    validation = bind_rows(purrr::map(results, "validation")),
    ntc_selection = bind_rows(purrr::map(results, "ntc_selection")),
    skipped = skipped
  )
}

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

  plate_list <- list()
  export_manifest <- purrr::map_dfr(seq_len(nrow(merged_rows)), function(i) {
    merged_row <- merged_rows[i, ]
    physical_wells <- split_merged_wells(merged_row$merged_wells[[1]])
    if (length(physical_wells) == 0L) {
      stop("LoD merged well has no physical source wells: ", run_id, " ", merged_row$well[[1]], call. = FALSE)
    }

    umbrella_name <- make_umbrella_id(run_id, merged_row$well[[1]])
    well_record <- combine_umbrella_physical_wells(
      archive_dir = archive_dir,
      wells = physical_wells,
      assay = assay
    )

    manifest_targets <- sample_manifest %>%
      filter(.data$well == !!merged_row$well[[1]]) %>%
      pull(target_normalised) %>%
      unique() %>%
      sort()
    archive_targets <- sort(c(well_record$ch1_target, well_record$ch2_target))
    if (!identical(manifest_targets, archive_targets)) {
      stop("LoD manifest targets disagree with archive metadata for ", run_id, " ", merged_row$well[[1]], call. = FALSE)
    }

    plate_list[[umbrella_name]] <<- well_record$data
    write_umbrella_export_record(
      dataset = "lod",
      assay = assay,
      run_id = run_id,
      sample = merged_row$sample[[1]],
      source_wells = physical_wells,
      expected_accepted = as.integer(merged_row$accepted_droplets[[1]]),
      out_dir = out_dir,
      well_record = well_record,
      umbrella_name = umbrella_name
    )
  })

  ntc_names <- save_umbrella_run_objects(out_dir, export_manifest, plate_list)
  validation <- validate_umbrella_manifest(export_manifest)

  readr::write_csv(export_manifest, file.path(out_dir, "umbrella_input_manifest.csv"))
  readr::write_csv(validation, file.path(out_dir, "umbrella_input_validation.csv"))

  list(
    export_manifest = export_manifest,
    validation = validation,
    ntc_selection = tibble(
      dataset = "lod",
      assay = assay,
      run_id = run_id,
      ntc_count = length(ntc_names),
      ntc_names = paste(ntc_names, collapse = "|")
    )
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

  message("Exporting LoD Umbrella inputs for ", nrow(manifest_index), " assay runs")
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
    validation = bind_rows(purrr::map(results, "validation")),
    ntc_selection = bind_rows(purrr::map(results, "ntc_selection")),
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
combined_validation <- bind_rows(purrr::map(all_results, "validation"))
combined_ntc_selection <- bind_rows(purrr::map(all_results, "ntc_selection"))
combined_skipped <- bind_rows(purrr::map(all_results, "skipped"))

readr::write_csv(combined_manifest, file.path(validation_root, "umbrella_input_manifest.csv"))
readr::write_csv(combined_validation, file.path(validation_root, "umbrella_input_validation.csv"))
readr::write_csv(combined_ntc_selection, file.path(validation_root, "ntc_selection.csv"))
readr::write_csv(combined_skipped, file.path(validation_root, "skipped_runs.csv"))

if (any(!combined_validation$valid)) {
  bad_path <- file.path(validation_root, "umbrella_input_validation.csv")
  stop("Umbrella input validation failed. See ", bad_path, call. = FALSE)
}

if (any(combined_ntc_selection$ntc_count == 0L)) {
  warning("Some Umbrella input runs have no NTC partition set. See ntc_selection.csv.")
}

message("Wrote Umbrella input manifest: ", file.path(validation_root, "umbrella_input_manifest.csv"))
message("Wrote Umbrella input validation: ", file.path(validation_root, "umbrella_input_validation.csv"))
message("Vendored Umbrella scripts: ", vendor_root)
