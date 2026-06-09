library(dplyr)
library(jsonlite)
library(purrr)
library(readr)
library(stringr)
library(tibble)

# Umbrella uses the same raw ddPCR source material as the dPCP compatibility
# export, so reuse the low-level metadata parsers and naming helpers rather
# than maintaining a parallel vocabulary for assays, targets, channels, and
# repository-relative paths.
dpcp_helper_path <- file.path(getwd(), "src", "ddpcr", "dpcp", "dpcp_input_helpers.R")
if (!file.exists(dpcp_helper_path)) {
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) character(0))
  if (length(source_file) == 1L && !is.null(source_file) && nzchar(source_file)) {
    dpcp_helper_path <- file.path(dirname(dirname(normalizePath(source_file, winslash = "/", mustWork = TRUE))), "dpcp", "dpcp_input_helpers.R")
  }
}
if (!file.exists(dpcp_helper_path)) {
  stop("Could not locate dPCP helper dependency for Umbrella input conversion.", call. = FALSE)
}
source(dpcp_helper_path)

# Keep the Umbrella identifiers visually distinct in manifests while using the
# exact same sanitisation rules as dPCP. This makes cross-method comparisons
# straightforward: the same run/well pair maps to the same stable base ID.
make_umbrella_id <- make_dpcp_id

# Umbrella can use a prior four-level cluster column when it is present. The
# numeric labels here intentionally encode the physical channel state:
#   1 = Ch1-/Ch2-
#   2 = Ch1+/Ch2-
#   3 = Ch1-/Ch2+
#   4 = Ch1+/Ch2+
# Umbrella re-identifies clusters by median fluorescence, so the labels are not
# interpreted biologically by name, but this ordering is still useful for audit.
umbrella_cluster_code <- function(ch1, ch2) {
  if (ch1 == "Negative" && ch2 == "Negative") {
    return(1L)
  }
  if (ch1 == "Positive" && ch2 == "Negative") {
    return(2L)
  }
  if (ch1 == "Negative" && ch2 == "Positive") {
    return(3L)
  }
  if (ch1 == "Positive" && ch2 == "Positive") {
    return(4L)
  }
  NA_integer_
}

# Read one physical well and return a data frame in Umbrella's expected shape:
# the first two columns are fluorescence intensities and the optional Cluster
# column carries the accepted-droplet cluster call reconstructed from archive
# metadata. Gated/unassigned droplets are deliberately excluded, matching the
# manifest count contract used by the official ddPCR workflow.
read_umbrella_physical_well <- function(archive_dir, well, assay) {
  assay <- normalise_dpcp_assay(assay)
  metadata_path <- file.path(archive_dir, "PeakMetaData", paste0(well, ".ddmetajson"))
  peak_path <- file.path(archive_dir, "PeakData", paste0(well, ".ddpeakjson"))

  if (!file.exists(metadata_path)) {
    stop("Missing peak metadata: ", metadata_path, call. = FALSE)
  }
  if (!file.exists(peak_path)) {
    stop("Missing peak data: ", peak_path, call. = FALSE)
  }

  metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)
  peak <- jsonlite::fromJSON(peak_path, simplifyVector = FALSE)

  targets <- metadata_targets(metadata)
  if (is.null(targets)) {
    stop("Could not read target metadata for ", archive_dir, " ", well, call. = FALSE)
  }

  selected <- selected_target_indices(targets, assay)
  if (is.null(selected)) {
    stop("Could not select WT and mutant targets for ", archive_dir, " ", well, call. = FALSE)
  }

  target_names <- vapply(targets, target_name, character(1))
  target_clean <- clean_dpcp_target(target_names)
  channels <- vapply(targets, target_channel, integer(1))
  selected_indices <- c(selected$ref, selected$mut)
  ch1_idx <- selected_indices[channels[selected_indices] == 1L][1]
  ch2_idx <- selected_indices[channels[selected_indices] == 2L][1]

  if (is.na(ch1_idx) || is.na(ch2_idx)) {
    stop("Could not map selected targets to Ch1 and Ch2 for ", archive_dir, " ", well, call. = FALSE)
  }

  partition_counts <- c(
    `Ch1+Ch2+` = 0L,
    `Ch1+Ch2-` = 0L,
    `Ch1-Ch2+` = 0L,
    `Ch1-Ch2-` = 0L
  )
  accepted_indices <- integer(0)
  accepted_clusters <- integer(0)
  gated_or_unassigned <- 0L

  for (cluster in metadata$Clusters %||% list()) {
    droplets <- as.integer(cluster$Droplets %||% integer(0))
    if (length(droplets) == 0L) {
      next
    }

    results <- as.character(unlist(cluster$Results, use.names = FALSE))
    if (length(results) < max(selected_indices) || is_json_true(cluster$Unassigned)) {
      gated_or_unassigned <- gated_or_unassigned + length(droplets)
      next
    }

    ch1 <- results[[ch1_idx]]
    ch2 <- results[[ch2_idx]]
    if (!all(c(ch1, ch2) %in% c("Negative", "Positive"))) {
      gated_or_unassigned <- gated_or_unassigned + length(droplets)
      next
    }

    key <- paste0(
      ifelse(ch1 == "Positive", "Ch1+", "Ch1-"),
      ifelse(ch2 == "Positive", "Ch2+", "Ch2-")
    )
    cluster_code <- umbrella_cluster_code(ch1, ch2)

    partition_counts[[key]] <- partition_counts[[key]] + length(droplets)
    accepted_indices <- c(accepted_indices, droplets)
    accepted_clusters <- c(accepted_clusters, rep(cluster_code, length(droplets)))
  }

  if (length(accepted_indices) == 0L) {
    stop("No accepted droplets for ", archive_dir, " ", well, call. = FALSE)
  }
  if (anyDuplicated(accepted_indices)) {
    stop("Duplicate accepted droplet indices for ", archive_dir, " ", well, call. = FALSE)
  }

  amplitudes <- peak$PeakInfo$Amplitudes
  if (is.null(amplitudes) || length(amplitudes) < 2L) {
    stop("Peak data does not contain two amplitude channels: ", peak_path, call. = FALSE)
  }

  ch1_amplitudes <- as.numeric(amplitudes[[1]])
  ch2_amplitudes <- as.numeric(amplitudes[[2]])
  if (length(ch1_amplitudes) != length(ch2_amplitudes)) {
    stop("Amplitude channel lengths differ in ", peak_path, call. = FALSE)
  }

  # Metadata droplet identifiers are zero-based indices into the PeakData
  # arrays. Sort once and carry the cluster labels through the same ordering so
  # the optional Umbrella prior clusters remain aligned to the amplitudes.
  order_idx <- order(accepted_indices)
  accepted_indices <- accepted_indices[order_idx]
  accepted_clusters <- accepted_clusters[order_idx]
  r_indices <- accepted_indices + 1L
  if (max(r_indices) > length(ch1_amplitudes)) {
    stop("Accepted droplet index exceeds amplitude array length in ", peak_path, call. = FALSE)
  }

  data <- tibble(
    `Ch1.Amplitude` = ch1_amplitudes[r_indices],
    `Ch2.Amplitude` = ch2_amplitudes[r_indices],
    Cluster = accepted_clusters
  )

  list(
    well = well,
    data = data,
    accepted_count = nrow(data),
    partition_counts = partition_counts,
    gated_or_unassigned = gated_or_unassigned,
    ch1_target = target_clean[[ch1_idx]],
    ch2_target = target_clean[[ch2_idx]],
    ch1_raw_target = target_names[[ch1_idx]],
    ch2_raw_target = target_names[[ch2_idx]],
    metadata_path = metadata_path,
    peak_path = peak_path
  )
}

# Merge LoD physical wells into a single Umbrella partition set. The merge is
# row-binding accepted droplets only; no averaging or summary-count conversion
# is performed, because Umbrella models the droplet-level fluorescence data.
combine_umbrella_physical_wells <- function(archive_dir, wells, assay) {
  records <- purrr::map(wells, read_umbrella_physical_well, archive_dir = archive_dir, assay = assay)
  first <- records[[1]]

  geometry <- purrr::map_chr(records, ~ paste(.x$ch1_target, .x$ch2_target, sep = "|"))
  if (length(unique(geometry)) != 1L) {
    stop("Merged wells do not share the same target/channel geometry: ", paste(wells, collapse = ", "), call. = FALSE)
  }

  list(
    wells = wells,
    data = bind_rows(purrr::map(records, "data")),
    accepted_count = sum(purrr::map_int(records, "accepted_count")),
    partition_counts = Reduce(`+`, purrr::map(records, "partition_counts")),
    gated_or_unassigned = sum(purrr::map_int(records, "gated_or_unassigned")),
    ch1_target = first$ch1_target,
    ch2_target = first$ch2_target,
    ch1_raw_target = first$ch1_raw_target,
    ch2_raw_target = first$ch2_raw_target,
    source_peak_files = paste(vapply(records, `[[`, character(1), "peak_path"), collapse = "|"),
    source_metadata_files = paste(vapply(records, `[[`, character(1), "metadata_path"), collapse = "|")
  )
}

write_umbrella_partition_csv <- function(data, path) {
  expected <- c("Ch1.Amplitude", "Ch2.Amplitude", "Cluster")
  if (!identical(names(data), expected)) {
    stop("Umbrella data must have Ch1.Amplitude, Ch2.Amplitude, and Cluster columns.", call. = FALSE)
  }
  if (!all(vapply(data[1:2], is.numeric, logical(1)))) {
    stop("Umbrella amplitude columns must be numeric.", call. = FALSE)
  }
  if (!is.integer(data$Cluster) && !is.numeric(data$Cluster)) {
    stop("Umbrella Cluster column must be numeric.", call. = FALSE)
  }
  readr::write_csv(data, path)
}

is_ntc_sample <- function(sample) {
  str_detect(sample %||% "", regex("NTC", ignore_case = TRUE))
}

validate_umbrella_manifest <- function(export_manifest) {
  export_manifest %>%
    mutate(
      missing_partition_file = !file.exists(partition_set_path),
      count_matches_manifest = accepted_droplets_exported == accepted_droplets_manifest,
      has_two_targets = !is.na(ch1_target) & !is.na(ch2_target) & ch1_target != ch2_target,
      cluster_levels_ok = str_detect(cluster_levels, "^[1-4](\\|[1-4])*$"),
      list_name_present = nzchar(umbrella_name),
      valid = !missing_partition_file &
        count_matches_manifest &
        has_two_targets &
        cluster_levels_ok &
        list_name_present
    )
}
