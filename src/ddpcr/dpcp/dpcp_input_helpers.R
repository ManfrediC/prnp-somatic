library(dplyr)
library(jsonlite)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

# Return the fallback value when an optional JSON field is absent.
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Keep the three SNV assay labels stable across source filenames, manifests,
# and older ddPCR metadata where the D178N assay was occasionally mistyped.
normalise_dpcp_assay <- function(x) {
  x <- toupper(as.character(x))
  x <- str_replace_all(x, "D1789N", "D178N")
  str_extract(x, "D178N|E200K|P102L")
}

# Match the target clean-up rules used by the main ddPCR pipeline so the dPCP
# compatibility export does not introduce a second vocabulary for the same
# biological targets.
clean_dpcp_target <- function(x) {
  out <- as.character(x)
  out <- str_replace_all(out, "-mut|_FAM1|_VIC2", "")
  out <- str_replace_all(out, "D1789N", "D178N")
  out[out == "PRNP"] <- "WT"
  out
}

# Convert run/sample/well identifiers into filesystem-safe dPCP well IDs.
# The date prefix is compacted to YYYYMMDD so generated filenames remain short
# but still preserve the run identity needed to avoid repeated A01/B01 names.
make_dpcp_id <- function(...) {
  raw <- paste(..., sep = "_")
  raw <- str_replace(raw, "^(\\d{4})-(\\d{2})-(\\d{2})", "\\1\\2\\3")
  raw <- str_replace_all(raw, "[^A-Za-z0-9]+", "_")
  raw <- str_replace_all(raw, "_+", "_")
  str_replace_all(raw, "^_|_$", "")
}

# Recover the repository root from either Rscript execution, source() calls,
# RStudio, or a working directory somewhere inside the checkout. This mirrors
# the path strategy used by the existing ddPCR scripts so the converter can be
# launched from the repository root or from an editor.
get_dpcp_project_root <- function() {
  normalise_existing_dir <- function(path) {
    if (length(path) != 1L || is.na(path) || !nzchar(path)) {
      return(character(0))
    }
    if (file.exists(path) && !dir.exists(path)) {
      path <- dirname(path)
    }
    if (!dir.exists(path)) {
      return(character(0))
    }
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }

  find_project_root <- function(start_dir) {
    current_dir <- normalise_existing_dir(start_dir)
    if (length(current_dir) == 0L) {
      return(character(0))
    }

    repeat {
      if (dir.exists(file.path(current_dir, "src", "ddpcr")) &&
          file.exists(file.path(current_dir, "env", "ddpcr.environment.yml"))) {
        return(current_dir)
      }
      parent_dir <- dirname(current_dir)
      if (identical(parent_dir, current_dir)) {
        return(character(0))
      }
      current_dir <- parent_dir
    }
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  source_files <- unlist(lapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) character(0) else frame$ofile
  }))

  rstudio_path <- character(0)
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    rstudio_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) character(0)
    )
  }

  candidate_dirs <- unique(c(
    normalise_existing_dir(sub("^--file=", "", file_arg[1])),
    unlist(lapply(source_files, normalise_existing_dir), use.names = FALSE),
    normalise_existing_dir(rstudio_path[1]),
    normalise_existing_dir(getwd())
  ))

  for (candidate_dir in candidate_dirs) {
    project_root <- find_project_root(candidate_dir)
    if (length(project_root) == 1L) {
      return(project_root)
    }
  }

  stop(
    "Could not determine project root. ",
    "Run from inside the prnp-somatic repository or source the script file directly.",
    call. = FALSE
  )
}

# Metadata values can be stored as logical TRUE or as string values depending
# on QX Manager/QuantaSoft export vintage.
is_json_true <- function(x) {
  isTRUE(x) || identical(tolower(as.character(x)), "true")
}

# Extract one scalar metadata value without assuming every archive uses the
# same exact target object layout.
target_value <- function(target, key, default = NA_character_) {
  value <- target[[key]]
  if (is.null(value) || length(value) == 0L) {
    return(default)
  }
  as.character(value[[1]])
}

# The explicit TargetName field is preferred, but some archive metadata stores
# the same label under Name.
target_name <- function(target) {
  value <- target_value(target, "TargetName")
  if (is.na(value)) {
    value <- target_value(target, "Name")
  }
  value
}

# Pull the ddPCR fluorescence channel from either the newer single Dye object
# or the older Dyes list representation.
target_channel <- function(target) {
  dye <- target[["Dye"]]
  if (!is.null(dye) && !is.null(dye[["Channel"]])) {
    return(as.integer(dye[["Channel"]]))
  }

  dyes <- target[["Dyes"]]
  if (!is.null(dyes)) {
    channels <- purrr::map_int(dyes, function(dye_entry) {
      if (is.null(dye_entry) || is.null(dye_entry[["Channel"]])) {
        return(NA_integer_)
      }
      as.integer(dye_entry[["Channel"]])
    })
    channels <- channels[!is.na(channels)]
    if (length(channels) > 0L) {
      return(channels[[1]])
    }
  }

  NA_integer_
}

# The target list is stored inside one of the metadata clusters. We keep the
# first cluster with at least two target definitions, which is the structure
# used by the extracted QX Manager archives in this repository.
metadata_targets <- function(metadata) {
  clusters <- metadata$Clusters %||% list()
  for (cluster in clusters) {
    targets <- cluster$Targets
    if (!is.null(targets) && length(targets) >= 2L) {
      return(targets)
    }
  }
  NULL
}

# Select the WT and assay-specific mutant target entries, requiring them to sit
# on opposite channels. dPCP needs channel amplitudes, while the source archive
# stores target calls by target index, so this target-to-channel bridge is the
# critical contract.
selected_target_indices <- function(targets, assay) {
  names <- vapply(targets, target_name, character(1))
  cleaned <- clean_dpcp_target(names)
  channels <- vapply(targets, target_channel, integer(1))

  mut_candidates <- which(cleaned == assay)
  ref_candidates <- which(cleaned == "WT")
  if (length(mut_candidates) == 0L || length(ref_candidates) == 0L) {
    return(NULL)
  }

  for (mut_idx in mut_candidates) {
    ref_idx <- ref_candidates[channels[ref_candidates] != channels[mut_idx]][1]
    if (!is.na(ref_idx)) {
      return(list(ref = ref_idx, mut = mut_idx))
    }
  }

  NULL
}

# Read one physical well from extracted archive contents and return the
# droplet-level amplitudes that dPCP expects. Accepted droplets are reconstructed
# from the metadata cluster assignments, not from the summary CSV, so the script
# can later prove that the row counts agree with the manifest.
read_dpcp_physical_well <- function(archive_dir, well, assay) {
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

  # dPCP interprets the raw amplitude CSV positionally: the first column is the
  # y-axis channel and the second column is the x-axis channel. Preserve the
  # Bio-Rad Channel 1/Channel 2 order even when WT and mutant swap dyes between
  # assays, then put the matching target names into the sample table later.
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
    partition_counts[[key]] <- partition_counts[[key]] + length(droplets)
    accepted_indices <- c(accepted_indices, droplets)
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

  # The archive stores droplet indices as zero-based offsets into the amplitude
  # arrays. Convert once to R's one-based indexing and keep the original droplet
  # order stable for reproducibility.
  accepted_indices <- sort(accepted_indices)
  r_indices <- accepted_indices + 1L
  if (max(r_indices) > length(ch1_amplitudes)) {
    stop("Accepted droplet index exceeds amplitude array length in ", peak_path, call. = FALSE)
  }

  data <- tibble(
    `Ch1.Amplitude` = ch1_amplitudes[r_indices],
    `Ch2.Amplitude` = ch2_amplitudes[r_indices]
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

# Bind multiple physical wells into one dPCP pseudo-well. This is used only for
# the LoD material, whose analysis CSV and manifest define merged wells such as
# M05 from several physical wells.
combine_dpcp_physical_wells <- function(archive_dir, wells, assay) {
  records <- purrr::map(wells, read_dpcp_physical_well, archive_dir = archive_dir, assay = assay)
  first <- records[[1]]

  # All wells in a merged LoD pseudo-well must describe the same assay geometry.
  # If this fails, merging would create a biologically incoherent dPCP file.
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

# dPCP requires this exact nine-column table shape. Column names are intentionally
# human-readable, matching the package template rather than tidyverse style.
dpcp_sample_table_columns <- c(
  "Sample name",
  "Chip ID/Well ID",
  "No of targets",
  "FAM target",
  "Target 3",
  "Target 4",
  "VIC/HEX target",
  "Reference",
  "Dilution"
)

make_dpcp_sample_table <- function(export_manifest, reference_filename = "") {
  out <- export_manifest %>%
    transmute(
      `Sample name` = sample,
      `Chip ID/Well ID` = dpcp_well_id,
      `No of targets` = 2L,
      `FAM target` = ch1_target,
      `Target 3` = "",
      `Target 4` = "",
      `VIC/HEX target` = ch2_target,
      Reference = reference_filename,
      Dilution = 1
    )

  out[dpcp_sample_table_columns]
}

# Pick a pragmatic DBSCAN reference candidate. The selected reference must have
# enough empty droplets and both single-positive clusters; otherwise dPCP cannot
# learn the geometry of the assay from that well alone.
score_reference_candidates <- function(export_manifest, min_partition_count = 50L) {
  export_manifest %>%
    mutate(
      min_required_cluster = pmin(ch1_only_droplets, ch2_only_droplets, empty_droplets),
      reference_candidate = min_required_cluster >= min_partition_count,
      reference_score = min_required_cluster + 0.001 * accepted_droplets_exported
    ) %>%
    arrange(desc(reference_candidate), desc(reference_score), dpcp_well_id)
}

choose_reference_filename <- function(reference_candidates) {
  selected <- reference_candidates %>%
    filter(reference_candidate) %>%
    slice_head(n = 1)

  if (nrow(selected) == 0L) {
    return("")
  }

  selected$amplitude_file[[1]]
}

split_merged_wells <- function(x) {
  str_split(x, ",")[[1]] %>%
    str_trim() %>%
    discard(~ !nzchar(.x))
}

write_dpcp_amplitude_csv <- function(data, path) {
  # Keep this guard close to the write call so future edits cannot accidentally
  # pass labels, cluster calls, or extra metadata columns into the dPCP raw-data
  # files. dPCP expects just the two fluorescence amplitude columns.
  if (!identical(names(data), c("Ch1.Amplitude", "Ch2.Amplitude"))) {
    stop("Amplitude data must have exactly Ch1.Amplitude and Ch2.Amplitude columns.", call. = FALSE)
  }
  if (!all(vapply(data, is.numeric, logical(1)))) {
    stop("Amplitude columns must be numeric.", call. = FALSE)
  }
  readr::write_csv(data, path)
}

validate_export_manifest <- function(export_manifest) {
  # This is a structural input check, not a clustering QC result. It verifies that
  # the conversion output is internally consistent and dPCP-shaped; DBSCAN/c-means
  # fit quality still needs to be reviewed in dPCP when the package is run.
  export_manifest %>%
    mutate(
      missing_amplitude_file = !file.exists(amplitude_path),
      count_matches_manifest = accepted_droplets_exported == accepted_droplets_manifest,
      has_two_targets = !is.na(ch1_target) & !is.na(ch2_target) & ch1_target != ch2_target,
      valid = !missing_amplitude_file & count_matches_manifest & has_two_targets
    )
}
