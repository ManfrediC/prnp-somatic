library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(jsonlite)
library(ggplot2)
library(openxlsx)
library(binom)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
experiment_root <- file.path(project_root, "scratch", "ddpcr_gating_method_comparison_v2")
raw_root <- file.path(project_root, "raw", "ddpcr")

source(file.path(project_root, "src", "ddpcr", "ddpcr_raw_import_helpers.R"))

mutation_order <- c("D178N", "E200K", "P102L")
target_class_levels <- c("NN", "WT", "MUT", "DP", "Rain", "Unclassified")
quanta_class_levels <- c(
  "Reference-only",
  "Mutant-only",
  "Double-positive",
  "Double-negative",
  "Gated/unassigned",
  "Rejected/unassigned"
)
accepted_quanta_classes <- quanta_class_levels[1:4]
lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)
fa_scale <- 100

path_v2 <- function(...) {
  file.path(experiment_root, ...)
}

setup_v2_dirs <- function() {
  dirs <- c(
    "code", "literature", "inputs/shared", "inputs/twoddpcr",
    "inputs/ddPCRclust", "inputs/dPCP", "inputs/definetherain",
    "inputs/ddpcRquant", "models/control_geometry", "models/twoddpcr",
    "models/ddPCRclust", "models/dPCP", "models/definetherain",
    "models/ddpcRquant", "models/bayesian_mixture",
    "models/polygon_gates", "data",
    "data/droplets", "tables",
    "plots/individual", "plots/panels", "report", "logs"
  )
  walk(file.path(experiment_root, dirs), dir.create, recursive = TRUE, showWarnings = FALSE)
}

safe_file_component <- function(x) {
  x <- str_squish(as.character(x))
  x <- str_replace_all(x, "[^A-Za-z0-9._-]+", "_")
  x <- str_replace_all(x, "_+", "_")
  str_replace_all(x, "^_+|_+$", "")
}

as_true_flag <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

sample_type_from_label <- function(sample) {
  sample_upper <- str_to_upper(sample)
  case_when(
    str_detect(sample_upper, "NTC") ~ "NTC",
    str_detect(sample_upper, "^WT|[_-]WT|WT[_-]") ~ "WT_control",
    str_detect(sample_upper, "MUT") ~ "positive_control",
    sample %in% c("E200K", "P102L", "D178N", "D178") ~ "positive_control",
    TRUE ~ "biological_sample"
  )
}

clean_analysis_sample <- function(sample) {
  sample <- as.character(sample)
  sample[grepl("NTC", sample, ignore.case = TRUE)] <- "NTC"
  sample[grepl("Mut", sample, ignore.case = TRUE)] <- "mutant_control"
  sample[grepl("WT", sample, ignore.case = TRUE)] <- "WT_control"
  sample[sample %in% c("E200K", "P102L", "D178N", "D178")] <- "mutant_control"

  sample %>%
    str_replace_all("CJD-|CJD|D178N|D178|E200K|P102L", "") %>%
    str_replace("_$", "") %>%
    str_replace_all("-bg", "_bg") %>%
    str_replace_all("-cb", "_cb") %>%
    str_replace_all("-fr", "_fr") %>%
    str_replace_all("-hc", "_hc") %>%
    str_replace_all("-sn", "_sn") %>%
    str_replace_all("-th", "_th") %>%
    str_replace_all("-pons", "_ps") %>%
    str_replace_all("_pons", "_ps") %>%
    str_replace_all("-cau", "_bg") %>%
    str_replace_all("-ce", "_cb") %>%
    str_replace_all("-fc", "_fr") %>%
    str_replace_all("-hip", "_hc") %>%
    str_replace_all("-mdb", "_sn")
}

normalise_sample_key <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("d1789n", "d178n") %>%
    str_replace_all("(?:[_-](?:d178n|e200k|p102l)(?:[_-]?mut)?)$", "") %>%
    str_replace_all("(?:^|[_-])(?:d178n|e200k|p102l)[_-]?mut$", "") %>%
    str_replace_all("^cjd[-_]?", "") %>%
    str_replace_all("pons", "ps") %>%
    str_replace_all("substantia[_ -]*nigra", "sn") %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

build_well_manifest <- function() {
  runs <- read_csv(file.path(raw_root, "manifests", "runs.csv"), show_col_types = FALSE) %>%
    filter(status == "active", experiment == "SNV") %>%
    select(run_id, archive_contents_relative_dir)

  read_csv(file.path(raw_root, "manifests", "sample_manifest.csv"), show_col_types = FALSE) %>%
    filter(experiment == "SNV", assay %in% mutation_order) %>%
    group_by(run_id, run_date, assay, well, sample) %>%
    summarise(accepted_droplets = max(accepted_droplets, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      run_date = as.Date(run_date),
      sample_type = sample_type_from_label(sample),
      sample_key = normalise_sample_key(sample)
    ) %>%
    left_join(runs, by = "run_id") %>%
    arrange(match(assay, mutation_order), run_date, well)
}

raw_paths_for_well <- function(row) {
  extract_path <- file.path(raw_root, row$archive_contents_relative_dir)
  list(
    peak_path = file.path(extract_path, "PeakData", paste0(row$well, ".ddpeakjson")),
    metadata_path = file.path(extract_path, "PeakMetaData", paste0(row$well, ".ddmetajson"))
  )
}

audit_raw_exports <- function(well_manifest) {
  well_manifest %>%
    mutate(
      extract_path = file.path(raw_root, archive_contents_relative_dir),
      peak_path = file.path(extract_path, "PeakData", paste0(well, ".ddpeakjson")),
      metadata_path = file.path(extract_path, "PeakMetaData", paste0(well, ".ddmetajson")),
      has_peak_data = file.exists(peak_path),
      has_peak_metadata = file.exists(metadata_path)
    ) %>%
    select(-extract_path)
}

threshold_entry_value <- function(entry) {
  if (is.null(entry) || length(entry) == 0L) {
    return(NA_real_)
  }
  if (is.list(entry) && !is.null(entry[[1]]) && is.list(entry[[1]]) &&
      !is.null(entry[[1]][["ThresholdValue"]])) {
    return(as.numeric(entry[[1]][["ThresholdValue"]]))
  }
  if (is.list(entry) && !is.null(entry[["ThresholdValue"]])) {
    return(as.numeric(entry[["ThresholdValue"]]))
  }
  if (is.numeric(entry) && length(entry) > 0L) {
    return(as.numeric(entry[[1]]))
  }
  NA_real_
}

threshold_rows_for_channel <- function(metadata, channel) {
  manual <- threshold_entry_value((metadata$ThresholdValues %||% list())[[channel]])
  auto <- threshold_entry_value((metadata$AutoThresholdValues %||% list())[[channel]])
  tibble(
    channel = paste0("Ch", channel),
    threshold_type = c("manual", "auto"),
    threshold = c(manual, auto)
  ) %>%
    filter(!is.na(threshold), is.finite(threshold), threshold > 0)
}

thresholds_from_metadata <- function(metadata) {
  bind_rows(
    threshold_rows_for_channel(metadata, 1L),
    threshold_rows_for_channel(metadata, 2L)
  )
}

class_from_calls <- function(ref_result, mut_result) {
  if (!all(c(ref_result, mut_result) %in% c("Negative", "Positive"))) {
    return("Gated/unassigned")
  }
  if (ref_result == "Positive" && mut_result == "Positive") {
    return("Double-positive")
  }
  if (ref_result == "Positive" && mut_result == "Negative") {
    return("Reference-only")
  }
  if (ref_result == "Negative" && mut_result == "Positive") {
    return("Mutant-only")
  }
  "Double-negative"
}

quanta_to_target_class <- function(x) {
  recode(
    as.character(x),
    "Double-negative" = "NN",
    "Reference-only" = "WT",
    "Mutant-only" = "MUT",
    "Double-positive" = "DP",
    "Gated/unassigned" = "Rain",
    "Rejected/unassigned" = "Unclassified",
    .default = "Unclassified"
  )
}

target_to_quanta_class <- function(x) {
  recode(
    as.character(x),
    "NN" = "Double-negative",
    "WT" = "Reference-only",
    "MUT" = "Mutant-only",
    "DP" = "Double-positive",
    "Rain" = "Gated/unassigned",
    "Unclassified" = "Rejected/unassigned",
    .default = "Rejected/unassigned"
  )
}

class_from_channel_calls <- function(ch1_call, ch2_call, ch1_role, ch2_role) {
  out <- rep("Rain", length(ch1_call))
  valid <- ch1_call %in% c("Negative", "Positive") &
    ch2_call %in% c("Negative", "Positive")
  if (!any(valid)) {
    return(out)
  }
  ref_call <- if (ch1_role == "reference") ch1_call else ch2_call
  mut_call <- if (ch1_role == "mutant") ch1_call else ch2_call

  out[valid] <- case_when(
    ref_call[valid] == "Positive" & mut_call[valid] == "Positive" ~ "DP",
    ref_call[valid] == "Positive" & mut_call[valid] == "Negative" ~ "WT",
    ref_call[valid] == "Negative" & mut_call[valid] == "Positive" ~ "MUT",
    ref_call[valid] == "Negative" & mut_call[valid] == "Negative" ~ "NN",
    TRUE ~ "Rain"
  )
  out
}

channel_state_for_target_class <- function(target_class, ch1_role, ch2_role) {
  ref_positive <- target_class %in% c("WT", "DP")
  mut_positive <- target_class %in% c("MUT", "DP")
  ch1_positive <- if (ch1_role == "reference") ref_positive else mut_positive
  ch2_positive <- if (ch2_role == "reference") ref_positive else mut_positive

  paste0(ifelse(ch1_positive, "P", "N"), ifelse(ch2_positive, "P", "N"))
}

target_class_from_twoddpcr <- function(channel_class, ch1_role, ch2_role) {
  channel_class <- as.character(channel_class)
  ch1_call <- ifelse(substr(channel_class, 1, 1) == "P", "Positive", "Negative")
  ch2_call <- ifelse(substr(channel_class, 2, 2) == "P", "Positive", "Negative")
  out <- class_from_channel_calls(ch1_call, ch2_call, ch1_role, ch2_role)
  out[!channel_class %in% c("NN", "NP", "PN", "PP")] <- "Rain"
  out
}

twoddpcr_class_from_target <- function(target_class, ch1_role, ch2_role) {
  channel_state_for_target_class(target_class, ch1_role, ch2_role)
}

read_well_droplets <- function(well_row) {
  paths <- raw_paths_for_well(well_row)
  if (!file.exists(paths$peak_path)) {
    stop("Missing peak data: ", paths$peak_path)
  }
  if (!file.exists(paths$metadata_path)) {
    stop("Missing peak metadata: ", paths$metadata_path)
  }

  peak <- jsonlite::fromJSON(paths$peak_path, simplifyVector = FALSE)
  metadata <- jsonlite::fromJSON(paths$metadata_path, simplifyVector = FALSE)
  amplitudes <- peak$PeakInfo$Amplitudes
  if (is.null(amplitudes) || length(amplitudes) < 2L) {
    stop("Peak data does not contain two amplitude channels: ", paths$peak_path)
  }

  targets <- metadata_targets(metadata)
  selected <- selected_target_indices(targets, well_row$assay)
  if (is.null(targets) || is.null(selected)) {
    stop("Could not select WT and mutant targets for ", well_row$run_id, " ", well_row$well)
  }

  target_names <- vapply(targets, target_name, character(1))
  target_clean <- clean_ddpcr_target(target_names)
  channels <- vapply(targets, target_channel, integer(1))
  selected_indices <- c(selected$ref, selected$mut)
  ch1_idx <- selected_indices[channels[selected_indices] == 1L][1]
  ch2_idx <- selected_indices[channels[selected_indices] == 2L][1]
  if (is.na(ch1_idx) || is.na(ch2_idx)) {
    stop("Could not map selected targets to Ch1 and Ch2 for ", well_row$run_id, " ", well_row$well)
  }

  ch1 <- as.numeric(amplitudes[[1]])
  ch2 <- as.numeric(amplitudes[[2]])
  droplet_count <- min(length(ch1), length(ch2))
  droplets <- tibble(
    droplet_index = seq.int(0L, droplet_count - 1L),
    droplet_id = paste(well_row$run_id, well_row$assay, well_row$well, droplet_index, sep = "::"),
    run_id = well_row$run_id,
    run_date = as.Date(well_row$run_date),
    assay = well_row$assay,
    well = well_row$well,
    sample = well_row$sample,
    sample_type = well_row$sample_type,
    sample_key = well_row$sample_key,
    ch1_amplitude = ch1[seq_len(droplet_count)],
    ch2_amplitude = ch2[seq_len(droplet_count)],
    Ch1.Amplitude = ch1[seq_len(droplet_count)],
    Ch2.Amplitude = ch2[seq_len(droplet_count)],
    quantaSoft_class = "Rejected/unassigned"
  )

  for (cluster in metadata$Clusters %||% list()) {
    droplet_indices <- as.integer(cluster$Droplets %||% integer(0))
    if (length(droplet_indices) == 0L) {
      next
    }
    row_indices <- droplet_indices + 1L
    row_indices <- row_indices[row_indices >= 1L & row_indices <= nrow(droplets)]
    if (length(row_indices) == 0L) {
      next
    }

    results <- as.character(unlist(cluster$Results, use.names = FALSE))
    if (length(results) < max(selected_indices) || is_true(cluster$Unassigned)) {
      droplets$quantaSoft_class[row_indices] <- "Gated/unassigned"
      next
    }

    ref_result <- results[[selected$ref]]
    mut_result <- results[[selected$mut]]
    droplets$quantaSoft_class[row_indices] <- class_from_calls(ref_result, mut_result)
  }

  droplets <- droplets %>%
    mutate(
      target_class = quanta_to_target_class(quantaSoft_class),
      ch1_role = ifelse(ch1_idx == selected$ref, "reference", "mutant"),
      ch2_role = ifelse(ch2_idx == selected$ref, "reference", "mutant"),
      ch1_target = target_clean[[ch1_idx]],
      ch2_target = target_clean[[ch2_idx]]
    )

  thresholds <- thresholds_from_metadata(metadata) %>%
    mutate(
      run_id = well_row$run_id,
      run_date = as.Date(well_row$run_date),
      assay = well_row$assay,
      well = well_row$well,
      sample = well_row$sample
    )

  list(
    droplets = droplets,
    thresholds = thresholds,
    metadata_path = paths$metadata_path,
    peak_path = paths$peak_path,
    run_id = well_row$run_id,
    run_date = as.Date(well_row$run_date),
    assay = well_row$assay,
    well = well_row$well,
    sample = well_row$sample,
    sample_type = well_row$sample_type,
    sample_key = well_row$sample_key,
    ch1_role = droplets$ch1_role[[1]],
    ch2_role = droplets$ch2_role[[1]],
    ch1_target = droplets$ch1_target[[1]],
    ch2_target = droplets$ch2_target[[1]]
  )
}

droplet_count_row <- function(parsed, method_id, method_variant, target_class,
                              status = "ok", message = NA_character_,
                              parameters = list(), weights = NULL) {
  target_class <- factor(as.character(target_class), levels = target_class_levels)
  counts <- table(target_class)
  hard_mut <- unname(counts[["MUT"]]) + unname(counts[["DP"]])
  hard_wt <- unname(counts[["WT"]]) + unname(counts[["DP"]])

  if (is.null(weights)) {
    n_mut_expected <- hard_mut
    n_wt_expected <- hard_wt
  } else {
    n_mut_expected <- sum(weights$MUT %||% 0, weights$DP %||% 0, na.rm = TRUE)
    n_wt_expected <- sum(weights$WT %||% 0, weights$DP %||% 0, na.rm = TRUE)
  }

  tibble(
    run_id = parsed$run_id,
    assay = parsed$assay,
    well = parsed$well,
    sample = parsed$sample,
    sample_type = parsed$sample_type,
    sample_key = parsed$sample_key,
    run_date = parsed$run_date,
    source_peak_data_json = parsed$peak_path,
    source_peak_metadata_json = parsed$metadata_path,
    method_id = method_id,
    method_variant = method_variant,
    n_total_droplets = length(target_class),
    n_accepted_droplets = unname(counts[["NN"]]) + unname(counts[["WT"]]) +
      unname(counts[["MUT"]]) + unname(counts[["DP"]]),
    n_nn = unname(counts[["NN"]]),
    n_wt_only = unname(counts[["WT"]]),
    n_mut_only = unname(counts[["MUT"]]),
    n_double_positive = unname(counts[["DP"]]),
    n_rain = unname(counts[["Rain"]]),
    n_unclassified = unname(counts[["Unclassified"]]),
    n_mut_expected = n_mut_expected,
    n_wt_expected = n_wt_expected,
    classification_status = status,
    classification_message = message,
    method_parameters_json = jsonlite::toJSON(parameters, auto_unbox = TRUE)
  )
}

write_method_droplets <- function(droplets, method_id) {
  out_dir <- path_v2("data", "droplets")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(droplets, file.path(out_dir, paste0(method_id, ".rds")))
}

load_shared_droplets <- function() {
  readRDS(path_v2("data", "shared_droplets.rds"))
}

load_parsed_wells <- function() {
  readRDS(path_v2("data", "parsed_wells.rds"))
}

parse_complete_wells <- function(raw_audit) {
  complete <- raw_audit %>%
    filter(has_peak_data, has_peak_metadata)

  parsed <- vector("list", nrow(complete))
  names(parsed) <- paste(complete$run_id, complete$assay, complete$well, sep = "::")

  for (i in seq_len(nrow(complete))) {
    parsed[[i]] <- read_well_droplets(complete[i, ])
  }

  parsed
}

shared_droplet_table <- function(parsed_wells) {
  bind_rows(lapply(parsed_wells, function(parsed) {
    parsed$droplets %>%
      mutate(
        source_peak_data_json = parsed$peak_path,
        source_peak_metadata_json = parsed$metadata_path
      )
  }))
}

shared_threshold_table <- function(parsed_wells) {
  bind_rows(lapply(parsed_wells, function(parsed) parsed$thresholds))
}

quanta_cluster_number <- function(target_class) {
  recode(
    as.character(target_class),
    "NN" = 1L,
    "WT" = 2L,
    "MUT" = 3L,
    "DP" = 4L,
    "Rain" = 5L,
    "Unclassified" = 6L,
    .default = 6L
  )
}

write_amplitude_csv <- function(path, droplets, include_cluster = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  out <- droplets %>%
    transmute(
      `Ch1 Amplitude` = ch1_amplitude,
      `Ch2 Amplitude` = ch2_amplitude
    )
  if (include_cluster) {
    out$Cluster <- quanta_cluster_number(droplets$target_class)
  }
  write_csv(out, path)
  tibble(path = path, rows = nrow(out), include_cluster = include_cluster)
}

write_single_channel_csv <- function(path, values) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  out <- tibble(Amplitude = as.numeric(values))
  write_csv(out, path)
  tibble(path = path, rows = nrow(out))
}

write_ddpcrclust_template <- function(path, parsed_rows) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  template <- bind_rows(lapply(parsed_rows, function(parsed) {
    tibble(
      Well = parsed$well,
      `Sample type` = parsed$sample_type,
      `# of markers` = 2L,
      `Marker 1` = parsed$ch1_target,
      `Marker 2` = parsed$ch2_target,
      `Marker 3` = "",
      `Marker 4` = ""
    )
  })) %>%
    distinct(Well, .keep_all = TRUE)

  header <- paste0(
    "> ",
    parsed_rows[[1]]$run_id,
    ", channel1=",
    parsed_rows[[1]]$ch1_target,
    ", channel2=",
    parsed_rows[[1]]$ch2_target,
    ", PRNP SNV v2"
  )
  writeLines(header, path)
  write.table(
    template,
    path,
    append = TRUE,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  tibble(path = path, rows = nrow(template))
}

reference_for_dpcp_job <- function(parsed_rows) {
  candidates <- bind_rows(lapply(parsed_rows, function(parsed) {
    tibble(
      well = parsed$well,
      sample = parsed$sample,
      sample_type = parsed$sample_type,
      accepted = nrow(parsed$droplets),
      filename = paste0(parsed$well, "_Amplitude.csv")
    )
  }))

  preferred <- candidates %>%
    filter(sample_type == "positive_control") %>%
    arrange(desc(accepted), well) %>%
    slice_head(n = 1)

  if (nrow(preferred) == 0L) {
    preferred <- candidates %>%
      filter(sample_type == "WT_control") %>%
      arrange(desc(accepted), well) %>%
      slice_head(n = 1)
  }

  if (nrow(preferred) == 0L) {
    preferred <- candidates %>%
      arrange(desc(accepted), well) %>%
      slice_head(n = 1)
  }

  preferred
}

write_dpcp_sample_table <- function(path, parsed_rows) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  reference <- reference_for_dpcp_job(parsed_rows)
  reference_file <- reference$filename[[1]]
  reference_well <- reference$well[[1]]

  sample_table <- bind_rows(lapply(parsed_rows, function(parsed) {
    tibble(
      `Sample name` = parsed$sample,
      `Chip ID/Well ID` = parsed$well,
      `No of targets` = 2L,
      `FAM target` = parsed$ch1_target,
      `Target 3` = "NA",
      `Target 4` = "NA",
      `VIC/HEX target` = parsed$ch2_target,
      Reference = ifelse(parsed$well == reference_well, "", reference_file),
      Dilution = 1
    )
  })) %>%
    distinct(`Chip ID/Well ID`, .keep_all = TRUE)

  write_csv(sample_table, path, na = "NA")
  tibble(
    path = path,
    rows = nrow(sample_table),
    reference_well = reference_well,
    reference_file = reference_file,
    reference_sample_type = reference$sample_type[[1]]
  )
}

export_twoddpcr_inputs <- function(parsed_wells) {
  groups <- split(parsed_wells, vapply(parsed_wells, function(parsed) {
    paste(parsed$assay, parsed$run_id, sep = "::")
  }, character(1)))

  bind_rows(lapply(groups, function(group) {
    assay <- group[[1]]$assay
    run_id <- group[[1]]$run_id
    out_dir <- path_v2("inputs", "twoddpcr", assay, safe_file_component(run_id))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    droplets <- bind_rows(lapply(group, function(parsed) parsed$droplets))
    manifest <- bind_rows(lapply(group, function(parsed) {
      tibble(
        run_id = parsed$run_id,
        assay = parsed$assay,
        well = parsed$well,
        sample = parsed$sample,
        sample_type = parsed$sample_type,
        droplet_rows = nrow(parsed$droplets)
      )
    }))

    saveRDS(droplets, file.path(out_dir, "droplets.rds"))
    write_csv(manifest, file.path(out_dir, "well_manifest.csv"))
    tibble(
      package = "twoddpcr",
      assay = assay,
      run_id = run_id,
      path = out_dir,
      files = 2L,
      rows = nrow(droplets)
    )
  }))
}

export_ddpcrclust_inputs <- function(parsed_wells) {
  groups <- split(parsed_wells, vapply(parsed_wells, function(parsed) {
    paste(parsed$assay, parsed$run_id, sep = "::")
  }, character(1)))

  bind_rows(lapply(groups, function(group) {
    assay <- group[[1]]$assay
    run_id <- group[[1]]$run_id
    out_dir <- path_v2("inputs", "ddPCRclust", assay, safe_file_component(run_id))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    file_rows <- bind_rows(lapply(group, function(parsed) {
      write_amplitude_csv(
        file.path(out_dir, paste0(parsed$well, "_Amplitude.csv")),
        parsed$droplets,
        include_cluster = FALSE
      ) %>%
        mutate(well = parsed$well)
    }))
    template_row <- write_ddpcrclust_template(file.path(out_dir, "ddPCR_run_template.csv"), group)

    tibble(
      package = "ddPCRclust",
      assay = assay,
      run_id = run_id,
      path = out_dir,
      files = nrow(file_rows) + 1L,
      rows = sum(file_rows$rows),
      template_path = template_row$path
    )
  }))
}

export_dpcp_inputs <- function(parsed_wells) {
  groups <- split(parsed_wells, vapply(parsed_wells, function(parsed) {
    paste(parsed$assay, parsed$run_id, sep = "::")
  }, character(1)))

  bind_rows(lapply(groups, function(group) {
    assay <- group[[1]]$assay
    run_id <- group[[1]]$run_id
    out_dir <- path_v2("inputs", "dPCP", assay, safe_file_component(run_id))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    file_rows <- bind_rows(lapply(group, function(parsed) {
      write_amplitude_csv(
        file.path(out_dir, paste0(parsed$well, "_Amplitude.csv")),
        parsed$droplets,
        include_cluster = TRUE
      ) %>%
        mutate(well = parsed$well)
    }))
    sample_table <- write_dpcp_sample_table(file.path(out_dir, "sample_table.csv"), group)

    tibble(
      package = "dPCP",
      assay = assay,
      run_id = run_id,
      path = out_dir,
      files = nrow(file_rows) + 1L,
      rows = sum(file_rows$rows),
      sample_table_path = sample_table$path,
      reference_well = sample_table$reference_well,
      reference_file = sample_table$reference_file,
      reference_sample_type = sample_table$reference_sample_type
    )
  }))
}

write_ddpcrquant_channel_dir <- function(out_dir, parsed_rows, channel) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  existing <- list.files(out_dir, full.names = TRUE)
  if (length(existing) > 0L) {
    unlink(existing)
  }

  manifest <- bind_rows(lapply(parsed_rows, function(parsed) {
    amp_path <- file.path(out_dir, paste0(parsed$well, "_Amplitude.csv"))
    amp <- parsed$droplets %>%
      transmute(
        `Ch1 Amplitude` = ch1_amplitude,
        `Ch2 Amplitude` = ch2_amplitude
      )
    write_csv(amp, amp_path)
    tibble(
      Well = parsed$well,
      Sample = parsed$sample,
      TypeAssay = ifelse(channel == "ch1", "Ch1", "Ch2"),
      Assay = paste(parsed$assay, channel, sep = "_"),
      rows = nrow(amp)
    )
  }))

  write_csv(
    manifest %>% select(Well, Sample, TypeAssay, Assay),
    file.path(out_dir, "summary.csv")
  )

  tibble(
    package = "ddpcRquant",
    assay = parsed_rows[[1]]$assay,
    run_id = parsed_rows[[1]]$run_id,
    path = out_dir,
    files = nrow(manifest) + 1L,
    rows = sum(manifest$rows),
    channel = channel
  )
}

export_one_channel_inputs <- function(parsed_wells) {
  definetherain_rows <- bind_rows(lapply(parsed_wells, function(parsed) {
    safe_run <- safe_file_component(parsed$run_id)
    rows <- list()
    for (channel in c("ch1", "ch2")) {
      values <- if (channel == "ch1") parsed$droplets$ch1_amplitude else parsed$droplets$ch2_amplitude
      out_dir <- path_v2("inputs", "definetherain", parsed$assay, safe_run, channel)
      out_path <- file.path(out_dir, paste0(parsed$well, "_amplitude.csv"))
      rows[[length(rows) + 1L]] <- write_single_channel_csv(out_path, values) %>%
        mutate(
          package = "definetherain",
          assay = parsed$assay,
          run_id = parsed$run_id,
          well = parsed$well,
          channel = channel,
          files = 1L,
          .before = 1
        )
    }
    bind_rows(rows)
  }))

  groups <- split(parsed_wells, vapply(parsed_wells, function(parsed) {
    paste(parsed$assay, parsed$run_id, sep = "::")
  }, character(1)))
  ddpc_rows <- bind_rows(lapply(groups, function(group) {
    safe_run <- safe_file_component(group[[1]]$run_id)
    bind_rows(lapply(c("ch1", "ch2"), function(channel) {
      out_dir <- path_v2("inputs", "ddpcRquant", group[[1]]$assay, safe_run, channel)
      write_ddpcrquant_channel_dir(out_dir, group, channel)
    }))
  }))

  bind_rows(
    definetherain_rows %>%
      transmute(package, assay, run_id, path, files, rows, well, channel),
    ddpc_rows %>%
      mutate(well = NA_character_) %>%
      select(package, assay, run_id, path, files, rows, well, channel)
  )
}

validate_package_inputs <- function(package_manifest) {
  validations <- list()

  ddpcrclust_row <- package_manifest %>%
    filter(package == "ddPCRclust") %>%
    slice_head(n = 1)
  if (nrow(ddpcrclust_row) == 1L) {
    validations[["ddPCRclust"]] <- tryCatch({
      files <- list.files(ddpcrclust_row$path, pattern = "_Amplitude[.]csv$", full.names = TRUE)
      template <- file.path(ddpcrclust_row$path, "ddPCR_run_template.csv")
      suppressPackageStartupMessages(library(ddPCRclust))
      read_files <- ddPCRclust::readFiles(files[1])
      read_template <- ddPCRclust::readTemplate(template)
      tibble(
        package = "ddPCRclust",
        validation = "readFiles_readTemplate",
        status = "ok",
        message = NA_character_,
        rows = nrow(read_files[[1]]),
        path = ddpcrclust_row$path
      )
    }, error = function(e) {
      tibble(
        package = "ddPCRclust",
        validation = "readFiles_readTemplate",
        status = "failed",
        message = conditionMessage(e),
        rows = NA_integer_,
        path = ddpcrclust_row$path
      )
    })
  }

  dpcp_row <- package_manifest %>%
    filter(package == "dPCP") %>%
    slice_head(n = 1)
  if (nrow(dpcp_row) == 1L) {
    validations[["dPCP"]] <- tryCatch({
      suppressPackageStartupMessages(library(dPCP))
      sample_table <- dPCP::read_sampleTable(
        file = dpcp_row$sample_table_path,
        system = "bio-rad",
        file.location = dpcp_row$path
      )
      tibble(
        package = "dPCP",
        validation = "read_sampleTable",
        status = "ok",
        message = NA_character_,
        rows = nrow(sample_table),
        path = dpcp_row$path
      )
    }, error = function(e) {
      tibble(
        package = "dPCP",
        validation = "read_sampleTable",
        status = "failed",
        message = conditionMessage(e),
        rows = NA_integer_,
        path = dpcp_row$path
      )
    })
  }

  ddpc_dir <- package_manifest %>%
    filter(package == "ddpcRquant") %>%
    arrange(assay, run_id, channel) %>%
    slice_head(n = 1)
  if (nrow(ddpc_dir) == 1L) {
    validations[["ddpcRquant"]] <- tryCatch({
      suppressPackageStartupMessages(library(dpcR))
      read_result <- dpcR::ddpcRquant(
        ddpc_dir$path,
        threshold.int = 0.999,
        reps = 1,
        blocks = 20,
        threshold.manual = FALSE
      )
      tibble(
        package = "ddpcRquant",
        validation = "ddpcRquant_smoke",
        status = "ok",
        message = NA_character_,
        rows = length(read_result),
        path = ddpc_dir$path
      )
    }, error = function(e) {
      tibble(
        package = "ddpcRquant",
        validation = "ddpcRquant_smoke",
        status = "failed",
        message = conditionMessage(e),
        rows = NA_integer_,
        path = ddpc_dir$path
      )
    })
  }

  bind_rows(validations)
}

matrix_to_cov_row <- function(assay, target_class, cov_mat, source, n_points) {
  tibble(
    assay = assay,
    target_class = target_class,
    cov_11 = cov_mat[1, 1],
    cov_12 = cov_mat[1, 2],
    cov_21 = cov_mat[2, 1],
    cov_22 = cov_mat[2, 2],
    covariance_source = source,
    covariance_n = n_points
  )
}

regularise_covariance <- function(points, pooled_cov, min_points = 30L, shrink = 0.25) {
  if (nrow(points) >= min_points) {
    qs <- apply(points, 2, quantile, probs = c(0.01, 0.99), na.rm = TRUE)
    keep <- points[, 1] >= qs[1, 1] & points[, 1] <= qs[2, 1] &
      points[, 2] >= qs[1, 2] & points[, 2] <= qs[2, 2]
    trimmed <- points[keep, , drop = FALSE]
    if (nrow(trimmed) >= min_points) {
      cov_mat <- cov(trimmed, use = "complete.obs")
      source <- "trimmed_class_covariance"
      n_used <- nrow(trimmed)
    } else {
      cov_mat <- pooled_cov
      source <- "pooled_covariance"
      n_used <- nrow(points)
    }
  } else {
    cov_mat <- pooled_cov
    source <- "pooled_covariance"
    n_used <- nrow(points)
  }

  if (any(!is.finite(cov_mat))) {
    cov_mat <- pooled_cov
    source <- "pooled_covariance_after_nonfinite"
  }

  diag_target <- diag(diag(cov_mat), 2)
  cov_mat <- (1 - shrink) * cov_mat + shrink * diag_target
  cov_mat <- cov_mat + diag(c(100, 100), 2)
  list(cov = cov_mat, source = source, n = n_used)
}

fit_control_geometry <- function(shared_droplets, min_class_points = 30L) {
  controls <- shared_droplets %>%
    filter(sample_type %in% c("NTC", "WT_control", "positive_control"))

  centre_rows <- list()
  cov_rows <- list()
  radius_rows <- list()
  validation_rows <- list()

  baseline_reference <- controls %>%
    filter(target_class == "NN") %>%
    group_by(assay) %>%
    summarise(
      reference_nn_ch1 = median(ch1_amplitude, na.rm = TRUE),
      reference_nn_ch2 = median(ch2_amplitude, na.rm = TRUE),
      .groups = "drop"
    )

  baseline_shifts <- controls %>%
    filter(target_class == "NN") %>%
    group_by(assay, run_id, well, sample, sample_type) %>%
    summarise(
      well_nn_ch1 = median(ch1_amplitude, na.rm = TRUE),
      well_nn_ch2 = median(ch2_amplitude, na.rm = TRUE),
      n_nn = n(),
      .groups = "drop"
    ) %>%
    left_join(baseline_reference, by = "assay") %>%
    mutate(
      shift_ch1 = reference_nn_ch1 - well_nn_ch1,
      shift_ch2 = reference_nn_ch2 - well_nn_ch2
    )

  for (assay_name in mutation_order) {
    assay_controls <- controls %>% filter(assay == assay_name)
    if (nrow(assay_controls) == 0L) {
      next
    }

    roles <- assay_controls %>%
      summarise(ch1_role = first(ch1_role), ch2_role = first(ch2_role), .groups = "drop")
    ch1_role <- roles$ch1_role[[1]]
    ch2_role <- roles$ch2_role[[1]]

    observed <- assay_controls %>%
      filter(target_class %in% c("NN", "WT", "MUT", "DP")) %>%
      group_by(target_class) %>%
      summarise(
        observed_n = n(),
        centre_ch1 = median(ch1_amplitude, na.rm = TRUE),
        centre_ch2 = median(ch2_amplitude, na.rm = TRUE),
        .groups = "drop"
      )

    centre_for <- function(class_name, coordinate) {
      value <- observed %>%
        filter(target_class == class_name) %>%
        pull(all_of(coordinate))
      if (length(value) == 1L && is.finite(value)) value else NA_real_
    }
    n_for <- function(class_name) {
      value <- observed %>% filter(target_class == class_name) %>% pull(observed_n)
      if (length(value) == 1L) value else 0L
    }

    nn <- c(centre_for("NN", "centre_ch1"), centre_for("NN", "centre_ch2"))
    wt <- c(centre_for("WT", "centre_ch1"), centre_for("WT", "centre_ch2"))
    mut <- c(centre_for("MUT", "centre_ch1"), centre_for("MUT", "centre_ch2"))
    dp_observed <- c(centre_for("DP", "centre_ch1"), centre_for("DP", "centre_ch2"))

    if (any(!is.finite(nn))) {
      nn <- c(
        median(assay_controls$ch1_amplitude, na.rm = TRUE),
        median(assay_controls$ch2_amplitude, na.rm = TRUE)
      )
    }
    if (any(!is.finite(wt))) {
      wt <- nn
      wt[ifelse(ch1_role == "reference", 1, 2)] <-
        quantile(assay_controls[[ifelse(ch1_role == "reference", "ch1_amplitude", "ch2_amplitude")]],
                 0.995, na.rm = TRUE, names = FALSE)
    }
    if (any(!is.finite(mut))) {
      mut <- nn
      mut[ifelse(ch1_role == "mutant", 1, 2)] <-
        quantile(assay_controls[[ifelse(ch1_role == "mutant", "ch1_amplitude", "ch2_amplitude")]],
                 0.995, na.rm = TRUE, names = FALSE)
    }

    dp_vector <- wt + mut - nn
    dp <- if (n_for("DP") >= min_class_points && all(is.finite(dp_observed))) {
      dp_observed
    } else {
      dp_vector
    }
    dp_source <- if (n_for("DP") >= min_class_points && all(is.finite(dp_observed))) {
      "observed_positive_control"
    } else {
      "vector_sum"
    }

    centres <- tibble(
      assay = assay_name,
      target_class = c("NN", "WT", "MUT", "DP"),
      centre_ch1 = c(nn[1], wt[1], mut[1], dp[1]),
      centre_ch2 = c(nn[2], wt[2], mut[2], dp[2]),
      centre_source = c(
        ifelse(n_for("NN") >= min_class_points, "observed_control", "fallback_median"),
        ifelse(n_for("WT") >= min_class_points, "observed_control", "channel_quantile"),
        ifelse(n_for("MUT") >= min_class_points, "observed_control", "channel_quantile"),
        dp_source
      ),
      observed_n = c(n_for("NN"), n_for("WT"), n_for("MUT"), n_for("DP")),
      ch1_role = ch1_role,
      ch2_role = ch2_role
    ) %>%
      mutate(twoddpcr_class = channel_state_for_target_class(target_class, .env$ch1_role, .env$ch2_role))
    centre_rows[[assay_name]] <- centres

    pooled_points <- assay_controls %>%
      filter(target_class %in% c("NN", "WT", "MUT", "DP")) %>%
      select(ch1_amplitude, ch2_amplitude) %>%
      as.matrix()
    pooled_cov <- cov(pooled_points, use = "complete.obs")
    if (any(!is.finite(pooled_cov))) {
      pooled_cov <- diag(apply(pooled_points, 2, sd, na.rm = TRUE)^2, 2)
    }
    pooled_cov <- pooled_cov + diag(c(100, 100), 2)

    covs <- list()
    for (class_name in c("NN", "WT", "MUT", "DP")) {
      class_points <- assay_controls %>%
        filter(target_class == class_name) %>%
        select(ch1_amplitude, ch2_amplitude) %>%
        as.matrix()
      reg <- regularise_covariance(class_points, pooled_cov, min_points = min_class_points)
      covs[[class_name]] <- reg$cov
      cov_rows[[paste(assay_name, class_name, sep = "::")]] <-
        matrix_to_cov_row(assay_name, class_name, reg$cov, reg$source, reg$n)

      centre <- centres %>% filter(target_class == class_name)
      if (nrow(class_points) >= min_class_points) {
        d2 <- mahalanobis(class_points, center = c(centre$centre_ch1, centre$centre_ch2), cov = reg$cov)
        radius <- unname(quantile(d2[is.finite(d2)], 0.995, na.rm = TRUE))
        radius_source <- "empirical_control_0.995"
      } else {
        radius <- qchisq(0.995, df = 2)
        radius_source <- "chisq_0.995"
      }
      if (!is.finite(radius) || radius <= 0) {
        radius <- qchisq(0.995, df = 2)
        radius_source <- "chisq_0.995_after_invalid"
      }
      radius_rows[[paste(assay_name, class_name, sep = "::")]] <- tibble(
        assay = assay_name,
        target_class = class_name,
        radius_d2 = radius,
        radius_source = radius_source
      )
    }

    assigned <- classify_by_control_geometry(
      droplets = assay_controls,
      centres = centres,
      covariance_rows = bind_rows(cov_rows) %>% filter(assay == assay_name),
      gate_radii = bind_rows(radius_rows) %>% filter(assay == assay_name)
    )

    validation_rows[[assay_name]] <- assigned %>%
      count(assay, sample_type, quanta_class = target_class, assigned_class, name = "droplets") %>%
      arrange(assay, sample_type, quanta_class, assigned_class)
  }

  list(
    centres = bind_rows(centre_rows),
    covariance_rows = bind_rows(cov_rows),
    gate_radii = bind_rows(radius_rows),
    baseline_shifts = baseline_shifts,
    control_validation = bind_rows(validation_rows)
  )
}

covariance_matrix_for <- function(covariance_rows, assay_name, class_name) {
  row <- covariance_rows %>% filter(assay == assay_name, target_class == class_name) %>% slice(1)
  matrix(c(row$cov_11, row$cov_12, row$cov_21, row$cov_22), nrow = 2, byrow = TRUE)
}

classify_by_control_geometry <- function(droplets, centres, covariance_rows, gate_radii) {
  out <- droplets
  x <- as.matrix(out[, c("ch1_amplitude", "ch2_amplitude")])

  distances <- sapply(seq_len(nrow(centres)), function(i) {
    class_name <- centres$target_class[[i]]
    cov_mat <- covariance_matrix_for(covariance_rows, centres$assay[[i]], class_name)
    mahalanobis(
      x,
      center = c(centres$centre_ch1[[i]], centres$centre_ch2[[i]]),
      cov = cov_mat
    )
  })
  colnames(distances) <- centres$target_class
  winner <- max.col(-distances, ties.method = "first")
  assigned <- centres$target_class[winner]
  winning_d2 <- distances[cbind(seq_len(nrow(out)), winner)]
  radius <- gate_radii$radius_d2[match(assigned, gate_radii$target_class)]
  assigned[!is.finite(winning_d2) | winning_d2 > radius] <- "Rain"

  out$assigned_class <- assigned
  out$assigned_d2 <- winning_d2
  out
}

ellipse_points_from_geometry <- function(centres, covariance_rows, gate_radii, points = 120L) {
  bind_rows(lapply(seq_len(nrow(centres)), function(i) {
    centre <- centres[i, ]
    cov_mat <- covariance_matrix_for(covariance_rows, centre$assay, centre$target_class)
    radius <- gate_radii$radius_d2[
      gate_radii$assay == centre$assay & gate_radii$target_class == centre$target_class
    ][1]
    eig <- eigen(cov_mat, symmetric = TRUE)
    theta <- seq(0, 2 * pi, length.out = points)
    circle <- rbind(cos(theta), sin(theta))
    coords <- c(centre$centre_ch1, centre$centre_ch2) +
      t(eig$vectors %*% diag(sqrt(pmax(eig$values, 0) * radius), 2) %*% circle)
    tibble(
      assay = centre$assay,
      target_class = centre$target_class,
      ch1_amplitude = coords[, 1],
      ch2_amplitude = coords[, 2],
      ellipse_id = paste(centre$assay, centre$target_class, sep = "::")
    )
  }))
}

save_control_geometry_plots <- function(shared_droplets, geometry) {
  out_dir <- path_v2("plots", "individual", "control_geometry")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ellipses <- ellipse_points_from_geometry(
    geometry$centres,
    geometry$covariance_rows,
    geometry$gate_radii
  )

  plot_rows <- list()
  for (assay_name in mutation_order) {
    controls <- shared_droplets %>%
      filter(assay == assay_name, sample_type %in% c("NTC", "WT_control", "positive_control")) %>%
      group_by(sample_type, target_class) %>%
      group_modify(function(.x, .y) {
        .x %>% slice_sample(n = min(nrow(.x), 1200L))
      }) %>%
      ungroup()
    centres <- geometry$centres %>% filter(assay == assay_name)
    assay_ellipses <- ellipses %>% filter(assay == assay_name)

    plot <- ggplot(controls, aes(ch1_amplitude, ch2_amplitude, colour = target_class)) +
      geom_point(size = 0.18, alpha = 0.35, stroke = 0) +
      geom_path(
        data = assay_ellipses,
        aes(ch1_amplitude, ch2_amplitude, group = ellipse_id, colour = target_class),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_point(
        data = centres,
        aes(centre_ch1, centre_ch2, fill = target_class),
        inherit.aes = FALSE,
        shape = 21,
        colour = "black",
        size = 2
      ) +
      facet_wrap(~sample_type, nrow = 1) +
      coord_cartesian(
        xlim = c(0, quantile(controls$ch1_amplitude, 0.999, na.rm = TRUE) * 1.05),
        ylim = c(0, quantile(controls$ch2_amplitude, 0.999, na.rm = TRUE) * 1.05),
        expand = FALSE
      ) +
      labs(
        title = paste(assay_name, "control-derived geometry"),
        x = "Channel 1 amplitude",
        y = "Channel 2 amplitude",
        colour = "Class"
      ) +
      theme_bw(base_size = 8) +
      theme(panel.grid.minor = element_blank(), legend.position = "bottom")

    svg_path <- file.path(out_dir, paste0(assay_name, "_control_geometry.svg"))
    pdf_path <- file.path(out_dir, paste0(assay_name, "_control_geometry.pdf"))
    grDevices::svg(svg_path, width = 8.5, height = 3.2, onefile = FALSE)
    print(plot)
    grDevices::dev.off()
    grDevices::pdf(pdf_path, width = 8.5, height = 3.2, useDingbats = FALSE)
    print(plot)
    grDevices::dev.off()

    plot_rows[[assay_name]] <- tibble(
      assay = assay_name,
      svg_path = svg_path,
      pdf_path = pdf_path
    )
  }

  bind_rows(plot_rows)
}
