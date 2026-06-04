library(tidyverse)
library(binom)

# Diagnostic only. This script compares the current LoB blank model
# (WT_control + NTC) with a WT_control-only blank model. It does not alter the
# main analysis code; it writes comparison tables under analysis/validation.

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  }

  source_files <- unlist(lapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) {
      character(0)
    } else {
      frame$ofile
    }
  }))

  if (length(source_files) > 0) {
    return(normalizePath(source_files[length(source_files)], winslash = "/", mustWork = TRUE))
  }

  stop("Could not determine this script's path.")
}

script_dir <- dirname(get_script_path())
main_script_path <- file.path(script_dir, "create_snv_dataframe_corrected.R")

analysis_env <- new.env(parent = globalenv())
sys.source(main_script_path, envir = analysis_env)

build_lob_outputs <- function(blank_samples, label) {
  blanks <- analysis_env$data.mut.controls %>%
    filter(Target %in% analysis_env$mutation.list, Sample %in% blank_samples) %>%
    filter(AcceptedDroplets >= 10000) %>%
    transmute(
      plate = Date,
      assay = ExperimentType,
      blank_sample = Sample,
      n = AcceptedDroplets,
      x = Positives
    )

  if (nrow(blanks) == 0) {
    stop("No qualifying blank wells for model: ", label)
  }

  blank_pooled <- blanks %>%
    group_by(plate, assay) %>%
    summarise(
      x_blank = sum(x, na.rm = TRUE),
      n_blank = sum(n, na.rm = TRUE),
      n_wells_blank = n(),
      blank_samples = paste(sort(unique(blank_sample)), collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(p0_upper = binom.confint(x_blank, n_blank, methods = "exact")$upper)

  assay_fallback <- blank_pooled %>%
    group_by(assay) %>%
    summarise(
      x_blank = sum(x_blank),
      n_blank = sum(n_blank),
      .groups = "drop"
    ) %>%
    mutate(p0_upper_fallback = binom.confint(x_blank, n_blank, methods = "exact")$upper) %>%
    select(assay, p0_upper_fallback)

  sample_lob_p0 <- analysis_env$data.mut.brief %>%
    select(Sample, assay = ExperimentType, plate = Date) %>%
    distinct() %>%
    left_join(blank_pooled, by = c("plate", "assay")) %>%
    left_join(assay_fallback, by = "assay") %>%
    mutate(p0_per_plate = coalesce(p0_upper, p0_upper_fallback)) %>%
    mutate(p0_per_plate = ifelse(is.na(p0_per_plate), 0, p0_per_plate)) %>%
    group_by(Sample, assay) %>%
    slice_max(p0_per_plate, n = 1, with_ties = FALSE) %>%
    summarise(
      p0_use = first(p0_per_plate),
      plate_p0_max = first(plate),
      .groups = "drop"
    )

  sample_lob <- analysis_env$merged %>%
    rename(assay = ExperimentType, n_tot = Drop_pool, x_mut = Pos_pool) %>%
    left_join(sample_lob_p0, by = c("Sample", "assay")) %>%
    mutate(p0_use = ifelse(is.na(p0_use), 0, p0_use)) %>%
    rowwise() %>%
    mutate(
      LoB_count = qbinom(0.95, size = n_tot, prob = p0_use),
      LoB_FA = ifelse(n_tot > 0, analysis_env$fa_scale * LoB_count / n_tot, NA_real_),
      detected_LoB = x_mut > LoB_count
    ) %>%
    ungroup() %>%
    select(Sample, assay, p0_use, plate_p0_max, LoB_count, LoB_FA, detected_LoB)

  plate_mapping <- analysis_env$data.mut.brief %>%
    mutate(
      code = sub("_.*$", "", Sample),
      assay = ExperimentType,
      plate = Date
    ) %>%
    distinct(code, assay, plate)

  p0_per_participant <- plate_mapping %>%
    left_join(blank_pooled %>% select(plate, assay, p0_upper), by = c("plate", "assay")) %>%
    left_join(assay_fallback, by = "assay") %>%
    mutate(p0_per_plate = coalesce(p0_upper, p0_upper_fallback)) %>%
    mutate(p0_per_plate = ifelse(is.na(p0_per_plate), 0, p0_per_plate)) %>%
    group_by(code, assay) %>%
    slice_max(p0_per_plate, n = 1, with_ties = FALSE) %>%
    summarise(
      p0_max = first(p0_per_plate),
      plate_p0_max = first(plate),
      .groups = "drop"
    )

  participant_lob <- analysis_env$pooled_participant_assay %>%
    select(code, assay, pos_total, drop_total) %>%
    left_join(p0_per_participant, by = c("code", "assay")) %>%
    mutate(p0_max = ifelse(is.na(p0_max), 0, p0_max)) %>%
    mutate(
      LoB_count = qbinom(0.95, size = drop_total, prob = p0_max),
      LoB_FA = analysis_env$fa_scale * LoB_count / drop_total,
      detected_above_lob = pos_total > LoB_count
    ) %>%
    select(code, assay, p0_max, plate_p0_max, LoB_count, LoB_FA, detected_above_lob)

  list(
    label = label,
    blanks = blanks,
    blank_pooled = blank_pooled,
    assay_fallback = assay_fallback,
    sample_lob = sample_lob,
    participant_lob = participant_lob
  )
}

summarise_compare <- function(compare_data, level) {
  tibble(
    level = level,
    n_rows = nrow(compare_data),
    n_lob_count_changed = sum(compare_data$current_lob_count != compare_data$wt_only_lob_count, na.rm = TRUE),
    n_lob_count_lower_with_wt_only = sum(compare_data$wt_only_lob_count < compare_data$current_lob_count, na.rm = TRUE),
    n_lob_count_higher_with_wt_only = sum(compare_data$wt_only_lob_count > compare_data$current_lob_count, na.rm = TRUE),
    current_lob_positive = sum(compare_data$current_detected_lob, na.rm = TRUE),
    wt_only_lob_positive = sum(compare_data$wt_only_detected_lob, na.rm = TRUE),
    n_lob_flag_changed = sum(compare_data$current_detected_lob != compare_data$wt_only_detected_lob, na.rm = TRUE),
    current_lod_lob_positive = sum(compare_data$current_lod_lob_positive, na.rm = TRUE),
    wt_only_lod_lob_positive = sum(compare_data$wt_only_lod_lob_positive, na.rm = TRUE),
    n_lod_lob_flag_changed = sum(compare_data$current_lod_lob_positive != compare_data$wt_only_lod_lob_positive, na.rm = TRUE)
  )
}

current <- build_lob_outputs(c("WT_control", "NTC"), "current_wt_plus_ntc")
wt_only <- build_lob_outputs("WT_control", "wt_only")

sample_lod <- analysis_env$final %>%
  select(Sample, assay = ExperimentType, detected_LoD)

sample_compare <- current$sample_lob %>%
  rename(
    current_p0 = p0_use,
    current_plate_p0_max = plate_p0_max,
    current_lob_count = LoB_count,
    current_lob_fa = LoB_FA,
    current_detected_lob = detected_LoB
  ) %>%
  inner_join(
    wt_only$sample_lob %>%
      rename(
        wt_only_p0 = p0_use,
        wt_only_plate_p0_max = plate_p0_max,
        wt_only_lob_count = LoB_count,
        wt_only_lob_fa = LoB_FA,
        wt_only_detected_lob = detected_LoB
      ),
    by = c("Sample", "assay")
  ) %>%
  left_join(sample_lod, by = c("Sample", "assay")) %>%
  mutate(
    current_lod_lob_positive = current_detected_lob & detected_LoD,
    wt_only_lod_lob_positive = wt_only_detected_lob & detected_LoD
  )

participant_lod <- analysis_env$pooled_participant_assay %>%
  select(code, assay, detected_above_lod)

participant_compare <- current$participant_lob %>%
  rename(
    current_p0 = p0_max,
    current_plate_p0_max = plate_p0_max,
    current_lob_count = LoB_count,
    current_lob_fa = LoB_FA,
    current_detected_lob = detected_above_lob
  ) %>%
  inner_join(
    wt_only$participant_lob %>%
      rename(
        wt_only_p0 = p0_max,
        wt_only_plate_p0_max = plate_p0_max,
        wt_only_lob_count = LoB_count,
        wt_only_lob_fa = LoB_FA,
        wt_only_detected_lob = detected_above_lob
      ),
    by = c("code", "assay")
  ) %>%
  left_join(participant_lod, by = c("code", "assay")) %>%
  mutate(
    current_lod_lob_positive = current_detected_lob & detected_above_lod,
    wt_only_lod_lob_positive = wt_only_detected_lob & detected_above_lod
  )

p0_compare <- current$blank_pooled %>%
  select(plate, assay, current_x_blank = x_blank, current_n_blank = n_blank,
         current_n_wells_blank = n_wells_blank, current_p0_upper = p0_upper) %>%
  full_join(
    wt_only$blank_pooled %>%
      select(plate, assay, wt_only_x_blank = x_blank, wt_only_n_blank = n_blank,
             wt_only_n_wells_blank = n_wells_blank, wt_only_p0_upper = p0_upper),
    by = c("plate", "assay")
  )

summary_table <- bind_rows(
  summarise_compare(sample_compare, "sample_region"),
  summarise_compare(participant_compare, "participant_pooled")
)

validation_dir <- analysis_env$validation_dir
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(sample_compare, file.path(validation_dir, "wt_only_lob_sample_region_compare.csv"))
write_csv(participant_compare, file.path(validation_dir, "wt_only_lob_participant_compare.csv"))
write_csv(p0_compare, file.path(validation_dir, "wt_only_lob_blank_p0_compare.csv"))
write_csv(summary_table, file.path(validation_dir, "wt_only_lob_sensitivity_summary.csv"))

print(summary_table)
