library(tidyverse)
library(readxl)

bootstrap_script_path <- local({
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) > 0L) {
    return(sub("^--file=", "", script_arg[[1]]))
  }
  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) == 0L) {
    stop("Could not determine script path. Run with Rscript or source().")
  }
  tail(frame_files, n = 1)[[1]]
})

script_dir <- dirname(normalizePath(bootstrap_script_path, winslash = "/", mustWork = TRUE))
repo_root <- script_dir
repeat {
  if (dir.exists(file.path(repo_root, "manuscript")) &&
      dir.exists(file.path(repo_root, "results"))) {
    break
  }
  parent_dir <- dirname(repo_root)
  if (identical(parent_dir, repo_root)) {
    stop("Could not locate repo root from: ", bootstrap_script_path)
  }
  repo_root <- parent_dir
}

repo_path <- function(...) file.path(repo_root, ...)
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# ---- user inputs ----
xlsx_path <- repo_path("results", "ddPCR", "SNV_data_final.xlsx")
output_dir <- ensure_dir(repo_path("manuscript", "tables", "ddpcr_sample_number"))


sheet_name <- "Sheet 1"

# Optional: keep only a specific cohort (set to NULL to keep everything)
participants_keep <- NULL

# Optional: exclude pooled samples (recommended if pooled is not meant to count as a "sample")
exclude_pooled <- TRUE

# If you later decide you need a QC filter, put it here (currently "keep all")
# Example placeholder: qc_pass <- n_total_droplets >= 10000
qc_filter <- function(.df) .df

# ---- helper ----
safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

format_vaf <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    format(round(x, 3), nsmall = 3, scientific = FALSE, trim = TRUE)
  )
}

as_true_flag <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

# ---- read + basic filtering ----
ddpcr <- read_excel(xlsx_path, sheet = sheet_name) %>%
  mutate(
    participant = as.character(participant),
    group = as.character(group),
    mutation = as.character(mutation),
    brain_region = as.character(brain_region)
  )

required_lob_lod_cols <- c("detected_above_LoB", "detected_above_LoD")
missing_lob_lod_cols <- setdiff(required_lob_lod_cols, names(ddpcr))
if (length(missing_lob_lod_cols) > 0) {
  stop("SNV_data_final.xlsx is missing LoB/LoD columns: ",
       paste(missing_lob_lod_cols, collapse = ", "))
}

mutation_order <- c("D178N", "E200K", "P102L")

lob_lod_status_rows <- ddpcr %>%
  mutate(
    pass_lob_lod = as_true_flag(detected_above_LoB) & as_true_flag(detected_above_LoD),
    known_heterozygous_e200k = participant == "CJD30" & mutation == "E200K",
    .mutation_order = match(mutation, mutation_order),
    .status_label = if_else(known_heterozygous_e200k, paste0(mutation, "*"), mutation)
  ) %>%
  filter(pass_lob_lod)

participant_lob_lod_status <- lob_lod_status_rows %>%
  distinct(participant, group, .mutation_order, .status_label) %>%
  arrange(participant, .mutation_order, .status_label) %>%
  group_by(participant, group) %>%
  summarise(`LoB+LoD pass` = paste(.status_label, collapse = "; "), .groups = "drop")

group_lob_lod_status <- lob_lod_status_rows %>%
  distinct(group, .mutation_order, .status_label) %>%
  arrange(group, .mutation_order, .status_label) %>%
  group_by(group) %>%
  summarise(`LoB+LoD pass` = paste(.status_label, collapse = "; "), .groups = "drop")

# ---- per-participant: number of analysed samples (any assay) ----
analysed_samples <- ddpcr %>%
  distinct(participant, group, brain_region) %>%
  count(participant, group, name = "analysed_samples")

# ---- per-participant, per-mutation: counts + maxima ----
by_mut <- ddpcr %>%
  group_by(participant, group, mutation) %>%
  summarise(
    analysed = n_distinct(brain_region),
    max_maf = safe_max(fractional_abundance),
    .groups = "drop"
  )

counts_wide <- by_mut %>%
  select(participant, group, mutation, analysed) %>%
  pivot_wider(
    names_from = mutation,
    values_from = analysed,
    names_prefix = "analysed_"
  ) %>%
  mutate(across(starts_with("analysed_"), ~ replace_na(.x, 0L)))

max_wide <- by_mut %>%
  select(participant, group, mutation, max_maf) %>%
  pivot_wider(
    names_from = mutation,
    values_from = max_maf,
    names_prefix = "max_maf_"
  )

# ---- assemble patient-level summary ----
ddpcr_by_patient <- analysed_samples %>%
  left_join(counts_wide, by = c("participant", "group")) %>%
  left_join(max_wide, by = c("participant", "group")) %>%
  left_join(participant_lob_lod_status, by = c("participant", "group")) %>%
  mutate(
    name = participant,
    `LoB+LoD pass` = replace_na(`LoB+LoD pass`, "None")
  ) %>%
  rowwise() %>%
  mutate(
    max_maf_overall = if_else(
      name == "CJD30",
      safe_max(c(max_maf_D178N, max_maf_P102L)),
      safe_max(c(max_maf_D178N, max_maf_E200K, max_maf_P102L))
    )
  ) %>%
  ungroup() %>%
  select(
    name,
    analysed_samples,
    analysed_D178N, analysed_E200K, analysed_P102L,
    max_maf_D178N, max_maf_E200K, max_maf_P102L,
    max_maf_overall,
    `LoB+LoD pass`,
    group
  )

# ---- Exclude CJD30 from E200K analysis, as it's heterozygous (set to NA) ----

#ddpcr_by_patient$max_maf_E200K[ddpcr_by_patient$name == "CJD30"] <- NA_real_
# this avoids setting the maximum E200K to 50%


# ---- group totals ("All CJD", "All controls"), excluding E200K for CJD30 (as it's heterozygous) ----
group_totals <- ddpcr_by_patient %>%
  group_by(group) %>%
  summarise(
    analysed_samples = sum(analysed_samples, na.rm = TRUE),
    analysed_D178N = sum(analysed_D178N, na.rm = TRUE),
    analysed_E200K = sum(analysed_E200K, na.rm = TRUE),
    analysed_P102L = sum(analysed_P102L, na.rm = TRUE),
    max_maf_D178N = safe_max(max_maf_D178N),
    max_maf_E200K = safe_max(replace(max_maf_E200K, name == "CJD30", NA_real_)),
    max_maf_P102L = safe_max(max_maf_P102L),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    max_maf_overall = safe_max(c(max_maf_D178N, max_maf_E200K, max_maf_P102L)),
    name = case_when(
      group == "prion" ~ "All CJD",
      group == "control" ~ "All controls",
      TRUE ~ paste("All", group)
    )
  ) %>%
  ungroup() %>%
  left_join(group_lob_lod_status, by = "group") %>%
  mutate(`LoB+LoD pass` = replace_na(`LoB+LoD pass`, "None")) %>%
  select(
    name,
    analysed_samples,
    analysed_D178N, analysed_E200K, analysed_P102L,
    max_maf_D178N, max_maf_E200K, max_maf_P102L,
    max_maf_overall,
    `LoB+LoD pass`,
    group
  )

# ---- final output: order rows ----
ddpcr_summary <- bind_rows(
  ddpcr_by_patient %>%
    filter(group == "prion") %>%
    mutate(.ord = readr::parse_number(name)) %>%
    arrange(.ord, name) %>%
    select(-.ord),

  ddpcr_by_patient %>%
    filter(group == "control") %>%
    mutate(.ord = readr::parse_number(name)) %>%
    arrange(.ord, name) %>%
    select(-.ord),

  group_totals %>%
    arrange(match(group, c("prion", "control")))
) %>%
  select(-group)

# ---- round numbers to three decimal places ----
ddpcr_summary <- ddpcr_summary %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))


# Latex

library(knitr)

ddpcr_summary_tex <- ddpcr_summary %>%
  mutate(
    max_maf_D178N = format_vaf(max_maf_D178N),
    max_maf_E200K = if_else(
      name == "CJD30",
      paste0(format_vaf(max_maf_E200K), "*"),
      format_vaf(max_maf_E200K)
    ),
    max_maf_P102L = format_vaf(max_maf_P102L),
    max_maf_overall = format_vaf(max_maf_overall)
  )

kable(ddpcr_summary_tex,
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      align = "lrrrrrrrrl")

writeLines(
  kable(
    ddpcr_summary_tex,
    format = "latex",
    booktabs = FALSE,
    longtable = FALSE,
    align = "lrrrrrrrrl"
  ),
  file.path(output_dir, "ddpcr_sample_number_basic.tex")
)



# ---- export as CSV ----
write.csv(ddpcr_summary, file.path(output_dir, "ddPCR_sample_number.csv"), row.names = FALSE)
