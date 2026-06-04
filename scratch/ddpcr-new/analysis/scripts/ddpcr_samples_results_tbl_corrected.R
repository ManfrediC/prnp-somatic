library(tidyverse)
library(openxlsx)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
scratch_root <- file.path(project_root, "scratch", "ddpcr-new")
analysis_root <- file.path(scratch_root, "analysis")
input_path <- file.path(analysis_root, "results", "ddPCR", "SNV_data_final.xlsx")
current_input_path <- file.path(project_root, "results", "ddPCR", "SNV_data_final.xlsx")
output_dir <- file.path(scratch_root, "tables", "sample_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)
region_labels <- c(
  bg = "basal ganglia",
  cb = "cerebellum",
  fr = "frontal cortex",
  hc = "hippocampus",
  ps = "pons",
  sn = "substantia nigra",
  th = "thalamus"
)
region_list <- unname(region_labels[c("bg", "cb", "hc", "fr", "sn", "th")])

as_true_flag <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

participant_order <- function(x) {
  tibble(participant = unique(as.character(x))) %>%
    mutate(
      .group = case_when(
        grepl("^CJD", participant) ~ 1L,
        grepl("^Control", participant) ~ 2L,
        TRUE ~ 3L
      ),
      .number = readr::parse_number(participant)
    ) %>%
    arrange(.group, .number, participant) %>%
    pull(participant)
}

format_decimal <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    format(round(x, 3), nsmall = 3, scientific = FALSE, trim = TRUE)
  )
}

escape_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

write_html_table <- function(table_data, mask_data, output_path) {
  headers <- paste0("<th>", escape_html(names(table_data)), "</th>", collapse = "")
  body <- vapply(seq_len(nrow(table_data)), function(i) {
    cells <- vapply(names(table_data), function(col) {
      value <- table_data[[col]][i]
      mask_value <- FALSE
      if (col %in% names(mask_data)) {
        mask_value <- isTRUE(mask_data[[col]][i])
      }
      style <- if (mask_value) " style=\"background:#e5e7eb\"" else ""
      paste0("<td", style, ">", escape_html(ifelse(is.na(value), "", value)), "</td>")
    }, character(1))
    paste0("<tr>", paste0(cells, collapse = ""), "</tr>")
  }, character(1))

  html <- c(
    "<!doctype html>",
    "<html><head><meta charset=\"utf-8\">",
    "<style>",
    "body{font-family:Arial,sans-serif;margin:24px;color:#111827}",
    "table{border-collapse:collapse;font-size:12px}",
    "th,td{border:1px solid #d1d5db;padding:4px 6px;vertical-align:top}",
    "th{background:#f3f4f6;position:sticky;top:0}",
    "caption{caption-side:top;text-align:left;font-weight:700;margin-bottom:8px}",
    "</style></head><body>",
    "<table><caption>Corrected ddPCR sample-region results. Grey cells are above LoB.</caption>",
    paste0("<thead><tr>", headers, "</tr></thead>"),
    paste0("<tbody>", paste0(body, collapse = "\n"), "</tbody>"),
    "</table></body></html>"
  )

  writeLines(html, output_path)
}

df_raw <- openxlsx::read.xlsx(input_path) %>%
  as_tibble()

df <- df_raw %>%
  mutate(
    brain_region = recode(brain_region, !!!region_labels),
    detected_above_LoB = as_true_flag(detected_above_LoB),
    detected_above_LoD = as_true_flag(detected_above_LoD),
    is_pooled = as_true_flag(is_pooled),
    LoD = lod_cut[mutation]
  ) %>%
  rename(LoB = lob_fa) %>%
  relocate(LoD, .before = detected_above_LoD) %>%
  select(-group, -code, -lob_count)

df_for_workbook <- df %>%
  mutate(
    detected_above_LoB = if_else(detected_above_LoB, "Yes", "No"),
    detected_above_LoD = if_else(detected_above_LoD, "Yes", "No")
  )

openxlsx::write.xlsx(
  df_for_workbook,
  file.path(output_dir, "ddPCR_results_table_corrected.xlsx"),
  overwrite = TRUE
)

comparison_keys <- c("participant", "code", "brain_region", "mutation")
if (file.exists(current_input_path)) {
  current_df <- openxlsx::read.xlsx(current_input_path) %>%
    as_tibble() %>%
    select(
      all_of(comparison_keys),
      is_pooled_current = is_pooled,
      fractional_abundance_current = fractional_abundance,
      ci_low_current = ci_low,
      ci_high_current = ci_high,
      lob_fa_current = lob_fa,
      detected_above_LoD_current = detected_above_LoD
    ) %>%
    mutate(
      detected_above_LoD_current = as_true_flag(detected_above_LoD_current),
      is_pooled_current = as_true_flag(is_pooled_current)
    )

  changed_rows <- df_raw %>%
    left_join(current_df, by = comparison_keys) %>%
    mutate(
      detected_above_LoD = as_true_flag(detected_above_LoD),
      is_pooled = as_true_flag(is_pooled),
      lod_flag_changed = detected_above_LoD != detected_above_LoD_current
    ) %>%
    filter(is_pooled | lod_flag_changed) %>%
    arrange(mutation, participant, brain_region) %>%
    select(
      all_of(comparison_keys),
      is_pooled,
      fractional_abundance_current,
      fractional_abundance_corrected = fractional_abundance,
      ci_low_current,
      ci_low_corrected = ci_low,
      ci_high_current,
      ci_high_corrected = ci_high,
      lob_fa_current,
      lob_fa_corrected = lob_fa,
      detected_above_LoD_current,
      detected_above_LoD_corrected = detected_above_LoD,
      lod_flag_changed
    )

  readr::write_csv(
    changed_rows,
    file.path(output_dir, "sample_region_pooled_and_lod_changed_rows_corrected.csv")
  )
}

table_base <- df %>%
  filter(brain_region != "pons") %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3)),
    participant = factor(participant, levels = participant_order(participant)),
    fa_ci = paste0(
      format_decimal(fractional_abundance),
      " (",
      format_decimal(ci_low),
      "-",
      format_decimal(ci_high),
      ")"
    )
  )

mytable <- table_base %>%
  select(participant, histotype, mutation, brain_region, fa_ci) %>%
  pivot_wider(names_from = brain_region, values_from = fa_ci) %>%
  select(any_of(c("participant", "histotype", "mutation", region_list))) %>%
  arrange(participant, mutation)

lob_mask <- table_base %>%
  transmute(participant, histotype = NA_character_, mutation, brain_region,
            detected_above_LoB) %>%
  pivot_wider(names_from = brain_region, values_from = detected_above_LoB) %>%
  select(any_of(c("participant", "histotype", "mutation", region_list))) %>%
  arrange(participant, mutation) %>%
  mutate(across(all_of(intersect(region_list, names(.))), as.logical))

readr::write_csv(mytable, file.path(output_dir, "ddPCR_results_by_region_corrected.csv"))
readr::write_csv(lob_mask, file.path(output_dir, "ddPCR_results_by_region_mask_corrected.csv"))
write_html_table(
  mytable,
  lob_mask,
  file.path(output_dir, "ddPCR_results_by_region_corrected.html")
)

pooled_rows <- df_raw %>%
  mutate(is_pooled = as_true_flag(is_pooled)) %>%
  filter(is_pooled) %>%
  arrange(mutation, participant, brain_region)

readr::write_csv(
  pooled_rows,
  file.path(output_dir, "sample_region_pooled_rows_corrected.csv")
)
