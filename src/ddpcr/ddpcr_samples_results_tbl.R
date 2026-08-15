library(tidyverse)
library(openxlsx)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_path <- file.path(project_root, "results", "ddPCR", "SNV_data_final.xlsx")
previous_results_dir <- file.path(project_root, "legacy", "results", "ddpcr", "2026-06-finalisation", "previous_official")
current_input_path <- file.path(previous_results_dir, "SNV_data_final.xlsx")
default_output_dir <- file.path(project_root, "manuscript", "tables", "ddpcr_sample_results")
output_dir <- Sys.getenv("PRNP_SAMPLE_RESULTS_OUTPUT_DIR", unset = default_output_dir)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)
region_labels <- c(
  bg = "basal ganglia",
  bs = "brainstem",
  cb = "cerebellum",
  fr = "frontal cortex",
  hc = "hippocampus",
  th = "thalamus"
)
region_list <- unname(region_labels[c("bg", "bs", "cb", "hc", "fr", "th")])

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

write_region_xlsx <- function(table_data, mask_data, output_path) {
  wb <- openxlsx::createWorkbook()
  ws <- openxlsx::addWorksheet(wb, "ddPCR by region")
  openxlsx::writeData(wb, ws, table_data, headerStyle = openxlsx::createStyle(textDecoration = "bold"))
  grey_style <- openxlsx::createStyle(fgFill = "#e5e7eb")
  for (col in names(mask_data)) {
    if (!col %in% names(table_data)) next
    col_idx <- which(names(table_data) == col)
    masked_rows <- which(as_true_flag(mask_data[[col]]))
    if (length(masked_rows) > 0L) {
      openxlsx::addStyle(
        wb, ws, style = grey_style,
        rows = masked_rows + 1L,
        cols = rep.int(col_idx, length(masked_rows)),
        gridExpand = FALSE
      )
    }
  }
  openxlsx::freezePane(wb, ws, firstActiveRow = 2L, firstActiveCol = 1L)
  openxlsx::setColWidths(wb, ws, cols = seq_along(table_data), widths = "auto")
  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
}

escape_latex <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- gsub("\\", "\\textbackslash{}", x, fixed = TRUE)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("$", "\\$", x, fixed = TRUE)
  x <- gsub("#", "\\#", x, fixed = TRUE)
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x <- gsub("{", "\\{", x, fixed = TRUE)
  x <- gsub("}", "\\}", x, fixed = TRUE)
  x
}

write_region_tex_table <- function(table_data, mask_data, output_path, caption, label, compact = FALSE) {
  region_cols <- setdiff(names(table_data), c("participant", "histotype", "mutation"))
  region_headers <- c(
    `basal ganglia` = "Basal ganglia (VAF, 95\\% CI)",
    brainstem = "Brainstem (VAF, 95\\% CI)",
    cerebellum = "Cerebellum (VAF, 95\\% CI)",
    hippocampus = "Hippocampus (VAF, 95\\% CI)",
    `frontal cortex` = "Frontal cortex (VAF, 95\\% CI)",
    thalamus = "Thalamus (VAF, 95\\% CI)"
  )
  header <- c(
    "Participant",
    "Histotype",
    "Mutation",
    unname(region_headers[region_cols])
  )
  col_spec <- paste(c("l", "l", "l", rep("r", length(region_cols))), collapse = "")

  rows <- vapply(seq_len(nrow(table_data)), function(i) {
    cells <- c(
      escape_latex(table_data$participant[i]),
      escape_latex(table_data$histotype[i]),
      escape_latex(table_data$mutation[i])
    )
    region_cells <- vapply(region_cols, function(col) {
      value <- escape_latex(table_data[[col]][i])
      if (col %in% names(mask_data) && isTRUE(mask_data[[col]][i])) {
        paste0("\\cellcolor{gray!10}{", value, "}")
      } else {
        value
      }
    }, character(1))
    paste0(paste(c(cells, region_cells), collapse = " & "), "\\\\")
  }, character(1))

  lines <- c(
    "% Auto-generated by src/ddpcr/ddpcr_samples_results_tbl.R",
    "\\begin{table}[p]",
    "\\centering",
    paste0("\\caption{\\label{", label, "}", caption, "}"),
    if (compact) "\\fontsize{8}{10}\\selectfont" else character(),
    "\\resizebox{\\linewidth}{!}{%",
    paste0("\\begin{tabular}{", col_spec, "}"),
    "\\toprule",
    paste0(paste(header, collapse = " & "), "\\\\"),
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\end{table}",
    "",
    "\\clearpage"
  )
  writeLines(lines, output_path, useBytes = TRUE)
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
  file.path(output_dir, "ddPCR_results_table.xlsx"),
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
    file.path(output_dir, "sample_region_pooled_and_lod_changed_rows.csv")
  )
}

table_base <- df %>%
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

readr::write_csv(mytable, file.path(output_dir, "ddPCR_results_by_region.csv"))
readr::write_csv(lob_mask, file.path(output_dir, "ddPCR_results_by_region_mask.csv"))
write_html_table(
  mytable,
  lob_mask,
  file.path(output_dir, "ddPCR_results_by_region.html")
)

write_region_tex_table(
  mytable,
  lob_mask,
  file.path(output_dir, "ddpcr_results_by_region_table.tex"),
  "ddPCR results by brain region. Shaded cells show VAF above the limit of blank.",
  "tab:ddpcr_results_by_region"
)

supplement_xlsx_path <- file.path(
  project_root, "manuscript", "tables", "supplement",
  "ddPCR_results_by_region_S3.xlsx"
)
write_region_xlsx(mytable, lob_mask, supplement_xlsx_path)

write_region_tex_table(
  mytable,
  lob_mask,
  file.path(output_dir, "ddpcr_results_table.tex"),
  "ddPCR sample-region results. Shaded cells show VAF above the limit of blank.",
  "tab:ddpcr_results_complete",
  compact = TRUE
)

pooled_rows <- df_raw %>%
  mutate(is_pooled = as_true_flag(is_pooled)) %>%
  filter(is_pooled) %>%
  arrange(mutation, participant, brain_region)

readr::write_csv(
  pooled_rows,
  file.path(output_dir, "sample_region_pooled_rows.csv")
)

generated_outputs <- c(
  file.path(output_dir, "ddPCR_results_by_region.csv"),
  file.path(output_dir, "ddPCR_results_by_region_mask.csv"),
  file.path(output_dir, "ddPCR_results_by_region.html"),
  file.path(output_dir, "ddpcr_results_by_region_table.tex"),
  file.path(output_dir, "ddpcr_results_table.tex"),
  file.path(output_dir, "sample_region_pooled_rows.csv"),
  supplement_xlsx_path
)

missing_or_empty <- generated_outputs[
  !file.exists(generated_outputs) |
    is.na(file.info(generated_outputs)$size) |
    file.info(generated_outputs)$size <= 0
]
if (length(missing_or_empty) > 0L) {
  stop("Regional ddPCR output is missing or empty: ", paste(missing_or_empty, collapse = ", "))
}
