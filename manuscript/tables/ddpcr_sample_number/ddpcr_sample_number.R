library(tidyverse)
library(readxl)

# ---- user inputs ----
xlsx_path <- "results/ddPCR/SNV_data_final.xlsx"


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

# ---- read + basic filtering ----
ddpcr <- read_excel(xlsx_path, sheet = sheet_name) %>%
  mutate(
    participant = as.character(participant),
    group = as.character(group),
    mutation = as.character(mutation),
    brain_region = as.character(brain_region)
  )


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
  mutate(
    name = participant
  ) %>%
  rowwise() %>%
  mutate(
    max_maf_overall = safe_max(c(max_maf_D178N, max_maf_E200K, max_maf_P102L))
  ) %>%
  ungroup() %>%
  select(
    name,
    analysed_samples,
    analysed_D178N, analysed_E200K, analysed_P102L,
    max_maf_D178N, max_maf_E200K, max_maf_P102L,
    max_maf_overall,
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
  select(
    name,
    analysed_samples,
    analysed_D178N, analysed_E200K, analysed_P102L,
    max_maf_D178N, max_maf_E200K, max_maf_P102L,
    max_maf_overall,
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

kable(ddpcr_summary,
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      align = "lrrrrrrrr")



# ---- export as CSV ----
setwd("manuscript/tables/ddpcr_sample_number")

write.csv(ddpcr_summary, "ddPCR_sample_number.csv", row.names = FALSE)
