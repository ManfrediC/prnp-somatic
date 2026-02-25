# ddPCR results table for all the samples
# Manfredi Carta - 20 Oct 2025

library(readr)
library(tidyverse)
library(openxlsx)
library(gtools)
library(glue)
library(gt)

# import
setwd("manuscript/tables/ddpcr_sample_results")
df <- read.xlsx("SNV_data_final.xlsx")

# -------------------------------------
# improve legibility
# -------------------------------------

df <- df %>%
  mutate(
    brain_region = recode(
      brain_region,
      bg = "basal ganglia",
      cb = "cerebellum",
      fr = "frontal cortex",
      hc = "hippocampus",
      ps = "pons",
      sn = "substantia nigra",
      th = "thalamus"
    )
  )

# reduce columns
df <- df %>%
  select(-group, -code, -is_pooled, -lob_count)

# LoB and LoD as yes/no
df$detected_above_LoB <- ifelse(df$detected_above_LoB, "Yes", "No")
df$detected_above_LoD <- ifelse(df$detected_above_LoD, "Yes", "No")

# better colname for LoB
df <- df %>%
  rename(LoB = lob_fa)

# add LoD column
lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)

df <- df %>%
  mutate(
    LoD = lod_cut[mutation]
  ) %>%
  relocate(LoD, .before = detected_above_LoD)

# -------------------------------------
# summarised Excel for the supplement
# -------------------------------------

setwd("manuscript/tables/ddpcr_sample_results")
write.xlsx(df, "ddPCR_results_table.xlsx")


# -------------------------------------
# which samples surpassed LoD and how many? Remove CJD30 as it's heterozygous for E200K
# -------------------------------------

# number of samples assayed


# how many over LoD?

D178N_LoD <- df[df$mutation == "D178N" & df$fractional_abundance > lod_cut[["D178N"]],]
nrow(D178N_LoD)

E200K_LoD <- df[df$mutation == "E200K" & df$fractional_abundance > lod_cut[["E200K"]],]
E200K_LoD <- E200K_LoD[E200K_LoD$participant != "CJD30",]

P102L_LoD <- df[df$mutation == "P102L" & df$fractional_abundance > lod_cut[["P102L"]],]


mean_FA_above_LoD <- bind_rows(D178N_LoD, E200K_LoD, P102L_LoD)

paste0("In D178N, ", nrow(D178N_LoD)," samples had mean FA that surpassed LoD, but CI did not")
paste0("In E200K, ", nrow(E200K_LoD)," samples had mean FA that surpassed LoD, but CI did not")
paste0("In P102L, ", nrow(P102L_LoD)," samples had mean FA that surpassed LoD, but CI did not")


# -------------------bind_rows()# ----------------------- --------------
# table for the paper + mask to be used for shading
# -------------------------------------

# remove ps (pons) as there are only 3 samples
df <- df[df$brain_region != "pons",]

# round all numeric to 3 decimal places
df <- df %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# condense FA and confidence intervals to one column
df <- df %>% mutate(
  fa_ci = paste0(fractional_abundance, " (", ci_low, "-", ci_high, ")")
  ) %>%
  select(-fractional_abundance, -ci_low, -ci_high) %>%
  relocate(fa_ci, .after = mutation)

# remove less interesting columns
mytable <- df %>% 
  select(-n_mut_droplets, -n_total_droplets, -LoB, -LoD, -detected_above_LoB, -detected_above_LoD)

# pivot brain regions wider
mytable <- mytable %>%
  mutate(brain_region = factor(brain_region, levels = sort(unique(brain_region)))) %>%
  pivot_wider(names_from = brain_region, values_from = fa_ci)

# arrange columns

region_list <- c(
  "basal ganglia", "cerebellum", "hippocampus",
  "frontal cortex", "substantia nigra", "thalamus"
)
col_order <- c("participant", "histotype", "mutation", region_list)

mytable <- mytable %>%
  select(any_of(col_order))

# sort samples alphabetically
mytable <- mytable %>%
  mutate(participant = factor(participant, levels = mixedsort(unique(participant)))) %>%
  arrange(participant)

# LoB df for highlighting values above LoB

lob_df <- df %>% select(
  participant, brain_region, mutation, detected_above_LoB
)

# LoB mask with same orientation for shading cells
lob_mask <- lob_df %>%
  pivot_wider(names_from = brain_region, values_from = detected_above_LoB) %>%
  select(any_of(col_order)) %>% # order columns
  #order alphabetically
  mutate(participant = factor(participant, levels = mixedsort(unique(participant)))) %>%
  arrange(participant)

lob_mask[region_list] <- lapply(lob_mask[region_list], function(x) {
  ifelse(x == "Yes", TRUE,
         ifelse(x == "No", FALSE, NA))
})

# mask must have same columns as data
lob_mask$histotype <- NA

lob_mask <- lob_mask %>% 
  relocate(histotype, .after = participant)

glimpse(mytable)
glimpse(lob_df)



# -------------------------------------
# export table and mask
# -------------------------------------

setwd("manuscript/tables/ddpcr_sample_results")

write.csv(mytable, "ddPCR_results_by_region.csv", row.names = FALSE)
write.csv(lob_mask, "ddPCR_results_by_region_mask.csv", row.names = FALSE)

# -------------------------------------
# create Latex code
# -------------------------------------


# 1) Identify value columns
id_cols <- c("participant", "histotype", "mutation")
value_cols <- setdiff(names(mytable), id_cols)

# 2) Get coordinates to shade (where lob_mask == TRUE in value columns)
shade_coords <- lob_mask %>%
  mutate(.row = row_number()) %>%
  pivot_longer(
    cols = all_of(value_cols),
    names_to = ".col",
    values_to = ".mask"
  ) %>%
  filter(.mask %in% TRUE) %>%
  select(.row, .col)

# 3) Build gt table (you can tweak fonts/size later in LaTeX if desired)
gt_tbl <- mytable %>%
  gt() %>%
  fmt_missing(columns = everything(), missing_text = "—")

# 4) Apply grey shading to the selected cell coordinates
#    Using a simple loop for clarity and reliability
for (i in seq_len(nrow(shade_coords))) {
  gt_tbl <- gt_tbl %>%
    tab_style(
      style = cell_fill(color = "gray90"),
      locations = cells_body(
        columns = all_of(shade_coords$.col[i]),
        rows    = shade_coords$.row[i]
      )
    )
}

# 5) (Option A) Get LaTeX code as a string (to paste into your .tex)
latex_code <- gt::as_latex(gt_tbl)
cat(as.character(latex_code))
