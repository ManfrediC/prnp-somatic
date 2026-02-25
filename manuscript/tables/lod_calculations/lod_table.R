# ddPCR LoD experiment results table
# Manfredi Carta - 20 Oct 2025

library(readr)
library(tidyverse)
library(openxlsx)
library(knitr)
library(kableExtra)

# import results of LoD ddPCR
input_dir <- "manuscript/legacy/figures/ddPCR_LoD"

D178N <- read.csv(file.path(input_dir, "D178N_LOD.csv"))
E200K <- read.csv(file.path(input_dir, "E200K_LOD.csv"))
P102L <- read.csv(file.path(input_dir, "P102L_LOD.csv"))

# collate to one dataframe
lod_data <- rbind(D178N, E200K, P102L)

# keep required columns
lod_data <- lod_data %>%
  select(Target, Sample,
         Conc, TotalConfMax, TotalConfMin,
         FractionalAbundance, PoissonFractionalAbundanceMax,
         PoissonFractionalAbundanceMin,
         Accepted.Droplets, Positives, Negatives)

# remove assays targeting WT
lod_data <- lod_data[lod_data$Target != "WT",]

# harmonise WT, NTC and P102L sample names
lod_data$Sample[grepl("WT", lod_data$Sample)] <- "WT"
lod_data$Sample[grepl("NTC", lod_data$Sample)] <- "NTC"
lod_data$Sample <- gsub("P102L-mut", "P102L", lod_data$Sample)

# extract expected percentage
lod_data <- lod_data %>%
  mutate(expected_pct = str_extract(Sample, "(?<=_)[\\d.]+$")) %>%
  relocate(expected_pct, .after = Target)

# set expected_pct to 0 for NTC and WT
lod_data$expected_pct[lod_data$Sample %in% c("NTC", "WT")] <- 0

# set FA to 0 if it's NA (NA results if there are no droplets)
lod_data$FractionalAbundance[is.na(lod_data$FractionalAbundance) & lod_data$Positives == 0] <- 0
lod_data$PoissonFractionalAbundanceMax[is.na(lod_data$PoissonFractionalAbundanceMax) & lod_data$Positives == 0] <- 0
lod_data$PoissonFractionalAbundanceMin[is.na(lod_data$PoissonFractionalAbundanceMin) & lod_data$Positives == 0] <- 0

# ---------------------------
# format dataframe
# ---------------------------

mutation_order <- c("D178N", "E200K", "P102L")

# round to 2 decimal places
num_col <- c("Conc", "TotalConfMax", "TotalConfMin",
             "FractionalAbundance", "PoissonFractionalAbundanceMax", "PoissonFractionalAbundanceMin",
             "Accepted.Droplets", "Positives")

lod_data <- lod_data %>%
  mutate(across(all_of(num_col), ~ round(., 2)))

# sort rows
lod_data <- lod_data %>%
  arrange(Target) %>%
  group_by(Target) %>%
  arrange(expected_pct, .by_group = TRUE) %>%
  ungroup()

# rename samples

# lod_data <- lod_data %>%
#   mutate(
#     Sample = if_else(
#       Sample %in% c("WT", "NTC"),
#       Sample,
#       str_replace(Sample, "_", " ") %>% paste0("%")
#     )
#   )

lod_data <- lod_data %>%
  mutate(
    Sample = if_else(
      Sample %in% c("WT", "NTC"),
      Sample,
      str_replace(Sample, "(D178N_|E200K_|P102L_)", " ") %>% paste0("%")
    )
  )

# keep relevant columns
lod_data <- lod_data %>%
  select("Target", "Sample", 
         #"expected_pct",
         "FractionalAbundance", "PoissonFractionalAbundanceMin", "PoissonFractionalAbundanceMax",
         "Accepted.Droplets", "Positives", "Negatives")

# combine to Measured (95% CI) column, 0.17 (0.15–0.20)

lod_data <- lod_data %>%
  mutate(
    combined = paste0(FractionalAbundance, " (", 
                      PoissonFractionalAbundanceMin, "-", 
                      PoissonFractionalAbundanceMax, ")")
  ) %>%
  relocate(combined, .after = Sample) %>%
  select(-FractionalAbundance, -PoissonFractionalAbundanceMin, -PoissonFractionalAbundanceMax)

# rename the columns

# lod_data <- lod_data %>%
#   rename(
#     Mutation = Target,
#     `Expected (%)` = expected_pct,
#     `Measured (%)` = FractionalAbundance,
#     `Measured lower CI (%)` = PoissonFractionalAbundanceMin,
#     `Measured upper CI (%)` = PoissonFractionalAbundanceMax,
#     `Accepted droplets` = Accepted.Droplets,
#     `Positive droplets` = Positives,
#     `Negative droplets` = Negatives
#   )


lod_data <- lod_data %>%
  rename(
    Mutation = Target,
    `Fractional abundance (95% CI)`= combined,
    `Accepted droplets` = Accepted.Droplets,
    `Positive droplets` = Positives,
    `Negative droplets` = Negatives
  )

# ---------------------------
# export
# ---------------------------

output_dir <- "manuscript/tables/lod_calculations"

write.csv(lod_data, file.path(output_dir, "LoD_data.csv"), row.names = FALSE)


# ---------------------------
# create Latex table
# ---------------------------

kable(
  lod_data,
  format = "latex",
  booktabs = TRUE,
  digits = 2,
  caption = "ddPCR quantification of mutant allele frequencies in dilution series and controls."
)


# maybe nicer
lod_data %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    digits = 2,
    caption = "ddPCR quantification of mutant allele frequencies in dilution series and controls."
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position", "scale_down"),
    stripe_color = "gray!10"
  )
