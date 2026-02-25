# ================================
# LoD tables export (no plotting)
# ================================
# Outputs, per mutation:
#   <MUT>_lob_table.csv
#   <MUT>_fraction_CI.csv
# And combined across all mutations:
#   ALL_mutations_lob_table.csv
#   ALL_mutations_fraction_CI.csv

library(tidyverse)

# ---- paths ----
input_dir <- "manuscript/legacy/figures/ddPCR_LoD"

mutation.list <- c("D178N", "E200K", "P102L")

# helper to get ordered sample IDs for a mutation
make_sample_list <- function(mut) {
  c(paste0(mut, "_NTC"),
    paste0(mut, "_WT"),
    paste0(mut, "_0.05"),
    paste0(mut, "_0.1"),
    paste0(mut, "_0.33"),
    paste0(mut, "_1"))
}

# core worker: returns a list of two tibbles (LoB table, FA CI table) with a 'mutation' column
process_mutation <- function(mut) {
  infile <- paste0(mut, "_LOD.csv")
  dat <- read.csv(file.path(input_dir, infile), fileEncoding = "UTF-8-BOM") %>% as_tibble()
  
  # set sample order to keep outputs tidy/consistent
  samp_levels <- make_sample_list(mut)
  dat$Sample <- factor(dat$Sample, levels = samp_levels)
  
  # ---- pooled p0 from blanks on the mutation channel ----
  blanks <- dat %>%
    filter(Sample %in% c(paste0(mut, "_NTC"), paste0(mut, "_WT")),
           Target == mut)
  
  p0 <- sum(blanks$Positives) / sum(blanks$Accepted.Droplets)
  
  # per-blank LoB and LoBFA; keep a compact table
  lob_tbl <- blanks %>%
    rowwise() %>%
    mutate(
      LoB_count = qbinom(0.95, size = Accepted.Droplets, prob = p0),
      LoBFA     = LoB_count / Accepted.Droplets
    ) %>%
    ungroup() %>%
    select(Well, Sample, Accepted.Droplets, Positives, LoB_count, LoBFA) %>%
    mutate(mutation = mut, .before = 1)
  
  # per-sample FA with Poisson CI (mutation channel only)
  fa_ci_tbl <- dat %>%
    filter(Target == mut) %>%
    select(Sample,
           FractionalAbundance,
           CI_low  = PoissonFractionalAbundanceMin,
           CI_high = PoissonFractionalAbundanceMax) %>%
    mutate(mutation = mut, .before = 1)
  
  # write per-mutation CSVs
  write.csv(lob_tbl, file = paste0(mut, "_lob_table.csv"), row.names = FALSE)
  write.csv(fa_ci_tbl, file = paste0(mut, "_fraction_CI.csv"), row.names = FALSE)
  
  list(lob = lob_tbl, fa = fa_ci_tbl)
}

# ---- run for all mutations and also write combined summaries ----
all_results <- lapply(mutation.list, process_mutation)

all_lob <- bind_rows(lapply(all_results, `[[`, "lob"))
all_fa  <- bind_rows(lapply(all_results, `[[`, "fa"))

# export to tables directory

setwd("manuscript/tables/lod_calculations")

write.csv(all_lob, "ALL_mutations_lob_table.csv", row.names = FALSE)
write.csv(all_fa,  "ALL_mutations_fraction_CI.csv", row.names = FALSE)
