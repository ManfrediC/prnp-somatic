library(readr)
library(tidyverse)
library(openxlsx)
library(magrittr)
library(binom)

# -------------------------------------
# reproducible, repo-relative paths
# -------------------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
if (length(script_path) == 0) {
  stop("Could not determine script path. Please run with: Rscript src/ddPCR/create_snv_dataframe.R")
}

script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

input_dir <- file.path(project_root, "ddPCR")
sample_details_path <- file.path(input_dir, "sample_details.xlsx")
output_dir <- file.path(project_root, "results", "ddPCR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(sample_details_path)) {
  stop("Missing metadata file: ", sample_details_path)
}

# -------------------------------------
# import CSVs exported from QuantaSoft
# -------------------------------------

mutation.list <- c("D178N", "E200K", "P102L")

# list of CSVs to import
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop("No CSV files found in ", input_dir)
}

# import
bigdata <- purrr::map_dfr(files, function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  # date from first 10 chars (YYYY-MM-DD)
  df$Date <- as.Date(substr(basename(f), 1, 10), format = "%Y-%m-%d")
  
  # experiment type from filename: SNV_<MUT>; tolerate extra suffixes (_v2, _repA) before .csv
  # whitelist only your three valid targets; case-insensitive
  mut_str <- stringr::str_extract(basename(f), "(?i)(?<=SNV_)(D178N|E200K|P102L)")
  
  if (is.na(mut_str)) {
    warning("Could not parse ExperimentType from filename: ", f)
  }
  df$ExperimentType <- toupper(mut_str)  # normalise just in case
  
  # tidy front columns
  df <- dplyr::select(df, Sample, Date, Well, Target, ExperimentType, dplyr::everything())
  df
})

# -------------------------------------
# subset mutant ddPCR probe experiments
# -------------------------------------

#columns we need
data.subset <- bigdata %>%
  dplyr::select(Sample, Date, Well, ExperimentType, Target,
                `Accepted Droplets`, Positives, Negatives,
                `Fractional Abundance`,
                PoissonFractionalAbundanceMax,
                PoissonFractionalAbundanceMin) %>%
    rename(FractionalAbundance = `Fractional Abundance`) %>%
    rename(AcceptedDroplets    = `Accepted Droplets`)

#remove unwanted suffixes from Target column
data.subset$Target <- str_replace_all(data.subset$Target, "-mut|_FAM1|_VIC2", "")

# Target == PRNP refers to WT
data.subset$Target[data.subset$Target == "PRNP"] <- "WT"

# subset ddPCR data with mutant probes (WT probes are analysed in the QC section)
data.mut <- subset(data.subset, Target %in% mutation.list)

# -------------------------------------
# clean data for mutant ddPCR probe experiments
# -------------------------------------

#clearly label Non-temple control (NTC), Mutant control (Mut) and WT samples
data.mut$Sample[grepl("NTC", data.mut$Sample, ignore.case = TRUE)] <- "NTC"
data.mut$Sample[grepl("Mut", data.mut$Sample, ignore.case = TRUE)] <- "mutant_control"
data.mut$Sample[grepl("WT", data.mut$Sample, ignore.case = TRUE)] <- "WT_control"

# in some mutant controls, the "Mut" string is missing in the Sample name
mut_missing  <- c("E200K", "P102L", "D178N", "D178") #names of Samples where "mut" is missing
data.mut$Sample[data.mut$Sample %in% mut_missing] <- "mutant_control"

#remove unwanted characters in Sample column
data.mut$Sample <- data.mut$Sample %>%
  str_replace_all("CJD-|CJD|D178N|D178|E200K|P102L", "") %>%
  str_replace("_$", "") #remove final underscore

#remove messy sample
data.mut <- subset(data.mut, Sample != "17-31_hc mixed with 17-32_cb")

# -------------------------------------
# split into controls vs analysis set
# -------------------------------------

control_samples <- c("NTC", "mutant_control", "WT_control")
data.mut.controls <- subset(data.mut, Sample %in% control_samples)
data.mut          <- subset(data.mut, !(Sample %in% control_samples))

# -------------------------------------
# harmonise sample names
# -------------------------------------

#some Sample strings end in "-bg" instead of "_bg" etc, replace hyphen (and "pons")
#also change Bologna sample names to our naming system
data.mut$Sample %<>%
  gsub("-bg", "_bg", .) %>%
  gsub("-cb", "_cb", .) %>%
  gsub("-fr", "_fr", .) %>%
  gsub("-hc", "_hc", .) %>%
  gsub("-sn", "_sn", .) %>%
  gsub("-th", "_th", .) %>%
  gsub("-pons", "_ps", .) %>%
  gsub("_pons", "_ps", .) %>%
  gsub("-cau", "_bg", .) %>%
  gsub("-ce", "_cb", .) %>%
  gsub("-fc", "_fr", .) %>%
  gsub("-hip", "_hc", .) %>%
  gsub("-mdb", "_sn", .)

# -------------------------------------
# if concentration is 0, quantasoft gives NA for fract.abund. -> replace with 0 to allow scaling
# -------------------------------------
fa_cols <- c("FractionalAbundance",
             "PoissonFractionalAbundanceMax",
             "PoissonFractionalAbundanceMin")

na_rows <- rowSums(is.na(data.mut[fa_cols])) > 0
na_data <- data.mut[na_rows, ]

# replace NA with 0 only in the FA columns
data.mut[fa_cols] <- lapply(data.mut[fa_cols], replace_na, 0)

# --------------------------------------------------------
# pool droplets (only affects samples with multiple runs)
# --------------------------------------------------------

# remove columns that aren't needed for pooling
data.mut.brief <- data.mut %>%
  select(-Well, -Target, -Negatives)

# --------------------------------------------------------
# Summarise droplet counts per Sample × Assay
# --------------------------------------------------------

counts <- data.mut.brief %>%
  group_by(Sample, ExperimentType) %>%
  summarise(
    Pos_pool   = sum(Positives,         na.rm = TRUE),
    Drop_pool  = sum(AcceptedDroplets,  na.rm = TRUE),
    n_wells    = n(),                       # how many rows pooled
    Date       = first(Date),
    .groups = "drop"
  )

# flag whether more than one well/run was combined
counts$pooled <- counts$n_wells > 1

# --------------------------------------------------------
# Compute pooled FA and its binomial CI (Clopper-Pearson)
# --------------------------------------------------------
# (these will be used only when pooled == TRUE)

counts$FA_pool <- counts$Pos_pool / counts$Drop_pool

ci <- binom.confint(
  x = counts$Pos_pool,
  n = counts$Drop_pool,
  methods = "exact"
)

counts$FA_pool_lower <- ci$lower
counts$FA_pool_upper <- ci$upper

# --------------------------------------------------------
# Combine pooled and non-pooled
# --------------------------------------------------------

# Identify Sample × Assay combinations with only ONE well
single_keys <- counts %>%                       # from step 1
  filter(!pooled) %>%                           # keeps only rows where pooled == FALSE (i.e. n_rows == 1)
  select(Sample, ExperimentType)

# Pull the Quantasoft FA and CI only for those singletons
qs <- data.mut.brief %>%
  semi_join(single_keys, by = c("Sample","ExperimentType")) %>%   # keep singles
  select(
    Sample, ExperimentType,
    FA_QS    = FractionalAbundance,
    QS_lower = PoissonFractionalAbundanceMin,
    QS_upper = PoissonFractionalAbundanceMax
  )

# Merge the two tables
merged <- left_join(counts, qs, by = c("Sample","ExperimentType"))

# --------------------------------------------------------
# >>> >>> LIMIT OF BLANK (LoB) SECTION <<< <<<
# - Use WT + NTC controls from the same plate (Date) and assay
# - Conservative p0: CP upper 95% bound on pooled blank proportion
# - Fallback: assay-wide p0 if a plate lacks blanks
# --------------------------------------------------------

# 1) Build blank table (QC: ≥10,000 droplets) from your controls object
blanks <- data.mut.controls %>%
  filter(Target %in% mutation.list, Sample %in% c("WT_control","NTC")) %>%
  filter(AcceptedDroplets >= 10000) %>%
  transmute(plate = Date,
            assay = ExperimentType,
            n = AcceptedDroplets,
            x = Positives)

# 2) Per-plate × assay pooled blanks
blank_pooled <- blanks %>%
  group_by(plate, assay) %>%
  summarise(x_blank = sum(x, na.rm = TRUE),
            n_blank = sum(n, na.rm = TRUE),
            n_wells_blank = n(),
            .groups = "drop") %>%
  mutate(p0_upper = binom.confint(x_blank, n_blank, methods = "exact")$upper)

# 3) Assay-wide fallback p0
assay_fallback <- blank_pooled %>%
  group_by(assay) %>%
  summarise(x_blank = sum(x_blank), n_blank = sum(n_blank), .groups = "drop") %>%
  mutate(p0_upper_fallback = binom.confint(x_blank, n_blank, methods = "exact")$upper) %>%
  select(assay, p0_upper_fallback)

# 4) Attach plate and assay to each sample’s pooled counts, then compute LoB
counts_lob <- merged %>%
  rename(assay = ExperimentType, n_tot = Drop_pool, x_mut = Pos_pool, plate = Date) %>%
  left_join(blank_pooled, by = c("plate","assay")) %>%
  left_join(assay_fallback, by = "assay") %>%
  mutate(p0_use = dplyr::coalesce(p0_upper, p0_upper_fallback)) %>%
  # final guard if coalesce still NA (shouldn’t happen if we had any blanks)
  mutate(p0_use = ifelse(is.na(p0_use), 0, p0_use)) %>%
  rowwise() %>%
  mutate(
    LoB_count = qbinom(0.95, size = n_tot, prob = p0_use),
    LoB_FA    = ifelse(n_tot > 0, LoB_count / n_tot, NA_real_),
    detected_LoB  = x_mut > LoB_count
  ) %>%
  ungroup() %>%
  select(Sample, assay, LoB_count, LoB_FA, detected_LoB)


# --------------------------------------------------------
# Choose the appropriate FA estimate & CI and add LoB outputs
# --------------------------------------------------------

merged <- merged %>%
  mutate(
    FA_estimate = ifelse(pooled, FA_pool, FA_QS),
    CI_lower    = ifelse(pooled, FA_pool_lower, QS_lower),
    CI_upper    = ifelse(pooled, FA_pool_upper, QS_upper)
  ) %>%
  left_join(counts_lob, by = c("Sample" = "Sample", "ExperimentType" = "assay"))

# --------------------------------------------------------
# Passed Limit of Detection?
# --------------------------------------------------------

lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)

final <- merged %>%
  mutate(
    detected_LoD = CI_lower > lod_cut[ExperimentType]
  )

# --------------------------------------------------------
# Keep what you need for the dot-plot (+ LoB fields)
# --------------------------------------------------------

ddpcr.plot.data <- final %>%
  select(Sample, ExperimentType, pooled,
         Pos_pool, Drop_pool,
         FA_estimate, CI_lower, CI_upper,
         LoB_count, LoB_FA, 
         detected_LoB, detected_LoD) %>%
  rename(
    sample_id            = Sample,
    mutation             = ExperimentType,
    is_pooled            = pooled,
    n_mut_droplets       = Pos_pool,
    n_total_droplets     = Drop_pool,
    fractional_abundance = FA_estimate,
    ci_low               = CI_lower,
    ci_high              = CI_upper,
    lob_count            = LoB_count,
    lob_fa               = LoB_FA,
    detected_above_LoB   = detected_LoB,
    detected_above_LoD   = detected_LoD
  )

# --------------------------------------------------------
# format dataframe, add participant metadata
# --------------------------------------------------------

ddpcr.plot.data <- ddpcr.plot.data %>%
  separate(sample_id, into = c("code", "brain_region"), sep = "_",
           fill = "right", remove = FALSE)

patient.names <- read.xlsx(sample_details_path) %>%
  as_tibble() %>%
  mutate(histotype = str_squish(histotype)) %>%
  distinct(code, .keep_all = TRUE) %>%
  rename(participant = new_name)

patient.names$histotype[is.na(patient.names$histotype)] <- "control"

ddpcr.plot.data <- ddpcr.plot.data %>%
  left_join(patient.names, by = "code") %>%
  select(participant, group, histotype, code, brain_region, mutation,
         is_pooled, n_mut_droplets, n_total_droplets,
         fractional_abundance, ci_low, ci_high,
         lob_count, lob_fa, 
         detected_above_LoB, detected_above_LoD)

# --------------------------------------------------------
# export
# --------------------------------------------------------

# dataset
write.xlsx(ddpcr.plot.data, file.path(output_dir, "SNV_data_final.xlsx"))

# assay-wide fallback p0
write.csv(assay_fallback, file.path(output_dir, "p0_fallback.csv"), row.names = FALSE)

###############################################################################################################

# --------------------------------------------------------
# NEW DATAFRAME - POOLED BY PARTICIPANT
# --------------------------------------------------------

# fix column names
merged.new <- merged %>%
  rename(assay = ExperimentType) %>% # rename assay
  mutate(code = sub("_.*$", "", Sample)) %>% # divide Sample to retrieve patient code
  relocate(code, .after = Sample)

# droplet sums
pooled_participant_assay <- merged.new %>%
  group_by(code, assay) %>%
  summarise(
    pos_total  = sum(Pos_pool,  na.rm = TRUE),
    drop_total = sum(Drop_pool, na.rm = TRUE),
    .groups = "drop"
  )

# drop 14-2 E200K as this sample was found to be heterozygous mutant E200K
pooled_participant_assay <- pooled_participant_assay %>%
  filter(!(code == "14-2" & assay == "E200K"))

# compute fractional abundance
pooled_participant_assay$FA <- pooled_participant_assay$pos_total / pooled_participant_assay$drop_total

# compute confidence intervals
ci_pooled <- binom.confint(
  x = pooled_participant_assay$pos_total,
  n = pooled_participant_assay$drop_total,
  methods = "exact"
)

pooled_participant_assay$ci_lower <- ci_pooled$lower
pooled_participant_assay$ci_upper <- ci_pooled$upper

# --------------------------------------------------------
# create mapping dataframe of plates x assays with respective p0 
# --------------------------------------------------------

# plate mapping (we need to know which plates the droplets come from)
plate_mapping <- merged.new %>%
  select(Sample, code, plate = Date, assay, Drop_pool)

# sum droplets (so we know how many came from each plate)
plate_mapping <- plate_mapping %>%
  group_by(code, assay, plate) %>%
  summarise(
    drop_from_plate  = sum(Drop_pool, na.rm = TRUE),
    .groups = "drop"
  )

# join plate-specific blank rates
plate_mapping <- plate_mapping %>%
  left_join(blank_pooled %>% select(plate, assay, p0_upper), # join with columns we need
            by = c("assay", "plate"))

# For plates without blanks, use the assay-wide p0_upper_fallback.
plate_mapping <- plate_mapping %>%
  left_join(assay_fallback, by = "assay") %>%
  mutate(p0_per_plate = coalesce(p0_upper, p0_upper_fallback)) # if p0_upper is NA, use fallback

# --------------------------------------------------------
# LoB for code x assay with conservative max-of-plates approach
# --------------------------------------------------------

# one row per code × assay with the conservative p0
p0_per_participant <- plate_mapping %>%
  group_by(code, assay) %>%
  slice_max(p0_per_plate, n = 1, with_ties = FALSE) %>%  # plate that set the threshold
  summarise(p0_max = first(p0_per_plate),
            plate_p0_max = first(plate),
            .groups = "drop")

# Join to the pooled counts (one row per code × assay) and calculate LoB
pooled_participant_assay <- pooled_participant_assay %>%
  left_join(p0_per_participant, 
            by = c("code", "assay")) %>%
  mutate(
    LoB_count          = qbinom(0.95, size = drop_total, prob = p0_max),
    LoB_FA             = LoB_count / drop_total,
    detected_above_lob = pos_total > LoB_count
  )

# --------------------------------------------------------
# detected above LoD?
# --------------------------------------------------------

pooled_participant_assay <- pooled_participant_assay %>%
  mutate(
    detected_above_lod = ci_lower > lod_cut[assay]
  )

# --------------------------------------------------------
# participant metadata, export
# --------------------------------------------------------

pooled_plot_data <- pooled_participant_assay %>%
  left_join(patient.names, by = "code") %>%
  transmute(
    participant,
    code,
    group,
    histotype,
    brain_region        = "pooled",        # or NA_character_ if you prefer
    mutation            = assay,           # if you already have a 'mutation' col, drop/adjust this
    n_mut_droplets      = pos_total,
    n_total_droplets    = drop_total,
    fractional_abundance = FA,
    ci_low              = ci_lower,
    ci_high             = ci_upper,
    lob_count           = LoB_count,
    lob_fa              = LoB_FA,
    detected_above_lob,
    detected_above_lod
  )

# export
write.xlsx(pooled_plot_data, file.path(output_dir, "SNV_pooled_participant.xlsx"))



