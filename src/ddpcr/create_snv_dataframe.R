library(readr)
library(tidyverse)
library(openxlsx)
library(magrittr)
library(binom)

# -------------------------------------
# reproducible, repo-relative paths
# -------------------------------------

get_ddpcr_project_root <- function() {
  # Convert a file or directory path into an existing directory path.
  normalise_existing_dir <- function(path) {
    if (length(path) != 1L || is.na(path) || !nzchar(path)) {
      return(character(0))
    }
    if (file.exists(path) && !dir.exists(path)) {
      path <- dirname(path)
    }
    if (!dir.exists(path)) {
      return(character(0))
    }
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }

  # Walk upwards from a starting directory until the repository markers are found.
  find_project_root <- function(start_dir) {
    current_dir <- normalise_existing_dir(start_dir)
    if (length(current_dir) == 0L) {
      return(character(0))
    }

    # The ddPCR source directory and conda environment file identify this repo.
    repeat {
      if (dir.exists(file.path(current_dir, "src", "ddpcr")) &&
          file.exists(file.path(current_dir, "env", "ddpcr.environment.yml"))) {
        return(current_dir)
      }

      # Stop once dirname() can no longer move up the directory tree.
      parent_dir <- dirname(current_dir)
      if (identical(parent_dir, current_dir)) {
        return(character(0))
      }
      current_dir <- parent_dir
    }
  }

  # When run with Rscript, prefer the path supplied by the --file argument.
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  # When sourced, sys.frames() can expose the source file path via frame$ofile.
  source_files <- unlist(lapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) {
      character(0)
    } else {
      frame$ofile
    }
  }))

  # When run interactively in RStudio, use the active editor path if available.
  rstudio_path <- character(0)
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    rstudio_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) character(0)
    )
  }

  # Try all plausible launch contexts, then de-duplicate before searching upward.
  candidate_dirs <- unique(c(
    normalise_existing_dir(sub("^--file=", "", file_arg[1])),
    unlist(lapply(source_files, normalise_existing_dir), use.names = FALSE),
    normalise_existing_dir(rstudio_path[1]),
    normalise_existing_dir(getwd())
  ))

  # Return the first candidate that resolves to this project's repository root.
  for (candidate_dir in candidate_dirs) {
    project_root <- find_project_root(candidate_dir)
    if (length(project_root) == 1L) {
      return(project_root)
    }
  }

  # Fail visibly if none of the launch contexts can be tied back to the repo.
  stop(
    "Could not determine project root. ",
    "Set the working directory to the repository root, ",
    "or source the script file directly."
  )
}

# Anchor all subsequent inputs and outputs to the detected repository root.
project_root <- get_ddpcr_project_root()

# Define the input metadata, raw database, validation, and output locations.
input_dir <- file.path(project_root, "raw", "ddpcr")
sample_details_path <- file.path(input_dir, "sample_details.xlsx")
raw_database_root <- file.path(project_root, "raw", "ddpcr")
output_dir <- file.path(project_root, "results", "ddPCR")
validation_dir <- file.path(output_dir, "validation")

# Ensure generated result and validation folders exist before writing files.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

# Express fractional abundance outputs as percentages.
fa_scale <- 100

# The sample metadata workbook is required before any ddPCR import can proceed.
if (!file.exists(sample_details_path)) {
  stop("Missing metadata file: ", sample_details_path)
}

# Load helpers that import and validate the raw ddPCR database exports.
source(file.path(project_root, "src", "ddpcr", "ddpcr_raw_import_helpers.R"))

# -------------------------------------
# import active wells from the raw ddPCR database
# -------------------------------------

mutation.list <- c("D178N", "E200K", "P102L")

# load raw ddPCR data and CSV export data from Quantasoft
bigdata <- read_ddpcr_raw_bigdata_from_database(
  raw_root = raw_database_root,
  validation_dir = validation_dir
)

# -------------------------------------
# subset mutant ddPCR probe experiments
# -------------------------------------

#columns we need
data.subset <- bigdata %>%
  dplyr::select(Sample, Date, Well, ExperimentType, Target,
                `Accepted Droplets`, Positives, Negatives,
                `Ch1+Ch2+`, `Ch1+Ch2-`, `Ch1-Ch2+`, `Ch1-Ch2-`) %>%
  rename(AcceptedDroplets = `Accepted Droplets`)

#remove unwanted suffixes from Target column
data.subset$Target <- str_replace_all(data.subset$Target, "-mut|_FAM1|_VIC2", "")

# Target == PRNP refers to WT (according to our droplet reader settings)
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

#some Sample strings end in "-bg" instead of "_bg" etc, replace hyphen
#also change Bologna brain region designators to our naming system
#midbrain, substantia nigra, and pons are all grouped as brainstem (bs)
data.mut$Sample %<>%
  gsub("-bg", "_bg", .) %>%
  gsub("-cb", "_cb", .) %>%
  gsub("-fr", "_fr", .) %>%
  gsub("-hc", "_hc", .) %>%
  gsub("-th", "_th", .) %>%
  gsub("[-_](sn|mb|mdb|pons)$", "_bs", .) %>%
  gsub("-cau", "_bg", .) %>%
  gsub("-ce", "_cb", .) %>%
  gsub("-fc", "_fr", .) %>%
  gsub("-hip", "_hc", .)

# --------------------------------------------------------
# pool droplets (only affects samples with multiple runs)
# --------------------------------------------------------

# remove columns that aren't needed for pooling
data.mut.brief <- data.mut %>%
  select(-Well, -Target, -Negatives)

data.mut.brief <- data.mut.brief %>%
  mutate(
    n_double_positive_droplets = as.numeric(`Ch1+Ch2+`),
    n_ref_only_droplets = case_when(
      ExperimentType == "P102L" ~ as.numeric(`Ch1-Ch2+`),
      TRUE ~ as.numeric(`Ch1+Ch2-`)
    ),
    n_mut_only_droplets = case_when(
      ExperimentType == "P102L" ~ as.numeric(`Ch1+Ch2-`),
      TRUE ~ as.numeric(`Ch1-Ch2+`)
    ),
    n_ref_positive_droplets = n_double_positive_droplets + n_ref_only_droplets,
    n_mut_positive_droplets = n_double_positive_droplets + n_mut_only_droplets
  )

# --------------------------------------------------------
# Summarise droplet counts per Sample × Assay
# --------------------------------------------------------

counts <- data.mut.brief %>%
  group_by(Sample, ExperimentType) %>%
  summarise(
    Pos_pool     = sum(n_mut_positive_droplets, na.rm = TRUE),
    Ref_pos_pool = sum(n_ref_positive_droplets, na.rm = TRUE),
    Drop_pool    = sum(AcceptedDroplets,        na.rm = TRUE),
    n_wells    = n(),                       # how many rows pooled
    Date       = first(Date),
    .groups = "drop"
  )

# flag whether more than one well/run was combined
counts$pooled <- counts$n_wells > 1

# --------------------------------------------------------
# calculate mutant fractional abundance from Poisson-corrected molecules
# --------------------------------------------------------

# dMIQE-style lambda is estimated from accepted positive/negative partitions.
# The CI follows the Dube/Fieller dPCR concentration-ratio approach, applied to
# lambda_mut / lambda_wt and transformed to mutant fractional abundance.
count_fa <- purrr::pmap_dfr(
  list(
    Sample = counts$Sample,
    ExperimentType = counts$ExperimentType,
    ref_positive = counts$Ref_pos_pool,
    ref_negative = counts$Drop_pool - counts$Ref_pos_pool,
    mut_positive = counts$Pos_pool,
    mut_negative = counts$Drop_pool - counts$Pos_pool,
    total = counts$Drop_pool
  ),
  function(Sample, ExperimentType, ref_positive, ref_negative, mut_positive, mut_negative, total) {
    fractional_abundance(
      ref_positive = ref_positive,
      ref_negative = ref_negative,
      mut_positive = mut_positive,
      mut_negative = mut_negative,
      total = total
    ) %>%
      mutate(Sample = Sample, ExperimentType = ExperimentType, .before = 1)
  }
)

merged <- counts %>%
  left_join(count_fa, by = c("Sample", "ExperimentType")) %>%
  rename(
    FA_estimate = fractional_abundance,
    CI_lower = ci_low,
    CI_upper = ci_high
  )

# --------------------------------------------------------
# >>> >>> LIMIT OF BLANK (LoB) SECTION <<< <<<
# - Use WT genomic controls from the same plate (Date) and assay
# - NTCs are contamination controls; they are not used to model rare SNV
#   background in wild-type genomic DNA
# - Conservative p0: CP upper 95% bound on pooled blank proportion
# - Fallback: assay-wide p0 if a plate lacks blanks
# --------------------------------------------------------

# WT controls are the biological blanks for the LoB calculation because they
# contain genomic DNA background without the targeted mutation.
lob_blank_samples <- "WT_control"

# 1) Build blank table (QC: >=10,000 droplets) from WT control wells
blanks <- data.mut.controls %>%
  filter(Target %in% mutation.list, Sample %in% lob_blank_samples) %>% # WT genomic blank wells
  filter(AcceptedDroplets >= 10000) %>% # only wells with at least 10,000 accepted droplets
  transmute(plate = Date,
            assay = ExperimentType,
            n = AcceptedDroplets,
            x = Positives)

# 2) Pool blanks by plate and assay
blank_pooled <- blanks %>%
  group_by(plate, assay) %>%
  summarise(x_blank = sum(x, na.rm = TRUE),
            n_blank = sum(n, na.rm = TRUE),
            n_wells_blank = n(),
            .groups = "drop") %>%
  mutate(p0_upper = binom.confint(x_blank, n_blank, methods = "exact")$upper)
    # p0: upper 95% Clopper-Pearson bound for the blank positive rate

# 3) Assay-wide fallback p0
assay_fallback <- blank_pooled %>%
  group_by(assay) %>%
  summarise(x_blank = sum(x_blank), n_blank = sum(n_blank), .groups = "drop") %>%
  mutate(p0_upper_fallback = binom.confint(x_blank, n_blank, methods = "exact")$upper) %>%
  select(assay, p0_upper_fallback)

# 4) For each sample x assay, use the highest blank rate from any plate that
# contributed droplets to that pooled sample. This is conservative when a
# sample was run on more than one plate.
sample_lob_p0 <- data.mut.brief %>%
  select(Sample, assay = ExperimentType, plate = Date) %>%
  distinct() %>%
  left_join(blank_pooled, by = c("plate","assay")) %>%
  left_join(assay_fallback, by = "assay") %>%
  mutate(p0_per_plate = dplyr::coalesce(p0_upper, p0_upper_fallback)) %>%
  # final guard if coalesce still NA (shouldn’t happen if we had any blanks)
  mutate(p0_per_plate = ifelse(is.na(p0_per_plate), 0, p0_per_plate)) %>%
  group_by(Sample, assay) %>%
  slice_max(p0_per_plate, n = 1, with_ties = FALSE) %>%
  summarise(p0_use = first(p0_per_plate),
            plate_p0_max = first(plate),
            .groups = "drop")

# 5) Attach the conservative blank rate to each sample’s pooled counts, then compute LoB
counts_lob <- merged %>%
  rename(assay = ExperimentType, n_tot = Drop_pool, x_mut = Pos_pool) %>%
  left_join(sample_lob_p0, by = c("Sample","assay")) %>%
  mutate(p0_use = ifelse(is.na(p0_use), 0, p0_use)) %>%
  rowwise() %>%

  # This asks: given this sample’s total droplets and the blank/background probability,
  # what mutant-positive count is still within the 95th percentile of expected blank noise?
  mutate(
    LoB_count = qbinom(0.95, size = n_tot, prob = p0_use),
    # fractional-abundance scale
    LoB_FA    = ifelse(n_tot > 0, fa_scale * LoB_count / n_tot, NA_real_),
    # flag detections above LoB
    detected_LoB  = x_mut > LoB_count
  ) %>%
  ungroup() %>%
  # Keep only LoB outputs
  select(Sample, assay, LoB_count, LoB_FA, detected_LoB)


# --------------------------------------------------------
# add LoB outputs
# --------------------------------------------------------

merged <- merged %>%
  left_join(counts_lob, by = c("Sample" = "Sample", "ExperimentType" = "assay"))

# --------------------------------------------------------
# Passed Limit of Detection?
# --------------------------------------------------------

# empirically determined LoD
lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)

final <- merged %>%
  mutate(
    detected_LoD = CI_lower > lod_cut[ExperimentType]
  )

# --------------------------------------------------------
# Keep what you need for the dot-plot (+ LoB fields)
# --------------------------------------------------------

# main cleaned sample-region × mutation ddPCR results table
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
# partition-count denominator table for downstream genome-equivalent summaries
# --------------------------------------------------------

# Build an auxiliary denominator table from the same cleaned, pooled sample set
# used for SNV_data_final.xlsx. This table keeps the underlying Ch1/Ch2
# partition counts needed by the haploid-genome supplementary script.

partition_rows <- data.mut.brief %>%
  mutate(
    # Double-positive and signal-negative droplets
    # can be read directly from the channel classes.
    n_double_positive_droplets = as.numeric(`Ch1+Ch2+`),
    n_signal_negative_droplets = as.numeric(`Ch1-Ch2-`),

    # Convert the single-positive partitions to REF-only and MUT-only counts
    # using each assay's channel orientation:
    #   D178N/E200K: Ch1 = REF, Ch2 = MUT
    #   P102L:       Ch1 = MUT, Ch2 = REF
    n_ref_only_droplets = case_when(
      ExperimentType == "P102L" ~ as.numeric(`Ch1-Ch2+`),
      TRUE ~ as.numeric(`Ch1+Ch2-`)
    ),
    n_mut_only_droplets = case_when(
      ExperimentType == "P102L" ~ as.numeric(`Ch1+Ch2-`),
      TRUE ~ as.numeric(`Ch1-Ch2+`)
    ),

    # A double-positive droplet is positive for both target channels, so it
    # contributes to both the REF-positive and MUT-positive droplet totals.
    # Single-positive droplets contribute only to the target on their channel.
    n_ref_positive_droplets = n_double_positive_droplets + n_ref_only_droplets,
    n_mut_positive_droplets = n_double_positive_droplets + n_mut_only_droplets,
    n_signal_positive_droplets =
      n_double_positive_droplets + n_ref_only_droplets + n_mut_only_droplets
  )

# Before pooling wells, check that the reconstructed partition counts match the
# source counts for each row:
#   1. positive partitions + signal-negative droplets = accepted droplets
#   2. double-positive + mutant-only droplets = QuantaSoft Positives
partition_row_errors <- partition_rows %>%
  filter(
    n_signal_positive_droplets + n_signal_negative_droplets != AcceptedDroplets |
      n_mut_positive_droplets != Positives
  ) %>%
  select(Sample, ExperimentType, AcceptedDroplets, Positives,
         n_signal_positive_droplets, n_signal_negative_droplets,
         n_mut_positive_droplets)

if (nrow(partition_row_errors) > 0) {
  stop("Partition-count validation failed before aggregation:\n",
       paste(utils::capture.output(print(partition_row_errors)), collapse = "\n"))
}

# Use one row for each Sample and mutation assay. If the same Sample/assay was
# run in multiple wells or on multiple dates, sum those droplet counts here.
partition_counts <- partition_rows %>%
  group_by(Sample, ExperimentType) %>%
  summarise(
    n_wells = n(),
    n_accepted_droplets = sum(AcceptedDroplets, na.rm = TRUE),
    n_double_positive_droplets = sum(n_double_positive_droplets, na.rm = TRUE),
    n_ref_only_droplets = sum(n_ref_only_droplets, na.rm = TRUE),
    n_mut_only_droplets = sum(n_mut_only_droplets, na.rm = TRUE),
    n_signal_positive_droplets = sum(n_signal_positive_droplets, na.rm = TRUE),
    n_signal_negative_droplets = sum(n_signal_negative_droplets, na.rm = TRUE),
    n_ref_positive_droplets = sum(n_ref_positive_droplets, na.rm = TRUE),
    n_mut_positive_droplets = sum(n_mut_positive_droplets, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(sample_id = Sample, mutation = ExperimentType) %>%
  separate(sample_id, into = c("code", "brain_region"), sep = "_",
           fill = "right", remove = FALSE) %>%
  left_join(patient.names, by = "code") %>%
  select(participant, group, histotype, code, brain_region, sample_id, mutation,
         n_wells, n_accepted_droplets,
         n_double_positive_droplets, n_ref_only_droplets, n_mut_only_droplets,
         n_signal_positive_droplets, n_signal_negative_droplets,
         n_ref_positive_droplets, n_mut_positive_droplets)

# Check that the denominator table has the same sample-region x mutation rows
# as ddpcr.plot.data, and that its total droplets and mutant-positive droplets
# match the values exported in SNV_data_final.xlsx.
partition_counts_check <- partition_counts %>%
  # Bring in the already-exported totals for the same sample-region x mutation.
  left_join(
    ddpcr.plot.data %>%
      select(code, brain_region, mutation, n_mut_droplets, n_total_droplets),
    by = c("code", "brain_region", "mutation")
  )

# Rows with no joined total exist in the denominator table but not in
# SNV_data_final.xlsx.
partition_missing <- partition_counts_check %>%
  filter(is.na(n_total_droplets))

# Stop if the denominator table contains any sample-region x mutation that the
# final SNV workbook does not contain.
if (nrow(partition_missing) > 0) {
  stop("Partition-count rows did not match SNV_data_final rows:\n",
       paste(utils::capture.output(print(partition_missing)), collapse = "\n"))
}

# Among matched rows, check that the denominator table reproduces the two
# exported count fields used by downstream summaries.
partition_mismatch <- partition_counts_check %>%
  filter(
    n_accepted_droplets != n_total_droplets |
      n_mut_positive_droplets != n_mut_droplets
  ) %>%
  select(participant, code, brain_region, mutation,
         n_accepted_droplets, n_total_droplets,
         n_mut_positive_droplets, n_mut_droplets)

# Stop if any matched row has different count totals, or if the two tables have
# different row counts after the missing-row check above.
if (nrow(partition_mismatch) > 0 || nrow(partition_counts) != nrow(ddpcr.plot.data)) {
  stop("Partition-count totals do not match SNV_data_final:\n",
       paste(utils::capture.output(print(partition_mismatch)), collapse = "\n"))
}

# --------------------------------------------------------
# export
# --------------------------------------------------------

# Raw partition-count denominator input for the downstream supplementary
# summary. The existing workbook exports below deliberately keep their original
# columns.
write.csv(partition_counts,
          file.path(output_dir, "ddpcr_partition_counts_by_sample_assay.csv"),
          row.names = FALSE)

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
    pos_total     = sum(Pos_pool,     na.rm = TRUE),
    ref_pos_total = sum(Ref_pos_pool, na.rm = TRUE),
    drop_total    = sum(Drop_pool,    na.rm = TRUE),
    .groups = "drop"
  )

# drop 14-2 E200K as this sample was found to be heterozygous mutant E200K
pooled_participant_assay <- pooled_participant_assay %>%
  filter(!(code == "14-2" & assay == "E200K"))

# compute fractional abundance and confidence intervals using the same
# molecule concentration-ratio method as the sample-region table.
participant_fa <- purrr::pmap_dfr(
  list(
    ref_positive = pooled_participant_assay$ref_pos_total,
    ref_negative = pooled_participant_assay$drop_total - pooled_participant_assay$ref_pos_total,
    mut_positive = pooled_participant_assay$pos_total,
    mut_negative = pooled_participant_assay$drop_total - pooled_participant_assay$pos_total,
    total = pooled_participant_assay$drop_total
  ),
  fractional_abundance
)

pooled_participant_assay$FA <- participant_fa$fractional_abundance
pooled_participant_assay$ci_lower <- participant_fa$ci_low
pooled_participant_assay$ci_upper <- participant_fa$ci_high

# --------------------------------------------------------
# create mapping dataframe of plates x assays with respective p0 
# --------------------------------------------------------

# plate mapping (we need to know every plate contributing droplets before
# sample-region rows were pooled)
plate_mapping <- data.mut.brief %>%
  mutate(code = sub("_.*$", "", Sample),
         assay = ExperimentType,
         plate = Date) %>%
  distinct(code, assay, plate)

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
    LoB_FA             = fa_scale * LoB_count / drop_total,
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
