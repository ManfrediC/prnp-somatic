#!/usr/bin/env Rscript
# ============================================================
# CJD+dilutions variant extraction + QC (with PoN)
# Adapted from legacy createTable_CJD.R with reproducible CLI I/O.
# ============================================================

fail <- function(msg) {
  stop(msg, call. = FALSE)
}

# Package check is explicit so failures happen before any partial output is written.
required_pkgs <- c("dplyr", "readr", "stringr", "tidyr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  fail(paste0(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "),
    ". Install them in your conda environment before running this step."
  ))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

parse_cli_args <- function(args) {
  # Lightweight --key value parser to keep runtime dependencies minimal.
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      fail(paste0("Unexpected argument format: ", key, " (expected --key value pairs)"))
    }
    if (i == length(args)) {
      fail(paste0("Missing value for argument: ", key))
    }
    value <- args[[i + 1L]]
    if (startsWith(value, "--")) {
      fail(paste0("Missing value for argument: ", key))
    }
    out[[sub("^--", "", key)]] <- value
    i <- i + 2L
  }
  out
}

get_arg <- function(parsed, name, default = NULL, required = FALSE) {
  if (!is.null(parsed[[name]])) {
    return(parsed[[name]])
  }
  if (required) {
    fail(paste0("Required argument missing: --", name))
  }
  default
}

parse_bool <- function(x, name) {
  x_low <- tolower(x)
  if (x_low %in% c("1", "true", "yes", "y")) return(TRUE)
  if (x_low %in% c("0", "false", "no", "n")) return(FALSE)
  fail(paste0("Invalid boolean value for ", name, ": ", x))
}

parse_num <- function(x, name) {
  out <- suppressWarnings(as.numeric(x))
  if (is.na(out)) {
    fail(paste0("Invalid numeric value for ", name, ": ", x))
  }
  out
}

first_numeric_from_list <- function(x) {
  # Handles multi-allelic comma lists by taking the first value, consistent with current workflow.
  vapply(
    x,
    function(value) {
      if (is.na(value) || value == "") {
        return(NA_real_)
      }
      token <- strsplit(value, ",", fixed = TRUE)[[1]][1]
      suppressWarnings(as.numeric(token))
    },
    numeric(1)
  )
}

# ----------------------------------------------------
# Parse CLI args
# ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
parsed <- parse_cli_args(args)

variant_dir <- get_arg(parsed, "variant-dir", required = TRUE)
metrics_dir <- get_arg(parsed, "metrics-dir", required = TRUE)
manual_freq_path <- get_arg(parsed, "manual-freq", required = TRUE)
output_dir <- get_arg(parsed, "output-dir", required = TRUE)

enable_aaf_filter <- parse_bool(get_arg(parsed, "enable-aaf-filter", "1"), "enable-aaf-filter")
aaf_threshold <- parse_num(get_arg(parsed, "aaf-threshold", "0.0081"), "aaf-threshold")

min_alt_count <- parse_num(get_arg(parsed, "min-alt-count", "10"), "min-alt-count")
min_dp <- parse_num(get_arg(parsed, "min-dp", "100"), "min-dp")
min_strand_alt <- parse_num(get_arg(parsed, "min-strand-alt", "3"), "min-strand-alt")
min_mean_bq <- parse_num(get_arg(parsed, "min-mean-bq", "20"), "min-mean-bq")
min_mean_mq <- parse_num(get_arg(parsed, "min-mean-mq", "20"), "min-mean-mq")
max_pop_freq <- parse_num(get_arg(parsed, "max-pop-freq", "0.001"), "max-pop-freq")
max_binom_p <- parse_num(get_arg(parsed, "max-binom-p", "1e-6"), "max-binom-p")

if (!dir.exists(variant_dir)) fail(paste0("Variant table directory does not exist: ", variant_dir))
if (!dir.exists(metrics_dir)) fail(paste0("Readcount metrics directory does not exist: ", metrics_dir))
if (!file.exists(manual_freq_path)) fail(paste0("Manual frequency TSV does not exist: ", manual_freq_path))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------
# 1) Import TSVs (VariantsToTable outputs)
# ----------------------------------------------------
# Keep the legacy exclusions so reruns do not accidentally re-import generated summaries.
file_list <- list.files(variant_dir, pattern = "\\.tsv$", full.names = TRUE)

exclude <- c(
  "summary_combined_variants.tsv", "withPoN_PRNP_PASS.tsv", "withPoN_TET2_PASS.tsv",
  "withPoN_TTN_PASS.tsv", "withPoN_PRNP_final.tsv", "final_withPoN_variants.tsv",
  "filtered_variants.tsv", "filtered_prnp_variants.tsv", "filter_counts.tsv", "run_settings.tsv"
)
file_list <- file_list[!(basename(file_list) %in% exclude)]

if (length(file_list) == 0) {
  fail(paste0("No variant TSV files found in: ", variant_dir))
}

summary_table <- NULL

for (file in file_list) {
  temp_table <- read_tsv(file, col_types = cols(.default = col_character()), show_col_types = FALSE)

  # rename columns to GT, DP and AD, removing the sample ID from the colname
  temp_table <- temp_table %>%
    rename_with(.fn = ~ str_remove(., "^.*\\."), .cols = ends_with(c(".GT", ".DP", ".AD", ".F1R2", ".F2R1", ".SB")))

  if (nrow(temp_table) == 0) next

  sample_name <- basename(file)
  sample_name <- str_remove(sample_name, "\\.withPoN\\.tsv$")
  sample_name <- str_remove(sample_name, "\\.tsv$")

  # Keep "sample" as the first column to match the legacy table layout.
  colnames(temp_table) <- gsub(paste0("^", sample_name, "\\."), "sample.", colnames(temp_table))

  temp_table$sample <- sample_name
  temp_table <- temp_table %>% select(sample, everything())

  if (is.null(summary_table)) {
    summary_table <- temp_table
  } else {
    summary_table <- bind_rows(summary_table, temp_table)
  }
}

if (is.null(summary_table) || nrow(summary_table) == 0) {
  fail("Variant import produced an empty summary table")
}

required_cols <- c("sample", "CHROM", "POS", "REF", "ALT", "FILTER", "FUNCOTATION", "DP", "AD", "SB")
missing_cols <- setdiff(required_cols, colnames(summary_table))
if (length(missing_cols) > 0) {
  fail(paste0("Variant table is missing required column(s): ", paste(missing_cols, collapse = ", ")))
}
if (!("GNOMAD_AF" %in% colnames(summary_table))) {
  summary_table$GNOMAD_AF <- NA_character_
}

# ----------------------------------------------------
# 2) Separate AD into REF_count and ALT_count
# ----------------------------------------------------
summary_table <- summary_table %>%
  separate(AD, into = c("REF_count", "ALT_count"), sep = ",", convert = TRUE, fill = "right", extra = "drop")

summary_table$DP <- as.numeric(summary_table$DP)
summary_table$POS <- as.numeric(summary_table$POS)
summary_table$REF_count <- as.numeric(summary_table$REF_count)
summary_table$ALT_count <- as.numeric(summary_table$ALT_count)
summary_table$AAF <- summary_table$ALT_count / summary_table$DP

# ----------------------------------------------------
# 3) Extract FUNCOTATION fields
# ----------------------------------------------------
# The parsing below intentionally follows the historical createTable workflow.

# Funcotator can include one entry per ALT; keep first entry for parsing.
funcotation_primary <- strsplit(coalesce(summary_table$FUNCOTATION, ""), ",", fixed = TRUE)
funcotation_primary <- vapply(funcotation_primary, function(x) if (length(x) == 0) "" else x[[1]], character(1))
funcotation_primary <- str_replace_all(funcotation_primary, "^\\[|\\]$", "")
summary_table$FUNCOTATION_PRIMARY <- funcotation_primary

# gene symbol
summary_table <- summary_table %>%
  mutate(gene = word(FUNCOTATION_PRIMARY, 1, sep = fixed("|")))

# remove all rows that are not PRNP, TTN or TET2
summary_table <- summary_table %>%
  filter(gene %in% c("PRNP", "TTN", "TET2"))

# location relative to gene
a <- summary_table %>%
  mutate(location_relative = word(FUNCOTATION_PRIMARY, 6, sep = fixed("|")))
summary_table <- a
rm(a)

# dbSNP ID
summary_table <- summary_table %>%
  mutate(dbsnp_id = str_extract(FUNCOTATION_PRIMARY, "rs\\d+"))

# mutation type
summary_table <- summary_table %>%
  mutate(mutation_type = str_extract(
    FUNCOTATION_PRIMARY,
    "\\b(SILENT|INTRON|MISSENSE|IGR|INTERGENIC|UPSTREAM|DOWNSTREAM|FIVE_PRIME_UTR|THREE_PRIME_UTR|SPLICE_SITE|SPLICE_DONOR|SPLICE_ACCEPTOR|FRAME_SHIFT|STOP_GAINED|STOP_LOST|START_LOST|START_GAINED|SYNONYMOUS|NONSENSE|INFRAME_DELETION|INFRAME_INSERTION|COMPLEX_SUBSTITUTION|TRANSCRIPT_ABLATION|TRANSCRIPT_AMPLIFICATION|REGULATORY_REGION_VARIANT|TF_BINDING_SITE_VARIANT)\\b"
  ))

# annotate pathogenic PRNP variants
summary_table <- summary_table %>%
  mutate(variant = case_when(
    CHROM == "chr20" & POS == 4699818 ~ "E200K_position",
    CHROM == "chr20" & POS == 4699525 ~ "P102L_position",
    CHROM == "chr20" & POS == 4699570 ~ "A117V_position",
    CHROM == "chr20" & POS == 4699752 ~ "D178N_position",
    CHROM == "chr20" & POS == 4699915 ~ "M232R_position",
    CHROM == "chr20" & POS == 4699758 ~ "V180I_position",
    CHROM == "chr20" & POS == 4699848 ~ "V210I_position",
    CHROM == "chr20" & POS == 4699842 ~ "R208C_position",
    CHROM == "chr20" & POS == 4699843 ~ "R208H_position",
    CHROM == "chr20" & POS == 4699534 ~ "P105L_position",
    CHROM == "chr20" & POS == 4699612 ~ "G131V_position",
    CHROM == "chr20" & POS == 4699618 ~ "A133V_position",
    CHROM == "chr20" & POS == 4699767 ~ "T183A_position",
    CHROM == "chr20" & POS == 4699812 ~ "F198V_position",
    CHROM == "chr20" & POS == 4699813 ~ "F198S_position",
    CHROM == "chr20" & POS == 4699870 ~ "G217R_position",
    CHROM == "chr20" & POS == 4699605 ~ "M129V_position",
    TRUE ~ NA_character_
  ))

# for PRNP: annotate whether in intron or protein coding region
summary_table <- summary_table %>%
  mutate(region = case_when(
    CHROM == "chr20" & POS <= 4686455 ~ "5' upstream",
    CHROM == "chr20" & POS >= 4686456 & POS <= 4686512 ~ "exon 1",
    CHROM == "chr20" & POS >= 4686513 & POS <= 4699210 ~ "intron",
    CHROM == "chr20" & POS >= 4699221 & POS <= 4699982 ~ "protein coding",
    CHROM == "chr20" & POS >= 4699983 & POS <= 4701588 ~ "exon 2 downstream of ORF",
    CHROM == "chr20" & POS >= 4701589 ~ "3' downstream",
    TRUE ~ NA_character_
  ))

# ----------------------------------------------------
# 4) Population frequencies (manual table + VCF GNOMAD_AF)
# ----------------------------------------------------

# rs996098774 appears to be an incorrect ID in the legacy workflow.
summary_table <- summary_table %>%
  mutate(dbsnp_id = if_else(dbsnp_id == "rs996098774", NA_character_, dbsnp_id))

# build variant ID used for manual lookup
summary_table <- summary_table %>%
  mutate(
    REF = toupper(REF),
    ALT = toupper(ALT),
    variant_id = str_c(str_remove(CHROM, "^chr"), POS, REF, ALT, sep = "-")
  )

manual_freq <- read_tsv(
  manual_freq_path,
  comment = "#",
  na = c("", "NA"),
  col_types = cols(
    variant_id = col_character(),
    dbsnp_id = col_character(),
    population_frequency = col_double(),
    source = col_character(),
    retrieved_date = col_character(),
    notes = col_character()
  ),
  show_col_types = FALSE
)

manual_by_dbsnp <- manual_freq %>%
  filter(!is.na(dbsnp_id), dbsnp_id != "") %>%
  select(dbsnp_id, population_frequency_gnomAD_manual_dbsnp = population_frequency) %>%
  distinct(dbsnp_id, .keep_all = TRUE)

manual_by_variant <- manual_freq %>%
  filter(!is.na(variant_id), variant_id != "") %>%
  select(variant_id, population_frequency_gnomAD_manual_variant = population_frequency) %>%
  distinct(variant_id, .keep_all = TRUE)

summary_table <- summary_table %>%
  left_join(manual_by_dbsnp, by = "dbsnp_id") %>%
  left_join(manual_by_variant, by = "variant_id") %>%
  mutate(
    population_frequency_gnomAD_manual = coalesce(
      population_frequency_gnomAD_manual_dbsnp,
      population_frequency_gnomAD_manual_variant
    )
  )

summary_table$GNOMAD_AF_table <- first_numeric_from_list(summary_table$GNOMAD_AF)

# Final population frequency used for filtering: prefer whichever frequency is higher if both exist.
summary_table <- summary_table %>%
  rowwise() %>%
  mutate(
    population_frequency = {
      # Use the highest available frequency as the conservative population estimate.
      vals <- c(GNOMAD_AF_table, population_frequency_gnomAD_manual)
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) NA_real_ else max(vals)
    }
  ) %>%
  ungroup()

# -----------------------------------------------------
# 5) Add base and read quality metrics (bam-readcount)
# -----------------------------------------------------
# Join key is sample_name + CHROM/POS/REF/ALT so base metrics stay variant-specific.
metric_files <- sort(list.files(metrics_dir, pattern = "_metrics\\.tsv$", full.names = TRUE))
if (length(metric_files) == 0) {
  fail(paste0("No *_metrics.tsv files found in: ", metrics_dir))
}

metrics_list <- vector("list", length(metric_files))
for (i in seq_along(metric_files)) {
  file <- metric_files[[i]]
  sample_name <- basename(file)
  sample_name <- sub("_metrics\\.tsv$", "", sample_name)

  temp <- read_tsv(
    file,
    col_types = cols(
      CHROM = col_character(),
      POS = col_double(),
      REF = col_character(),
      BASE = col_character(),
      COUNT = col_double(),
      MEAN_BQ = col_double(),
      MEAN_MQ = col_double(),
      .default = col_guess()
    ),
    show_col_types = FALSE
  )

  temp <- temp %>%
    rename(ALT = BASE) %>%
    mutate(sample_name = sample_name) %>%
    filter(COUNT != 0) %>%
    mutate(
      POS = as.numeric(POS),
      REF = toupper(REF),
      ALT = toupper(ALT)
    ) %>%
    select(sample_name, CHROM, POS, REF, ALT, MEAN_BQ, MEAN_MQ)

  metrics_list[[i]] <- temp
}
sample_basecount <- bind_rows(metrics_list)

summary_table$sample_name <- sub("\\.withPoN$", "", summary_table$sample)
summary_table$sample_name <- sub("\\.func\\.af$", "", summary_table$sample_name)
summary_table$REF <- toupper(summary_table$REF)
summary_table$ALT <- toupper(summary_table$ALT)

summary_table <- summary_table %>%
  left_join(sample_basecount, by = c("sample_name", "CHROM", "POS", "REF", "ALT"))

# -----------------------------------------------------
# 6) Apply QC criteria
# -----------------------------------------------------
# Criteria ordering is preserved so filter_counts reports remain easy to compare with legacy runs.

filtered_table <- summary_table %>%
  # 1) break SB into its four pieces
  separate(SB, into = c("SB_refF", "SB_refR", "SB_altF", "SB_altR"), sep = ",", convert = TRUE, fill = "right", extra = "drop") %>%
  # ALT >= min_alt_count and depth >= min_dp
  filter(ALT_count >= min_alt_count, DP >= min_dp) %>%
  # enforce strand balance: at least min_strand_alt ALT on each strand
  filter(SB_altF >= min_strand_alt, SB_altR >= min_strand_alt)

# minimum values for mapping and base quality
filtered_table <- filtered_table %>%
  filter(MEAN_BQ >= min_mean_bq & MEAN_MQ >= min_mean_mq)

# Binomial-test filter: remove calls with VAF ~= 0.5
filtered_table <- filtered_table %>%
  rowwise() %>%
  mutate(
    p_value = if (is.na(ALT_count) || is.na(DP) || DP <= 0) {
      NA_real_
    } else {
      binom.test(ALT_count, DP, p = 0.5)$p.value
    }
  ) %>%
  ungroup() %>%
  filter(!is.na(p_value), p_value <= max_binom_p)

# keep those with population frequency < max_pop_freq OR NA
filtered_table <- filtered_table %>%
  filter(is.na(population_frequency) | population_frequency < max_pop_freq)

# TEMPORARY (manual-check need): disable application of the AAF filter.
# if (enable_aaf_filter) {
#   filtered_final <- filtered_table %>%
#     filter(!is.na(AAF), AAF > aaf_threshold)
# } else {
#   filtered_final <- filtered_table
# }
filtered_final <- filtered_table

# inspect rows that were filtered out
filtered_out <- anti_join(
  summary_table,
  filtered_final,
  by = c("sample", "CHROM", "POS", "REF", "ALT")
)

# -----------------------------------------------------
# 6b) Sequential filter counts (for transparent reporting)
# -----------------------------------------------------
summary_split <- summary_table %>%
  separate(SB, into = c("SB_refF", "SB_refR", "SB_altF", "SB_altR"), sep = ",", convert = TRUE, fill = "right", extra = "drop")

n0 <- nrow(summary_split)

step1 <- summary_split %>% filter(ALT_count >= min_alt_count, DP >= min_dp)
n1 <- nrow(step1)

step2 <- step1 %>% filter(SB_altF >= min_strand_alt, SB_altR >= min_strand_alt)
n2 <- nrow(step2)

step3 <- step2 %>%
  filter(MEAN_BQ >= min_mean_bq, MEAN_MQ >= min_mean_mq) %>%
  rowwise() %>%
  mutate(p_value = if (is.na(ALT_count) || is.na(DP) || DP <= 0) NA_real_ else binom.test(ALT_count, DP, p = 0.5)$p.value) %>%
  ungroup() %>%
  filter(!is.na(p_value), p_value <= max_binom_p)
n3 <- nrow(step3)

step4 <- step3 %>%
  filter(is.na(population_frequency) | population_frequency < max_pop_freq)
n4 <- nrow(step4)

# TEMPORARY (manual-check need): keep AAF step unfiltered in sequential counts.
step5 <- step4
n5 <- n4
step5_label <- paste0("AAF filter temporarily disabled (configured threshold would be > ", aaf_threshold, ")")

filter_counts <- tibble(
  step = c(
    paste0("ALT >= ", min_alt_count, " and DP >= ", min_dp),
    paste0("strand balance >= ", min_strand_alt, " ALT reads per strand"),
    paste0("MEAN_BQ >= ", min_mean_bq, ", MEAN_MQ >= ", min_mean_mq, ", binomial p <= ", max_binom_p),
    paste0("population frequency < ", max_pop_freq, " or NA"),
    step5_label
  ),
  rows_before = c(n0, n1, n2, n3, n4),
  rows_after = c(n1, n2, n3, n4, n5),
  rows_removed = c(n0 - n1, n1 - n2, n2 - n3, n3 - n4, n4 - n5)
)
# filter_counts is used directly in methods/reporting write-ups.

# ------------------------------------------------
# 7) Save tables
# ------------------------------------------------
# Write both new reproducible outputs and legacy convenience files for continuity.
settings <- tibble(
  key = c(
    "enable_aaf_filter", "aaf_filter_applied", "aaf_threshold", "min_alt_count", "min_dp",
    "min_strand_alt", "min_mean_bq", "min_mean_mq", "max_pop_freq", "max_binom_p"
  ),
  value = c(
    as.character(enable_aaf_filter), "FALSE (temporary override)", as.character(aaf_threshold), as.character(min_alt_count), as.character(min_dp),
    as.character(min_strand_alt), as.character(min_mean_bq), as.character(min_mean_mq), as.character(max_pop_freq), as.character(max_binom_p)
  )
)

# primary outputs
write_tsv(summary_table, file.path(output_dir, "summary_combined_variants.tsv"), na = "")
write_tsv(filtered_final, file.path(output_dir, "filtered_variants.tsv"), na = "")
write_tsv(filtered_final %>% filter(gene == "PRNP"), file.path(output_dir, "filtered_prnp_variants.tsv"), na = "")
write_tsv(filtered_out, file.path(output_dir, "filtered_out_variants.tsv"), na = "")
write_tsv(filter_counts, file.path(output_dir, "filter_counts.tsv"), na = "")
write_tsv(settings, file.path(output_dir, "run_settings.tsv"), na = "")

# legacy-style convenience outputs
write_tsv(filtered_final, file.path(output_dir, "final_withPoN_variants.tsv"), na = "")
write_tsv(summary_table %>% filter(gene == "PRNP"), file.path(output_dir, "withPoN_PRNP_PASS.tsv"), na = "")
write_tsv(summary_table %>% filter(gene == "TET2"), file.path(output_dir, "withPoN_TET2_PASS.tsv"), na = "")
write_tsv(summary_table %>% filter(gene == "TTN"), file.path(output_dir, "withPoN_TTN_PASS.tsv"), na = "")

cat("Variant QC complete.\n")
cat("Summary table:  ", file.path(output_dir, "summary_combined_variants.tsv"), "\n", sep = "")
cat("Filtered table: ", file.path(output_dir, "filtered_variants.tsv"), "\n", sep = "")
cat("Filter counts:  ", file.path(output_dir, "filter_counts.tsv"), "\n", sep = "")
