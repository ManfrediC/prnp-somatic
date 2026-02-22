library(tidyverse)

# Function to process files from a directory and create a combined dataset
process_directory <- function(dir_path, pattern = "*.tsv") {
  setwd(dir_path)
  cat("Processing directory:", dir_path, "\n")
  
  col_types <- cols(
    CHROM = col_character(),
    POS = col_double(),
    REF = col_character(),
    BASE = col_character(),
    COUNT = col_double(),
    MEAN_BQ = col_double(),
    MEAN_MQ = col_double(),
    MEAN_POS = col_double(),
    FWD = col_double(),
    REV = col_double(),
    FRAC = col_double(),
    MMLQ = col_double(),
    MEAN_READ_POS = col_double(),
    FWD_BQ = col_double(),
    REV_BQ = col_double(),
    FRD = col_double(),
    SB_RATIO = col_double()
  )
  
  tsvs <- list.files(pattern = pattern)
  cat("Found", length(tsvs), "files\n")
  
  result_df <- tibble(
    Sample = character(),
    CHROM = character(),
    POS = double(),
    REF = character(),
    BASE = character(),
    COUNT = double(),
    MEAN_BQ = double(),
    MEAN_MQ = double(),
    MEAN_POS = double(),
    FWD = double(),
    REV = double(),
    FRAC = double(),
    MMLQ = double(),
    MEAN_READ_POS = double(),
    FWD_BQ = double(),
    REV_BQ = double(),
    FRD = double(),
    SB_RATIO = double()
  )
  
  total_files <- length(tsvs)
  
  for (i in seq_along(tsvs)) {
    file <- tsvs[i]
    cat(sprintf("  Processing file %d of %d: %s\n", i, total_files, file))
    
    tryCatch({
      tempdf <- read_tsv(file, col_names = TRUE, col_types = col_types, show_col_types = FALSE)
      sample_name <- gsub("_metrics.tsv", "", file)
      tempdf <- tempdf %>% mutate(Sample = sample_name) %>% select(Sample, everything())
      result_df <- bind_rows(result_df, tempdf)
    }, error = function(e) {
      cat("  Error processing file:", file, "\n")
      cat("  Error message:", conditionMessage(e), "\n")
    })
  }
  
  return(result_df)
}

# Output directory
output_dir <- "C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-05-10_CJD_QC_analysis/output_CSV"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Process samples directory
samples_dir <- "C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-05-10_CJD_QC_analysis/metrics"
samples_df <- process_directory(samples_dir)

# Set working directory to output directory
setwd(output_dir)

# Save the combined samples dataset
samples_output_file <- "readcounts_samples.csv"
write_csv(samples_df, samples_output_file)
cat("\nSamples dataset saved to:", file.path(output_dir, samples_output_file), "\n")
