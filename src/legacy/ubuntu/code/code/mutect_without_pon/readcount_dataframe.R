#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)

# Function to process files from a directory and create a combined dataset
process_directory <- function(dir_path, pattern = "*.tsv") {
  # Set working directory
  setwd(dir_path)
  cat("Processing directory:", dir_path, "\n")
  
  # Define column types explicitly to ensure consistency
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
  
  # Get list of all TSV files
  tsvs <- list.files(pattern = pattern)
  cat("Found", length(tsvs), "files\n")
  
  # Create empty data frame with proper column types
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
  
  # Add a progress counter
  total_files <- length(tsvs)
  
  # Loop through each file
  for (i in seq_along(tsvs)) {
    file <- tsvs[i]
    
    # Show progress
    cat(sprintf("  Processing file %d of %d: %s\n", i, total_files, file))
    
    # Try-catch block to handle potential errors in individual files
    tryCatch({
      # Read the TSV file with explicit column types
      tempdf <- read_tsv(file, col_names = TRUE, col_types = col_types, show_col_types = FALSE)
      
      # Extract the sample name from the file name
      sample_name <- gsub("_metrics.tsv", "", file)
      
      # Add a new column with the sample name as the first column
      tempdf <- tempdf %>% mutate(Sample = sample_name) %>%
        select(Sample, everything())
      
      # Append the data frame to the result data frame
      result_df <- bind_rows(result_df, tempdf)
      
    }, error = function(e) {
      cat("  Error processing file:", file, "\n")
      cat("  Error message:", conditionMessage(e), "\n")
    })
  }
  
  # Return the combined dataset
  return(result_df)
}

# Define paths
samples_dir <- "/add/results/no_PoN/metrics/"
dilutions_dir <- "/add/results/no_PoN/metrics_dil/"
output_dir <- "/add/results/no_PoN/readcount_summaries/"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Process samples directory
cat("\n=== Processing Samples ===\n")
samples_df <- process_directory(samples_dir)

# Process dilutions directory
cat("\n=== Processing Dilutions ===\n")
dilutions_df <- process_directory(dilutions_dir)

# Set working directory to output directory
setwd(output_dir)

# Save the combined samples dataset
samples_output_file <- "readcounts_samples.csv"
write_csv(samples_df, samples_output_file)
cat("\nSamples dataset saved to:", file.path(output_dir, samples_output_file), "\n")
cat("Total rows in samples dataset:", nrow(samples_df), "\n")
cat("Number of unique samples:", length(unique(samples_df$Sample)), "\n")

# Save the combined dilutions dataset
dilutions_output_file <- "readcounts_dilutions.csv"
write_csv(dilutions_df, dilutions_output_file)
cat("\nDilutions dataset saved to:", file.path(output_dir, dilutions_output_file), "\n")
cat("Total rows in dilutions dataset:", nrow(dilutions_df), "\n")
cat("Number of unique dilutions:", length(unique(dilutions_df$Sample)), "\n")

cat("\nProcessing complete!\n")