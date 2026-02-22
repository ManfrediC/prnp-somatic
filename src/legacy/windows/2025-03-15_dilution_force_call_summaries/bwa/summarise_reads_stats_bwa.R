library(readr)
library(tidyverse)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(viridis)
library(gtools)
library(broom)

input_dir = "C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-03-15_dilution_force_call_summaries/bwa/input"
output_dir = "C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-03-15_dilution_force_call_summaries/bwa/mutation_csv"

setwd(input_dir)

########################
# create summary of read counts for every mutation
########################

tsv_files <- list.files(pattern = "\\.tsv$")

### LOOP STARTS
for (tsv in tsv_files) {
  
  # Read the file as text
  setwd(input_dir)
  file_content <- readLines(tsv)
  
  # Replace consecutive spaces with a single tab character
  fixed_content <- gsub("\\s+", "\t", file_content)
  
  # Write the fixed content to a temporary file
  temp_file <- tempfile()
  writeLines(fixed_content, temp_file)
  
  # Import the fixed TSV file
  data <- read.delim(temp_file, header = TRUE)
  
  # Replace dots with 0 in AD.ALT and ensure columns are numeric
  data$AD.ALT[data$AD.ALT == "."] <- 0
  data$AD.ALT <- as.numeric(data$AD.ALT)
  data$FORMAT.DP <- as.numeric(data$FORMAT.DP)
  
  # Categorize the samples to CJD or Ctrl
  data$Group <- ifelse(grepl("CJD", data$Sample), "CJD",
                       ifelse(grepl("Ctrl", data$Sample), "Ctrl", NA))
  
  # Extract the patterns "CJD" or "Ctrl" followed by numbers from the Sample column
  pat_id <- gsub(".*-(CJD\\d+|Ctrl\\d+)_.*", "\\1", data$Sample)
  
  # Assign the extracted patterns to the PatID column
  data$PatID <- pat_id
  
  # Calculate AAF
  data$AAF <- data$AD.ALT / data$FORMAT.DP
  
  # Create the summary dataframe
  summary_data <- data %>%
    select(sample = PatID, REF, ALT, read_count = FORMAT.DP, REF_count = AD.REF, ALT_count = AD.ALT, AAF) %>%
    mutate(sample = factor(sample, levels = mixedsort(unique(sample)))) %>%
    arrange(sample)
  
  # Create the output filename by replacing "_summary.tsv" with "_summary_table.csv"
  output_filename <- gsub("_summary.tsv", "_summary_table.csv", tsv)
  
  # output directory
  setwd(output_dir)
  
  # Write the summary data to a CSV file
  write.csv(summary_data, file = output_filename, row.names = FALSE)
  
  ### LOOP ENDS
}

rm(list = ls())

########################
# combine mutation read counts to single dataframe
########################

input_dir = "C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-03-15_dilution_force_call_summaries/bwa/mutation_csv"
output_dir = "C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-03-15_dilution_force_call_summaries/bwa/summary_df"

setwd(input_dir)

csv_summaries <- list.files(pattern = "\\.csv$")

# Initialize all_counts to NULL outside the loop
all_counts <- NULL

# Loop over the csv_summaries list
for (csv_file in csv_summaries) {
  
  # Read the CSV file
  mycsv <- read.csv(csv_file)
  
  # Extract mutation name by removing "_summary_table.csv" from the file name
  mutation <- gsub("_summary_table\\.csv$", "", csv_file)
  
  # Select the relevant columns
  mycsv <- subset(mycsv, select = c(sample, read_count, REF_count, ALT_count, AAF))
  
  # Add the mutation name to the column names
  colnames(mycsv) <- c("sample",
                       paste(mutation, "read_count", sep = "_"),
                       paste(mutation, "REF_count", sep = "_"),
                       paste(mutation, "ALT_count", sep = "_"),
                       paste(mutation, "AAF", sep = "_"))
  
  # Join with cumulative all_counts
  if (is.null(all_counts)) {
    all_counts <- mycsv  # Initialize with the first mycsv dataframe
  } else {
    all_counts <- left_join(all_counts, mycsv, by = "sample")  # Join on the "sample" column
  }
  
  # End of loop
}

setwd(output_dir)

write.csv(all_counts, "all_mutation_counts_bwa.csv", row.names = FALSE)

