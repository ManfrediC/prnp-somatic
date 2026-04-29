# 2025-03-22
# this script compares the force call results of the original and bwa pipeline
# goal: which pipeline has higher read counts and better detection limit

library(tidyverse)

# bwa
setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-03-15_dilution_force_call_summaries/bwa/summary_df")
bwa <- read.csv("all_mutation_counts_bwa.csv")

# original 
setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-03-15_dilution_force_call_summaries/original/summary_df")
original <- read.csv("all_mutation_counts_original.csv")

# Order the rows by the sample column
bwa <- bwa[order(bwa$sample), ]
original <- original[order(original$sample), ]

# NA100_1to2 is somehow missing in original
bwa <- bwa %>% filter(!grepl("NA100_1to2", sample))

### subtract: bwa - original

# Merge the dataframes on the 'sample' column
merged_df <- merge(bwa, original, by = "sample", suffixes = c("_bwa", "_original"))

# Subtract the corresponding cells
result_df <- merged_df
for (col in names(bwa)) {
  if (col != "sample") {
    result_df[[col]] <- merged_df[[paste0(col, "_bwa")]] - merged_df[[paste0(col, "_original")]]
  }
}

# Keep only the necessary columns
result_df <- result_df[, names(bwa)]


# I concluded that the bwa pipeline is better. Continue with analysis



# Subset the data frame for samples of interest
subset_bwa <- subset(bwa, sample %in% c("A100_1to2", "NA100_1to10"))

# Identify columns containing "AAF" in their names
aaf_cols <- grep("AAF", names(subset_bwa), value = TRUE)

### gDNA AAFs
gDNA_NA_only <- subset(bwa, sample == "NA100_1to10")
gDNA_AAF <- gDNA_NA_only[aaf_cols]
gDNA_AAF_long <- pivot_longer(gDNA_AAF, cols = everything(), names_to = "mutation", values_to = "AAF")
gDNA_AAF_long$mutation <- sub("_AAF$", "", gDNA_AAF_long$mutation)

### A117V (100%) AAFs
A117V_only <- subset(bwa, sample == "A100_1to2")
A117V_AAF <- A117V_only[aaf_cols]
A117V_AAF_long <- pivot_longer(A117V_AAF, cols = everything(), names_to = "mutation", values_to = "AAF")
A117V_AAF_long$mutation <- sub("_AAF$", "", A117V_AAF_long$mutation)
