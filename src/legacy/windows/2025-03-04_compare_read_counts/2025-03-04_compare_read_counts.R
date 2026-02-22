
library(readr)
library(tidyverse)
library(ace_tools)

#-------------------------------------
# compare read count numbers
#-------------------------------------

# import summary CSVs

# bwa (aggressive trimming)
setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-02-23_bwa_read_count_summaries/summary_df")

bwa <- read.csv("all_mutation_counts_bwa.csv")

# old (minimap and normal trimming)
setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2024-10-19_read_count_summaries/summary_df")

old <- read.csv("all_mutation_counts.csv")

setdiff(old$sample, bwa$sample)  # Samples in 'old' but not in 'bwa'
setdiff(bwa$sample, old$sample)  # Samples in 'bwa' but not in 'old'

# for some reason, CJD26 is missing in "old"

### calculate difference

# Ensure both dataframes have the same samples before proceeding
common_samples <- intersect(old$sample, bwa$sample)

# Filter both dataframes to include only the common samples
old_filtered <- old[old$sample %in% common_samples, ]
bwa_filtered <- bwa[bwa$sample %in% common_samples, ]

# Ensure same order before subtraction
old_filtered <- old_filtered[match(common_samples, old_filtered$sample), ]
bwa_filtered <- bwa_filtered[match(common_samples, bwa_filtered$sample), ]

# Subtract old from bwa (excluding the sample column)
diff_df <- bwa_filtered
diff_df[-1] <- bwa_filtered[-1] - old_filtered[-1]

# Display results
print(diff_df)

# export
setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-03-04_compare_read_counts")
write.csv(diff_df, "read_count_diff.csv", row.names = FALSE)

#-------------------------------------
# count number of mutations over threshold
#-------------------------------------

mydf <- old
threshold = 0.001 # this needs to adapted based on the limit of detection derived from the dilution experiment

# columns ending in "AAF" -> subset
aaf_cols <- grep("AAF$", colnames(bwa), value = TRUE)
df_aaf <- mydf[, c("sample", aaf_cols)]

# Identify cells with values >= threshold
high_aaf <- df_aaf[, -1] >= threshold  # Exclude "sample" column

# Get row and column indices of matches
matches <- which(high_aaf, arr.ind = TRUE)

# Create a dataframe showing the sample names, column names, and values
result <- data.frame(
  sample = df_aaf$sample[matches[, 1]],
  column = colnames(df_aaf)[matches[, 2] + 1],  # +1 to adjust for sample column
  value = as.vector(df_aaf[, -1][matches])
)

print(result)








