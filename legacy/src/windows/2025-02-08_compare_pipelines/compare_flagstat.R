# plot_flagstat_comparison.R
# This script reads in a TSV file (e.g., summary_flagstat_consistent.tsv) with flagstat metrics,
# splits the "Sample" column into a base sample name and a Pipeline ("bwa" or "orig"),
# computes normalized metrics, and then creates side-by-side bar charts for comparison.

# Load required libraries
library(tidyverse)  # Includes dplyr, tidyr, and ggplot2

# Read the summary TSV file.
setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-02-08_compare_pipelines")

df <- read.delim("summary_flagstat_samples.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# The "Sample" column has names like "CJD11_bwa" or "CJD11_orig".
# We split it into two new columns: "SampleBase" and "Pipeline".
df <- df %>%
  mutate(
    Pipeline = ifelse(grepl("_bwa$", Sample), "bwa", "orig"),
    SampleBase = sub("_(bwa|orig)$", "", Sample)
  )

# Ensure the key columns are numeric.
# (The AWK script should have produced numeric values, but we ensure here.)
numeric_cols <- c("Total", "Mapped", "Paired", "ProperlyPaired", "Duplicates")
df[numeric_cols] <- lapply(df[numeric_cols], as.numeric)

# Compute normalized metrics:
# - PctMapped: percentage of total reads that are mapped.
# - PctProperlyPaired: percentage of paired reads that are properly paired.
# - PctDuplicates: percentage of total reads marked as duplicates.
df <- df %>%
  mutate(
    PctMapped = (Mapped / Total) * 100,
    PctProperlyPaired = (ProperlyPaired / Paired) * 100,
    PctDuplicates = (Duplicates / Total) * 100
  )

# For plotting, reshape the data to long format for these three metrics.
df_long <- df %>%
  select(SampleBase, Pipeline, PctMapped, PctProperlyPaired, PctDuplicates) %>%
  pivot_longer(cols = c(PctMapped, PctProperlyPaired, PctDuplicates), 
               names_to = "Metric", values_to = "Value")

# Create side-by-side bar charts comparing the pipelines for each metric.
plot <- ggplot(df_long, aes(x = SampleBase, y = Value, fill = Pipeline)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  labs(x = "Sample", y = "Percentage", title = "Normalized Flagstat Metrics by Pipeline") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

# Display the plot
print(plot)

# Optionally, save the plot to a file:
ggsave("flagstat_comparison.png", plot, width = 10, height = 12)
