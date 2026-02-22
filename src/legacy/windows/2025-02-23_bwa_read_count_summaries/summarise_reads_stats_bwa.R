
library(readr)
library(tidyverse)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(viridis)
library(gtools)
library(broom)

input_dir = "C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-02-23_bwa_read_count_summaries/input"

output_dir = "C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-02-23_bwa_read_count_summaries/mutation_csv"

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

input_dir = "C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-02-23_bwa_read_count_summaries/mutation_csv"
output_dir = "C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-02-23_bwa_read_count_summaries/summary_df"

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

########################
# plots
########################


# plot AAF
aaf.boxplot <- 
  
  ggplot(df, mapping = aes(x = Group, y = AAF, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.2) + #add horizontal bars to boxplot whiskers
  geom_boxplot(outlier.shape = NA) + #as outliers are seen in jitter
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  scale_color_viridis(discrete = TRUE, alpha=0.6) + 
  ylab(paste0("D178N alternate allele fraction")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank()) +
  stat_compare_means(method = "wilcox.test")

plot(aaf.boxplot)

#plot while excluding reads n < 100

df_filtered <- df[df$FORMAT.DP >= 100, ]

aaf.filtered.boxplot <- 
  
  ggplot(df_filtered, mapping = aes(x = Group, y = AAF, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.2) + #add horizontal bars to boxplot whiskers
  geom_boxplot(outlier.shape = NA) + #as outliers are seen in jitter
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  scale_color_viridis(discrete = TRUE, alpha=0.6) + 
  ylab(paste0("D178N alternate allele fraction")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank()) +
  stat_compare_means(method = "wilcox.test")

plot(aaf.filtered.boxplot)


# with t-test

aaf.filtered.boxplot.ttest <- 
  
  ggplot(df_filtered, mapping = aes(x = Group, y = AAF, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.2) + #add horizontal bars to boxplot whiskers
  geom_boxplot(outlier.shape = NA) + #as outliers are seen in jitter
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  scale_color_viridis(discrete = TRUE, alpha=0.6) + 
  ylab(paste0("D178N alternate allele fraction")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank()) +
  stat_compare_means(method = "wilcox.test")

plot(aaf.filtered.boxplot.ttest)




########################
# define mutated status + stats on mutated frequency -> I don't think that this is legitimate as sample sizes too small
########################

df <- df %>%
  mutate(mutated = ifelse(AD.ALT == 0, FALSE, TRUE))

# Calculate the frequency of TRUE in df$mutated for both groups
frequency_stats <- df %>%
  group_by(Group) %>%
  summarise(
    total = n(),                              # Total number of samples in each group
    mutated_true_count = sum(mutated),         # Number of TRUE in mutated column
    mutated_true_freq = mean(mutated)          # Frequency of TRUE (as a proportion)
  )

# Display the frequency statistics
print(frequency_stats)






# Summarise mutation frequency for each group
mutation_stats <- df %>%
  group_by(Group) %>%
  summarise(
    mutated_true_count = sum(mutated),
    total = n(),
    mutated_true_freq = mean(mutated)
  )

# Fisher's Exact Test or Chi-Square Test
# Creating a contingency table for the two groups
mutation_table <- table(df$Group, df$mutated)

# Fisher's Exact Test (preferred for smaller sample sizes)
fisher_test <- fisher.test(mutation_table)

# OR Chi-Square Test
chi_sq_test <- chisq.test(mutation_table)

# Bar plot of mutation frequency with 95% confidence intervals
ggplot(mutation_stats, aes(x = Group, y = mutated_true_freq)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
  geom_errorbar(aes(ymin = prop.test(mutated_true_count, total)$conf.int[1],
                    ymax = prop.test(mutated_true_count, total)$conf.int[2]),
                width = 0.2) +
  labs(title = "Frequency of Mutations in CJD vs Ctrl",
       x = "Group",
       y = "Proportion of Mutated Samples") +
  theme_minimal()

# Display the test result
print(fisher_test)
# If Chi-Square Test was used
print(chi_sq_test)






