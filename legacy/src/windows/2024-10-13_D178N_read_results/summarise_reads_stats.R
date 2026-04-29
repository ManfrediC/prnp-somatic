setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2024-10-13_D178N_read_results")

library(readr)
library(tidyverse)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(viridis)
library(gtools)
library(broom)

########################
# prepare dataframe
########################

# Read the file as text
file_content <- readLines("D178N_summary.tsv")

# Replace consecutive spaces with a single tab character
fixed_content <- gsub("\\s+", "\t", file_content)

# Write the fixed content to a temporary file
temp_file <- tempfile()
writeLines(fixed_content, temp_file)

# Import the fixed TSV file
data <- read.delim(temp_file, header = TRUE)

# categorise the samples to CJD or Ctrl
data$Group <- ifelse(grepl("CJD", data$Sample), "CJD", 
                     ifelse(grepl("Ctrl", data$Sample), "Ctrl", NA))

# Extract the patterns "CJD" or "Ctrl" followed by numbers from the Sample column
pat_id <- gsub(".*-(CJD\\d+|Ctrl\\d+)_.*", "\\1", data$Sample)

# Assign the extracted patterns to the PatID column
data$PatID <- pat_id

#calculate AAF
data$AAF <- data$AD.ALT/data$FORMAT.DP

df <- data


########################
# export as clean csv file
########################

D178N <- df %>%
  select(sample = PatID, REF, ALT, read_count = FORMAT.DP, REF_count = AD.REF, ALT_count = AD.ALT, AAF) %>%
  mutate(sample = factor(sample, levels = mixedsort(unique(sample)))) %>%
  arrange(sample)

write.csv(D178N, file = "D178N_summary_table.csv", row.names = FALSE)

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






