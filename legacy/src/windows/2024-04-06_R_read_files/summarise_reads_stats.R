setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2024-04-06_R_read_files")

library(readr)
#library(readxl)
library(tidyverse)
#library(janitor)
library(RColorBrewer)
#library(lmerTest)
#library(ggthemes)
#library(cowplot)
library(grid)
library(ggpubr)
library(viridis)

# Read the file as text
file_content <- readLines("E200K_summary.tsv")

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

# data w/o CJD30, who is heterozygous for E200K
df <- data[!(data$PatID == "CJD30"),]

# plot AAF
aaf.boxplot <- 
  
  ggplot(df, mapping = aes(x = Group, y = AAF, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.2) + #add horizontal bars to boxplot whiskers
  geom_boxplot(outlier.shape = NA) + #as outliers are seen in jitter
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  scale_color_viridis(discrete = TRUE, alpha=0.6) + 
  ylab(paste0("E200K alternate allele fraction")) +
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
  ylab(paste0("E200K alternate allele fraction")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank()) +
  stat_compare_means(method = "wilcox.test")

plot(aaf.filtered.boxplot)
