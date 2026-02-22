# created on 2025-04-05
# here, I compare the read counts in CJD vs Ctrl for the 4 highly pathogenic mutaions

library(tidyverse)
library(ggpubr)

setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-02-23_bwa_read_count_summaries/summary_df")

# import table
count_table <- read.csv("all_mutation_counts_bwa.csv")

# assign categories "CJD" or "Ctrl"
count_table$category <- ifelse(grepl("CJD", count_table$sample), "CJD", "Ctrl")

# only keep columns containing "AAF"
count_table <- count_table %>%
  select(sample, category, contains("AAF"))  %>%
  select(sample, category, contains(c("E200K", "D178N", "P102L", "A117V")))

e200k <- count_table %>%
  select(sample, category, E200K_AAF)

d178n <- count_table %>%
  select(sample, category, D178N_AAF)

p102l <- count_table %>%
  select(sample, category, P102L_AAF)

a117v <- count_table %>%
  select(sample, category, A117V_AAF)

# exclude the E200K heterozygous sample CJD30
e200k <- e200k %>%
  filter(sample != "CJD30")

# E200K: boxplot and t-test
e200k_plot <- ggplot(e200k, aes(x = category, y = E200K_AAF, fill = category)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = "E200K alternate allele fraction (force-called reads)", x = "Category", y = "E200K AAF") +
  theme_classic() +
  stat_compare_means(method = "t.test", label.x = 1.4)
e200k_plot

# D178N: boxplot and t-test
d178n_plot <- ggplot(d178n, aes(x = category, y = D178N_AAF, fill = category)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = "D178N alternate allele fraction (force-called reads)", x = "Category", y = "D178N AAF") +
  theme_classic() +
  stat_compare_means(method = "t.test", label.x = 1.4)
d178n_plot

# P102L: boxplot and t-test
p102l_plot <- ggplot(p102l, aes(x = category, y = P102L_AAF, fill = category)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = "P102L alternate allele fraction (force-called reads)", x = "Category", y = "P102L AAF") +
  theme_classic() +
  stat_compare_means(method = "t.test", label.x = 1.4)
p102l_plot

# A117V: boxplot and t-test
a117v_plot <- ggplot(a117v, aes(x = category, y = A117V_AAF, fill = category)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = "A117V alternate allele fraction (force-called reads)", x = "Category", y = "A117V AAF") +
  theme_classic() +
  stat_compare_means(method = "t.test", label.x = 1.4)
a117v_plot

# sum of AAF
count_table$AAF_sum <- rowSums(count_table[, 3:6], na.rm = TRUE)

# remove the heterozygous sample CJD30
count_table2 <- count_table %>%
  filter(sample != "CJD30")

# AAF sum: boxplot and t-test
sum_plot <- ggplot(count_table2, aes(x = category, y = AAF_sum, fill = category)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = "Sum of alternate allele fractions (force-called reads)", x = "Category", y = "Sum of AAF") +
  theme_classic() +
  stat_compare_means(method = "t.test", label.x = 1.4)
sum_plot


# save plots as pdf
setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/final results/plots")
ggsave("E200K_AAF_forcecall_boxplot.pdf", plot = e200k_plot, width = 6, height = 4)
ggsave("D178N_AAF_forcecall_boxplot.pdf", plot = d178n_plot, width = 6, height = 4)
ggsave("P102L_AAF_forcecall_boxplot.pdf", plot = p102l_plot, width = 6, height = 4)
ggsave("A117V_AAF_forcecall_boxplot.pdf", plot = a117v_plot, width = 6, height = 4)



