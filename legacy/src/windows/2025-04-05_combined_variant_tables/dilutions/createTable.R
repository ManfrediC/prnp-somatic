library(tidyverse)

setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-04-05_combined_variant_tables")

file_list <- list.files(pattern = "*.tsv$")

summary_table <- NULL

for(file in file_list){
  temp_table <- read_tsv(file, col_names = TRUE, col_types = cols(.default = col_character()))
  
  if(nrow(temp_table) == 0) next  # skip empty files
  
  sample_name <- str_remove(file, "\\.tsv$")
  
  colnames(temp_table) <- gsub(paste0(sample_name, "\\."), "sample.", colnames(temp_table))
  
  temp_table$sample <- sample_name
  
  temp_table <- temp_table %>%
    select(sample, everything())
  
  if(is.null(summary_table)){
    summary_table <- temp_table
  } else {
    summary_table <- bind_rows(summary_table, temp_table)
  }
}

# create column "variant" 
summary_table <- summary_table %>%
  add_column(variant = NA_character_, .after = 1)

test <- summary_table

# annotate pathogenic variants
summary_table <- summary_table %>%
  mutate(variant = case_when(
    CHROM == "chr20" & POS == "4699818" ~ "E200K position",
    CHROM == "chr20" & POS == "4699570" ~ "A117V position",
    CHROM == "chr20" & POS == "4699570" ~ "A117V",
    CHROM == "chr20" & POS == "4699752" ~ "D178N position",
    CHROM == "chr20" & POS == "4699525" ~ "P102L position",
    TRUE ~ NA_character_
  )) %>%
  relocate(variant, .after = sample)

# maybe also add the other pathogenic variants (you have a VCF with all listed on Ubuntu)

# position REF ALT
# 4699818	G	A -> E200K
# 4699570	C	T -> A117V
# 4699752	G	A -> D178N
# 4699525	C	T ->  P102L



# save to file
write_tsv(summary_table, "summary_combined_variants.tsv")
