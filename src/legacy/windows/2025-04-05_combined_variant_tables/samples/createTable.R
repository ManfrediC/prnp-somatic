library(tidyverse)

setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-04-05_combined_variant_tables/samples")

# list of TSV files in the working directory, excluding the summary file
file_list <- list.files(pattern = "*.tsv$")
file_list <- file_list[file_list != "summary_combined_variants.tsv"]

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
    CHROM == "chr20" & POS == "4699818" ~ "E200K_position",
    CHROM == "chr20" & POS == "4699525" ~ "P102L_position",
    CHROM == "chr20" & POS == "4699570" ~ "A117V_position",
    CHROM == "chr20" & POS == "4699752" ~ "D178N_position",
    CHROM == "chr20" & POS == "4699915" ~ "M232R_position",
    CHROM == "chr20" & POS == "4699758" ~ "V180I_position",
    CHROM == "chr20" & POS == "4699848" ~ "V210I_position",
    CHROM == "chr20" & POS == "4699842" ~ "R208C_position",
    CHROM == "chr20" & POS == "4699843" ~ "R208H_position",
    CHROM == "chr20" & POS == "4699534" ~ "P105L_position",
    CHROM == "chr20" & POS == "4699612" ~ "G131V_position",
    CHROM == "chr20" & POS == "4699618" ~ "A133V_position",
    CHROM == "chr20" & POS == "4699767" ~ "T183A_position",
    CHROM == "chr20" & POS == "4699812" ~ "F198V_position",
    CHROM == "chr20" & POS == "4699813" ~ "F198S_position",
    CHROM == "chr20" & POS == "4699870" ~ "G217R_position",
    CHROM == "chr20" & POS == "4699605" ~ "M129V_position",
    TRUE ~ NA_character_
  )) %>%
  relocate(variant, .after = sample)

# annotate whether in intron or protein coding region
summary_table <- summary_table %>%
  mutate(region = case_when(
    CHROM == "chr20" & POS <= 4686455 ~ "5' upstream",
    CHROM == "chr20" & POS >= 4686456 & POS <= 4686512 ~ "exon 1",
    CHROM == "chr20" & POS >= 4686513 & POS <= 4699210 ~ "intron",
    CHROM == "chr20" & POS >= 4699221 & POS <= 4699982 ~ "protein coding",
    CHROM == "chr20" & POS >= 4699983 & POS <= 4701588 ~ "exon 2 downstream of ORF",
    CHROM == "chr20" & POS >= 4701589 ~ "3' downstream"
  ))

summary_table <- summary_table %>%
  relocate(region, .after = colnames(summary_table)[2])
  



# variant	dbSNP id	position_GRCh38	REF	ALT
# E200K	rs28933385	4699818	G	A
# P102L	rs74315401 	4699525	C	T
# A117V	rs74315402	4699570	C	T
# D178N	rs74315403	4699752	G	A
# M232R	rs74315409	4699915	T	G
# V180I	rs74315408	4699758	G	A
# V210I	rs74315407	4699848	G	A
# R208C	rs55826236	4699842	C	T
# R208H	rs74315412	4699843	G	A
# P105L_A	rs11538758	4699534	C	A
# P105L_T	rs11538758	4699534	C	T
# G131V_A	rs74315410	4699612	G	A
# G131V_T	rs74315410	4699612	G	T
# A133V	rs74315415	4699618	C	T
# T183A	rs74315411	4699767	A	G
# F198V	rs55871421	4699812	T	G
# F198S	rs74315405	4699813	T	C
# G217R	rs74315406	4699870	A	G


# save to file
write_tsv(summary_table, "summary_combined_variants.tsv")
