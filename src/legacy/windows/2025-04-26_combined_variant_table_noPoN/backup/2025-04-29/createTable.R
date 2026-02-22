library(tidyverse)

setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-04-26_combined_variant_table_noPoN")

# list of TSV files in the working directory, excluding the summary file
file_list <- list.files(pattern = "*.tsv$")

# exclude from file_list
exclude <- c("summary_combined_variants.tsv", "noPoN_PRNP_PASS.tsv", "noPoN_TET2_PASS.tsv", "noPoN_TTN_PASS.tsv")
file_list <- setdiff(file_list, exclude)
rm(exclude)

summary_table <- NULL

# import TSVs
for(file in file_list){

  temp_table <- read_tsv(file, col_names = TRUE, col_types = cols(.default = col_character()))

  # rename columns to GT, DP and AD, removing the sample ID from the colname
  temp_table <- temp_table %>%
    rename_with(.fn = ~ str_remove(., "^.*\\."), .cols = ends_with(c(".GT", ".DP", ".AD")))

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

# separate AD into REF_count and ALT_count
summary_table <- summary_table %>%
  separate(AD, into = c("REF_count", "ALT_count"), sep = ",", convert = TRUE)

#convert DP to numeric
summary_table$DP <- as.numeric(summary_table$DP)

# calculate AAF (ALT_count / DP)
summary_table$AAF <- summary_table$ALT_count / summary_table$DP

#----------------------------------------------------
# extract FUNCOTATION fields
#----------------------------------------------------

# gene symbol
summary_table <- summary_table %>%
  mutate(gene = str_extract(FUNCOTATION, "(?<=\\[)[^|]+"))

# remove all rows that are not PRNP, TTN or TET2
summary_table <- summary_table %>%
  filter(gene %in% c("PRNP", "TTN", "TET2"))

# Location relative to gene
summary_table <- summary_table %>%
  mutate(
    location_relative = word(FUNCOTATION, 6, sep = fixed("|"))
  )

# dbSNP ID
summary_table <- summary_table %>%
  mutate(dbsnp_id = str_extract(FUNCOTATION, "rs\\d+"))

# mutation type
summary_table <- summary_table %>%
  mutate(mutation_type = str_extract(FUNCOTATION, "\\b(SILENT|INTRON|MISSENSE|IGR|INTERGENIC|UPSTREAM|DOWNSTREAM|FIVE_PRIME_UTR|
  THREE_PRIME_UTR|SPLICE_SITE|SPLICE_DONOR|SPLICE_ACCEPTOR|FRAME_SHIFT|STOP_GAINED|STOP_LOST|START_LOST|START_GAINED|SYNONYMOUS|
  NONSENSE|INFRAME_DELETION|INFRAME_INSERTION|COMPLEX_SUBSTITUTION|TRANSCRIPT_ABLATION|TRANSCRIPT_AMPLIFICATION|
                                     REGULATORY_REGION_VARIANT|TF_BINDING_SITE_VARIANT)\\b"))

# create column "variant" 
summary_table <- summary_table %>%
  add_column(variant = NA_character_, .after = 1)

# annotate pathogenic PRNP variants
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



# for PRNP: annotate whether in intron or protein coding region
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
  relocate(region, .after = colnames(summary_table)[4])

# add population frequencies manually for dbSNP IDs (I checked gnomAD v4.1.0 and dbSNP for the population frequencies of the variants)

summary_table <- summary_table %>%
  mutate(population_frequency_gnomAD_manual = case_when(
    dbsnp_id == "rs1012397146" ~ 0.00001318,
    dbsnp_id == "rs1931969" ~ 0.2146,
    dbsnp_id == "rs2093390" ~ 0.3055,
    dbsnp_id == "rs2093391" ~ 0.2583,
    dbsnp_id == "rs6037932" ~ 0.2872,
    dbsnp_id == "rs6052769" ~ 0.4115,
    dbsnp_id == "rs6052772" ~ 0.2361,
    dbsnp_id == "rs6107516" ~ 0.2271,
    dbsnp_id == "rs781755787" ~ 0.000001369,
    dbsnp_id == "rs996098774" ~ NA_real_,
    TRUE ~ NA_real_
  ))

# rs996098774 appears to be an incorrect ID, so replace it with NA in the table
summary_table <- summary_table %>%
  mutate(dbsnp_id = case_when(
    dbsnp_id == "rs996098774" ~ NA_character_,
    TRUE ~ dbsnp_id
  ))

# for the remaining variants: build a variant ID that can be used in gnomAD v4.1.0 search
summary_table <- summary_table %>% 
  mutate(variant_id = str_c(str_remove(CHROM, "^chr"), POS, REF, ALT, sep = "-"))

# now add the frequencies manually: none of them were found in gnomAD -> NA

summary_table <- summary_table %>%
  mutate(population_frequency_gnomAD_manual = case_when(
    variant_id %in% c(
      "2-178581749-C-T", "2-178585952-T-A", "2-178585984-A-T",
      "2-178649895-C-T", "2-178649897-T-A", "2-178698399-T-A",
      "2-178698412-G-A", "20-4691495-C-G", "20-4691920-G-A",
      "20-4693834-CA-AG", "20-4694249-T-C", "20-4694270-A-G",
      "20-4695237-T-C", "20-4695357-G-C", "4-105175251-A-G",
      "4-105175261-T-G", "4-105215846-T-G", "4-105261420-A-T",
      "4-105273277-G-A"
    ) ~ NA_real_,
    TRUE ~ population_frequency_gnomAD_manual
  ))


  

# save to file
write_tsv(summary_table, "summary_combined_variants.tsv")


### single genes

# PRNP
prnp <- summary_table %>%
  filter(gene == "PRNP")

write_tsv(prnp, "noPoN_PRNP_PASS.tsv")

# TET2
tet2 <- summary_table %>%
  filter(gene == "TET2")

write_tsv(prnp, "noPoN_TET2_PASS.tsv")

# TTN
ttn <- summary_table %>%
  filter(gene == "TTN")

write_tsv(prnp, "noPoN_TTN_PASS.tsv")

