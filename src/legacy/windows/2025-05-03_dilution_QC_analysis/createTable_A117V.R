

library(tidyverse)
library(gt)
library(pagedown)

setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-05-03_dilution_QC_analysis")

# list of TSV files in the working directory, excluding the summary file
file_list <- list.files(pattern = "*.tsv$")

# exclude from file_list
exclude <- c("summary_combined_variants.tsv", "noPoN_PRNP_PASS.tsv", "noPoN_TET2_PASS.tsv", 
             "noPoN_TTN_PASS.tsv", "noPoN_PRNP_final.tsv", "final_noPoN_variants.tsv")
file_list <- setdiff(file_list, exclude)
rm(exclude)

summary_table <- NULL

# import TSVs
for(file in file_list){

  temp_table <- read_tsv(file, col_names = TRUE, col_types = cols(.default = col_character()))

  # rename columns to GT, DP and AD, removing the sample ID from the colname
  temp_table <- temp_table %>%
    rename_with(.fn = ~ str_remove(., "^.*\\."), .cols = ends_with(c(".GT", ".DP", ".AD", ".F1R2", ".F2R1", ".SB")))

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
rm(file, file_list, sample_name, temp_table)

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

#-----------------------------------------------------
# add base and read quality metrics (derived from base-readcount)
#-----------------------------------------------------

# import base-readcount CSV summary
setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-05-02_readcount_TSV_noPoN_samples/output")

sample_basecount <- read.csv("readcounts_dilutions.csv", header = TRUE)

# rename columns in sample_basecount 
sample_basecount <- sample_basecount %>%
  rename(ALT = BASE) %>% # BASE to ALT
  rename(sample_name = Sample) # Sample to sample_name

# sample_basecount: remove empty rows (COUNT == 0)
sample_basecount <- sample_basecount %>%
  filter(COUNT != 0)

# sample_basecount: remove all variants that aren't A117V
sample_basecount <- sample_basecount %>%
  filter(CHROM == "chr20" & POS == 4699570 & REF == "C" & ALT == "T")

# POS to numeric in summary_table and sample_basecount
summary_table$POS <- as.numeric(summary_table$POS)
sample_basecount$POS <- as.numeric(sample_basecount$POS)

# in summary_table, sample_name without ".noPoN" in position 2
summary_table$sample_name <- sub(".noPoN", "", summary_table$sample)
summary_table <- summary_table %>% 
  relocate(sample_name, .after = sample)

# in sample_basecount, change REF case to uppercase for consistency
sample_basecount$REF <- toupper(sample_basecount$REF)

# left join sample_basecount quality values (MEAN_BQ and MEAN_MQ) to summary_table (matching according to sample and variant coordinates)
sample_basecount_selected <- sample_basecount %>%
  select(sample_name, CHROM, POS, REF, ALT, MEAN_BQ, MEAN_MQ)

# the spike in had the base substitution CA -> TG. As the quality values are for POS 4699570 (the pathogenic C->T) alone, 
# I'm omitting the REF and ALT columns from the join (A -> G is synonymous for Valine)
summary_table <- summary_table %>%
  left_join(sample_basecount_selected, by = c("sample_name", "CHROM", "POS"))

# keep only A117V
A117V <- summary_table %>%
  filter(variant == "A117V_position")

# fix REF and ALT
A117V <- A117V %>%
  rename(REF = REF.y, ALT = ALT.y) %>%
  select(-REF.x, -ALT.x)

# add dbSNP ID
A117V$dbsnp_id <- "rs74315402"

#-----------------------------------------------------
# apply QC criteria
#-----------------------------------------------------

filtered_table <- A117V %>%
  # 1) break SB into its four pieces
  separate(SB,
           into = c("SB_refF","SB_refR","SB_altF","SB_altR"),
           sep   = ",",
           convert = TRUE) %>%

  # ALT ≥ 10  depth ≥100
  filter(
    ALT_count  >= 10,
    DP         >= 100) %>%
  
  # enforce strand balance: at least 3 ALT on each strand
  filter(  
    SB_altF    >= 3,
    SB_altR    >= 3
  )

# minimum values of 20 for mapping and base quality
filtered_table <- filtered_table %>%
  filter(MEAN_BQ >= 20 & MEAN_MQ >= 20)


# Binomial‐test filter: remove calls with VAF ≃0.5 (p > 1e-6)
filtered_table <- filtered_table %>%
  # for each row, compute two-sided binomial p-value against p=0.5
  rowwise() %>%
  mutate(
    p_value = binom.test(ALT_count, DP, p = 0.5)$p.value
  ) %>%
  ungroup()

# keep only those highly unlikely to be germline (AAF = 0.5)
filtered_table <- filtered_table %>%
  filter(p_value <= 1e-6)

# keep those with gnomAD frequency < 0.001 OR frequency = NA
filtered_table <- filtered_table %>%
  filter(is.na(population_frequency_gnomAD_manual) | population_frequency_gnomAD_manual < 0.001)

# no AAF cutoff, as this is a dilution sample used to estimate limit of detection

# rename to A117V_raw
A117V_raw <- filtered_table

# AAF in percent
A117V_raw$AAF_percent <- A117V_raw$AAF * 100

# A117V variant
A117V_raw$variant[A117V_raw$variant == "A117V_position"] <- "A117V"

# rename samples
A117V_raw$sample[A117V_raw$sample == "NA995A05_undil.noPoN"] <- "A117V 0.5% spike-in"
A117V_raw$sample[A117V_raw$sample == "NA99A1_undil.noPoN"] <- "A117V 1% spike-in"


#------------------------------------------------
# editing for final table
#------------------------------------------------

# only keep select columns, rename
A117V_table <- A117V_raw %>%
  select(
    "Sample" = sample,
    "Variant" = variant,
    "dbSNP ID" = dbsnp_id,
    "Chromosome" = CHROM,
    "Position" = POS,
    "REF" = REF,
    "ALT" = ALT,
    "Mutation Type" = mutation_type,
    "Read Depth (DP)" = DP,
    "REF Count" = REF_count,
    "ALT Count" = ALT_count,
    "AAF (%)" = AAF_percent,
    "Base Quality" = MEAN_BQ,
    "Mapping Quality" = MEAN_MQ,
    "REF Forward Count" = SB_refF,
    "REF Reverse Count" = SB_refR,
    "ALT Forward Count" = SB_altF,
    "ALT Reverse Count" = SB_altR
  )


# table in gt

A117V_gt <- A117V_table %>%
  gt() %>%
  tab_header(
    title = "A117V Spike-In Samples"
  ) %>%
  fmt_number(
    columns = c("AAF (%)", "Base Quality", "Mapping Quality"),
    decimals = 1
  ) %>%
  fmt_number(
    columns = c(
      "Read Depth (DP)", "REF Count", "ALT Count",
      "REF Forward Count", "REF Reverse Count",
      "ALT Forward Count", "ALT Reverse Count"
    ),
    decimals    = 0,
    use_seps    = TRUE
  ) %>%
  cols_align(
    align   = "center",
    columns = everything()
  ) %>%
  # remove all internal lines
  opt_table_lines(extent = "none") %>%
  tab_options(
    ## line under the title (midrule below heading)
    heading.border.bottom.width = px(1),
    heading.border.bottom.style = "solid",
    heading.border.bottom.color = "black",
    
    ## top rule (\toprule above column labels)
    table.border.top.width    = px(2),
    table.border.top.style    = "solid",
    table.border.top.color    = "black",
    
    ## midrule under column labels
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.style = "solid",
    column_labels.border.bottom.color = "black",
    
    ## bottom rule (\bottomrule)
    table.border.bottom.width = px(2),
    table.border.bottom.style = "solid",
    table.border.bottom.color = "black",
    
    ## strip any body hlines
    table_body.hlines.width   = px(0),
    
    ## layout tweaks
    table.font.size  = "small",
    data_row.padding = px(4),
    table.width      = pct(100)
  ) %>%
  opt_table_font(
    font = list(google_font("IBM Plex Sans"), default_fonts())
  )

# save as csv

setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-05-03_dilution_QC_analysis")

write.csv(A117V_table,
          file      = "A117V_table.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")

