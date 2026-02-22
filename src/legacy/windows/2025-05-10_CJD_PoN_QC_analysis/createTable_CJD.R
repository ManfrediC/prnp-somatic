library(tidyverse)
library(gt)
library(pagedown)

# ─────────────────────────────────────────────────────────────
# 1.  INPUT ─ read all *.withPoN.tsv files
# ─────────────────────────────────────────────────────────────
setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-05-10_CJD_PoN_QC_analysis/variant_tables")

file_list <- list.files(pattern = "\\.withPoN\\.tsv$")

# drop any pre-existing summary/helper files
file_list <- setdiff(
  file_list,
  c("summary_combined_variants.tsv", "withPoN_PRNP_PASS.tsv",
    "withPoN_TET2_PASS.tsv",      "withPoN_TTN_PASS.tsv",
    "final_withPoN_variants.tsv")
)

summary_table <- NULL

for (file in file_list) {
  temp_table <- read_tsv(file,
                         col_names = TRUE,
                         col_types = cols(.default = col_character()))
  
  # Standardise FORMAT columns (GT, DP, AD, …)
  temp_table <- temp_table %>%
    rename_with(~ str_remove(., "^.*\\."), ends_with(c(".GT", ".DP", ".AD",
                                                       ".F1R2", ".F2R1", ".SB")))
  
  if (nrow(temp_table) == 0) next    # skip empty files
  
  sample_name <- str_remove(file, "\\.withPoN\\.tsv$")
  colnames(temp_table) <- str_replace(colnames(temp_table),
                                      paste0("^", sample_name, "\\."),
                                      "sample.")
  temp_table$sample <- sample_name
  temp_table <- select(temp_table, sample, everything())
  
  summary_table <- bind_rows(summary_table, temp_table)
}
rm(file, file_list, sample_name, temp_table)

# ─────────────────────────────────────────────────────────────
# 2.  REF AND ALT COUNTS
# ─────────────────────────────────────────────────────────────
summary_table <- summary_table %>%
  separate(AD, into = c("REF_count", "ALT_count"), sep = ",", convert = TRUE) %>%
  mutate(DP  = as.numeric(DP),
         AAF = ALT_count / DP)

# ─────────────────────────────────────────────────────────────
# 3.  FUNCOTATION PARSING
# ─────────────────────────────────────────────────────────────
summary_table <- summary_table %>%
  mutate(
    gene               = str_extract(FUNCOTATION, "(?<=\\[)[^|]+"),
    location_relative  = word(FUNCOTATION, 6, sep = fixed("|")),
    dbsnp_id           = str_extract(FUNCOTATION, "rs\\d+"),
    mutation_type      = str_extract(
      FUNCOTATION,
      "\\b(SILENT|INTRON|MISSENSE|IGR|INTERGENIC|UPSTREAM|DOWNSTREAM|FIVE_PRIME_UTR|
          THREE_PRIME_UTR|SPLICE_SITE|SPLICE_DONOR|SPLICE_ACCEPTOR|FRAME_SHIFT|
          STOP_GAINED|STOP_LOST|START_LOST|START_GAINED|SYNONYMOUS|NONSENSE|
          INFRAME_DELETION|INFRAME_INSERTION|COMPLEX_SUBSTITUTION|TRANSCRIPT_ABLATION|
          TRANSCRIPT_AMPLIFICATION|REGULATORY_REGION_VARIANT|
          TF_BINDING_SITE_VARIANT)\\b"),
    variant = case_when(                         # Pathogenic PRNP sites
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
    )
  ) %>%
  relocate(variant, .after = sample) %>%
  filter(gene %in% c("PRNP", "TTN", "TET2"))

# POS as numeric
summary_table <- summary_table %>%
  mutate(POS = as.numeric(POS))

# PRNP intron / coding region classification
summary_table <- summary_table %>%
  mutate(region = case_when(
    CHROM == "chr20" & POS <= 4686455                       ~ "5' upstream",
    CHROM == "chr20" & between(POS, 4686456, 4686512)       ~ "exon 1",
    CHROM == "chr20" & between(POS, 4686513, 4699210)       ~ "intron",
    CHROM == "chr20" & between(POS, 4699221, 4699982)       ~ "protein coding",
    CHROM == "chr20" & between(POS, 4699983, 4701588)       ~ "exon 2 downstream of ORF",
    CHROM == "chr20" & POS >= 4701589                       ~ "3' downstream",
    TRUE                                                    ~ NA_character_
  )) %>%
  relocate(region, .after = CHROM)

# ─────────────────────────────────────────────────────────────
# 4.  POPULATION-FREQUENCY ANNOTATIONS
# ─────────────────────────────────────────────────────────────
summary_table <- summary_table %>%
  mutate(population_frequency_gnomAD_manual = case_when(
    dbsnp_id == "rs1012397146" ~ 0.00001318,
    dbsnp_id == "rs1931969"    ~ 0.2146,
    dbsnp_id == "rs2093390"    ~ 0.3055,
    dbsnp_id == "rs2093391"    ~ 0.2583,
    dbsnp_id == "rs6037932"    ~ 0.2872,
    dbsnp_id == "rs6052769"    ~ 0.4115,
    dbsnp_id == "rs6052772"    ~ 0.2361,
    dbsnp_id == "rs6107516"    ~ 0.2271,
    dbsnp_id == "rs781755787"  ~ 0.000001369,
    TRUE ~ NA_real_
  )) %>%
  mutate(dbsnp_id = if_else(dbsnp_id == "rs996098774", NA_character_, dbsnp_id),
         variant_id = str_c(str_remove(CHROM, "^chr"), POS, REF, ALT, sep = "-")) %>%
  mutate(population_frequency_gnomAD_manual =
           if_else(variant_id %in%
                     c("2-178581749-C-T", "2-178585952-T-A", "2-178585984-A-T",
                       "2-178649895-C-T", "2-178649897-T-A", "2-178698399-T-A",
                       "2-178698412-G-A", "20-4691495-C-G", "20-4691920-G-A",
                       "20-4693834-CA-AG", "20-4694249-T-C", "20-4694270-A-G",
                       "20-4695237-T-C", "20-4695357-G-C", "4-105175251-A-G",
                       "4-105175261-T-G", "4-105215846-T-G", "4-105261420-A-T",
                       "4-105273277-G-A"),
                   NA_real_,
                   population_frequency_gnomAD_manual))

# ─────────────────────────────────────────────────────────────
# 5.  BASE/READ-QUALITY METRICS (from bam-readcount)
# ─────────────────────────────────────────────────────────────
setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-05-10_CJD_PoN_QC_analysis/output_CSV")

sample_basecount <- read_csv("readcounts_samples.csv", show_col_types = FALSE) %>%
  rename(ALT = BASE, sample_name = Sample) %>%      # harmonise names
  filter(COUNT != 0)                                # drop empty rows

### harmonise dataframes

# POS as numeric
sample_basecount <- sample_basecount %>%
  mutate(POS = as.numeric(POS))
summary_table <- summary_table %>%
  mutate(POS = as.numeric(POS))

# sample_name column
summary_table <- summary_table %>%
  mutate(sample_name = sample)

# REF and ALT must be all capitals
sample_basecount <- sample_basecount %>%
  mutate(REF = toupper(REF),
         ALT = toupper(ALT))

# Join MEAN_BQ and MEAN_MQ onto summary_table
summary_table <- summary_table %>%
  left_join(
    sample_basecount %>%
      select(sample_name, CHROM, POS, REF, ALT, MEAN_BQ, MEAN_MQ),
    by = c("sample_name", "CHROM", "POS", "REF", "ALT")
  )

# ─────────────────────────────────────────────────────────────
# 6.  QUALITY-CONTROL FILTERS
# ─────────────────────────────────────────────────────────────

# 6.1  Split strand-bias field into separate counts
qc_table <- summary_table %>% 
  separate(SB,
           into = c("SB_refF", "SB_refR", "SB_altF", "SB_altR"),
           sep   = ",",
           convert = TRUE)

# 6.2  Depth filter  (ALT ≥ 10; total depth ≥ 100)
qc_table <- qc_table %>% 
  filter(ALT_count >= 10,
         DP        >= 100)

# 6.3  Strand-balance filter  (≥ 3 ALT reads on each strand)
qc_table <- qc_table %>% 
  filter(SB_altF >= 3,
         SB_altR >= 3)

# 6.4  Base- and mapping-quality filter  (mean ≥ 20)
qc_table <- qc_table %>% 
  filter(MEAN_BQ >= 20,
         MEAN_MQ >= 20)

# 6.5  Binomial test against germline VAF ≈ 0.5
qc_table <- qc_table %>% 
  rowwise() %>% 
  mutate(p_value = binom.test(ALT_count, DP, p = 0.5)$p.value) %>% 
  ungroup() %>%
  filter(p_value <= 1e-6)

# 6.6  Population-frequency filter  (gnomAD < 0.001 or NA)
filtered_table <- qc_table %>% 
  filter(is.na(population_frequency_gnomAD_manual) |
           population_frequency_gnomAD_manual < 0.001)


# ─────────────────────────────────────────────────────────────
# 7.  OUTPUT
# ─────────────────────────────────────────────────────────────
setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2025-05-10_CJD_PoN_QC_analysis/output_CSV")

write_csv(filtered_table, "filtered_withPoN_variants.csv")
cat("\nFiltered variant table written to filtered_withPoN_variants.csv\n")

# ─────────────────────────────────────────────────────────────
# 8.  editing for final table
# ─────────────────────────────────────────────────────────────

# only keep select columns, rename
CJD_table <- filtered_table %>%
  select(
    "Sample" = sample,
    "Chromosome" = CHROM,
    "Position" = POS,
    "REF" = REF,
    "ALT" = ALT,
    "Mutation Type" = mutation_type,
    "Population Frequency" = population_frequency_gnomAD_manual,
    "Read Depth (DP)" = DP,
    "REF Count" = REF_count,
    "ALT Count" = ALT_count,
    "AAF (%)" = AAF,
    "Base Quality" = MEAN_BQ,
    "Mapping Quality" = MEAN_MQ,
    "REF Forward Count" = SB_refF,
    "REF Reverse Count" = SB_refR,
    "ALT Forward Count" = SB_altF,
    "ALT Reverse Count" = SB_altR
  )

# AAF in %
CJD_table <- CJD_table %>%
  mutate(`AAF (%)` = `AAF (%)` * 100)
CJD_table$`AAF (%)` <- round(CJD_table$`AAF (%)`, 2)

# if Population frequency is NA, set to "not described"
CJD_table$`Population Frequency` <- ifelse(
  is.na(CJD_table$`Population Frequency`),
  "not described",
  CJD_table$`Population Frequency`
)

# write table to CSV
write.csv(CJD_table,
          file      = "CJD_variants.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")
