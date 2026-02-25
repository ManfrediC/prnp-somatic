
# final analysis of SNV data

library(tidyverse)
library(readxl)
library(cowplot)
library(grid)
library(ggthemes)
library(ggpubr)
library(viridis)

# ---- paths ----
data_path <- "results/ddPCR/SNV_data_final.xlsx"
out_dir   <- "manuscript/figures/ddpcr_fractional_abundance"

snv.data <- readxl::read_excel(data_path)

# optional: exclude histological controls, keep only cases
# snv.data <- snv.data %>% filter(histotype != "control")

# colour-blind friendly palette for regions (same mapping as before)
cbbPalette <- c(
  bg = "#E69F00",  # basal ganglia
  cb = "#56B4E9",  # cerebellum
  fr = "#009E73",  # frontal cortex
  hc = "#0072B2",  # hippocampus
  sn = "#D55E00",  # midbrain (substantia nigra)
  th = "#CC79A7"   # thalamus
)

# nice labels for legend
region_labels <- c(
  bg = "basal ganglia",
  cb = "cerebellum",
  fr = "frontal cortex",
  hc = "hippocampus",
  sn = "midbrain",
  th = "thalamus"
)

# LoD cut-offs
lod_cut <- c(D178N = 0.05, E200K = 0.05, P102L = 0.10)

# ------------------------------------------------------
# mutation loop
# ------------------------------------------------------

mutation.list <- c("D178N", "E200K", "P102L")

# initialise plot list
plots <- list()

for (mut in mutation.list) {
  
# ------------------------------------------------------
# subset and light formatting
# ------------------------------------------------------

nomut <- snv.data %>%
  filter(mutation == mut) %>%
  drop_na(fractional_abundance) %>%
  mutate(
    # order patients: put any NTC/WT/mutant controls first if present, then others by code
    participant = factor(
      participant,
      levels = unique(participant)
    ),
    brain_region = factor(brain_region, levels = names(region_labels))
  )

# for E200K, we need to remove CJD30 as that sample is heterozygous
nomut <- nomut %>%
  filter(mutation != "E200K" | (mutation == "E200K" & participant != "CJD30"))

# ------------------------------------------------------
# order: CJD1..CJDn, then Control1..Controlm
# ------------------------------------------------------

# desired order
levels_order <- c(paste0("CJD", 1:31),
                  paste0("Control", 1:8))


# restrict to samples actually present in the data
levels_order <- intersect(levels_order, unique(nomut$participant))

# apply order
nomut$participant <- factor(nomut$participant, levels = levels_order)

# ------------------------------------------------------
# mylabel
# ------------------------------------------------------

mylabel <- grobTree(
  textGrob(mut, x = 0.9, y = 0.85, hjust = 0, gp = gpar(fontsize = 15))
)

# ------------------------------------------------------
# dodge
# ------------------------------------------------------

# define a single dodge to reuse everywhere
pd <- position_dodge(width = 0.6)

# ------------------------------------------------------
# the plot
# ------------------------------------------------------

# SNV plot
snv_plot <- ggplot(
  nomut,
  aes(x = participant,
      y = fractional_abundance,
      colour = brain_region)
) +
  
  # Points
  geom_point(
    aes(shape = detected_above_LoB, group = brain_region),
    position = pd, na.rm = TRUE
  ) +
  
  # Error bars
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high, group = brain_region),
    position = pd, width = 0, size = 0.3
  ) +
  
  # Define style of graph background
  coord_cartesian(ylim = c(0, 0.25)) +
  xlab("\nParticipant") +
  ylab(paste0("Fractional abundance of ", mut, " allele\n")) +
  theme_bw() +
  scale_colour_manual(
    name   = "Brain region",
    breaks = names(region_labels),
    labels = unname(region_labels),
    values = cbbPalette
  ) +

  # Define LoB legend
  scale_shape_manual(
    name   = "Above LoB",
    values = c(`FALSE` = 1, `TRUE` = 16),
    labels = c("No", "Yes")
  ) +
  
  # Add mutation in top right as lebel
  annotation_custom(mylabel) +
  
  # format legend
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  
  # add horizontal line for LoD
  geom_hline(yintercept = lod_cut[[mut]], linetype = "dashed", linewidth = 0.5) # add the LOD line

# ------------------------------------------------------
# save
# ------------------------------------------------------

# add plot to list (name it by mutation for clarity)
plots[[mut]] <- snv_plot

# save each panel
ggsave(
  filename = file.path(out_dir, paste0("SNV_", mut, "_panel.pdf")),
  plot = snv_plot, width = 9, height = 4.5, dpi = 300
)

# ------------------------------------------------------
# end of the loop
# ------------------------------------------------------

}

# ------------------------------------------------------
# combined row of all three
# ------------------------------------------------------

combined <- cowplot::plot_grid(plotlist = plots, ncol = 1, labels = "AUTO", align = "v")
ggsave(file.path(out_dir, "SNV_all_mutations.pdf"),
       combined, width = 10, height = 12, dpi = 300)