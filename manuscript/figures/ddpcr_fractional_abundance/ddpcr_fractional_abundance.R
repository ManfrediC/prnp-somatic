
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
lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)

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
# count samples
# ------------------------------------------------------

# total number of brain samples (without the single pons sample)
samples <- snv.data %>%
  distinct(participant, histotype, brain_region) %>%
  filter(brain_region != "ps")

# individuals
individuals <- snv.data$participant %>%
  unique()

# per brain region
table(samples$brain_region)

# per histological subtype
histo <- samples %>%
  distinct(participant, histotype)

table(histo$histotype)


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

# tall
# ggsave( filename = file.path(out_dir, paste0("SNV_", mut, "_panel.pdf")), 
#         plot = snv_plot, 
#         width = 9, 
#         height = 4.5, 
#         dpi = 300 )

# short
ggsave(
  filename = file.path(out_dir, paste0("SNV_", mut, "_panel.pdf")),
  plot = snv_plot + theme(plot.margin = margin(5.5, 20, 5.5, 5.5, "pt")),
  width = 9,
  height = 3,
  dpi = 300
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


# ------------------------------------------------------
# combined row of all three with legend at bottom
# ------------------------------------------------------

# Build one plot with bottom legend to extract
p_for_legend <- snv_plot +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(face = "bold", size = 9),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(4, "mm")) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3)),
         shape  = guide_legend(nrow = 1))

leg <- cowplot::get_legend(p_for_legend)

# Remove legends from individual panels
plots_noleg <- lapply(plots, function(p) p + theme(legend.position = "none"))

# Stack panels and add the single legend at the bottom
combined <- cowplot::plot_grid(plotlist = plots_noleg, ncol = 1, labels = , align = "v")
final_with_legend_bottom <- cowplot::plot_grid(combined, leg, ncol = 1, rel_heights = c(1, 0.08))

cowplot::save_plot(
  file.path(out_dir, "SNV_all_mutations_legend_bottom.pdf"),
  final_with_legend_bottom,
  base_width = 10, base_height = 8  # shorter overall
)