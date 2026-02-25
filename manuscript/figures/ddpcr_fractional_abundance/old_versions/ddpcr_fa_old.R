
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

mutation.list <- c("D178N", "E200K", "P102L")

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

# ------------------------------------------------
# loop per mutation
# ------------------------------------------------

plots <- vector("list", length(mutation.list))
names(plots) <- mutation.list

for (i in seq_along(mutation.list)) {
  mut <- mutation.list[i]
  
  # subset and light formatting
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
  
  # order: CJD1..CJDn, then Control1..Controlm
  cjd_levels <- nomut %>%
    filter(str_detect(participant, "^CJD\\d+$")) %>%
    mutate(n = as.integer(str_extract(participant, "\\d+"))) %>%
    arrange(n) %>%
    pull(participant) %>%
    unique()
  
  ctrl_levels <- nomut %>%
    filter(str_detect(participant, "^Control\\d+$")) %>%
    mutate(n = as.integer(str_extract(participant, "\\d+"))) %>%
    arrange(n) %>%
    pull(participant) %>%
    unique()
  
  level_order <- c(cjd_levels, ctrl_levels)
  
  nomut <- nomut %>%
    mutate(participant = factor(participant, levels = level_order))
  
  # label grob
  mylabel <- grobTree(
    textGrob(mut, x = 0.9, y = 0.85, hjust = 0, gp = gpar(fontsize = 15))
  )
  
  # define a single dodge to reuse everywhere
  pos <- position_dodge2(width = 0.6, preserve = "single")
  
  # --------------------------------
  # build plot
  # --------------------------------
  
  p <- ggplot(
    nomut,
    aes(x = participant,
        y = fractional_abundance,
        colour = brain_region)
  ) +
    # POINTS
    geom_point(
      aes(shape = detected_above_lob, group = brain_region),
      position = pos, na.rm = TRUE
    ) +
    # ERROR BARS
    geom_errorbar(
      aes(ymin = ci_low, ymax = ci_high, group = brain_region),
      width = 0, size = 0.3, position = pos
    ) +
    # LOB TICK
    geom_segment(
      aes(x = as.numeric(participant) - 0.2,
          xend = as.numeric(participant) + 0.2,
          y = lob_fa, yend = lob_fa),
      inherit.aes = FALSE,
      linewidth = 0.25,
      alpha = 0.4
    ) +
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
    scale_shape_manual(
      name   = "Above LoB",
      values = c(`FALSE` = 1, `TRUE` = 16),
      labels = c("No", "Yes")
    ) +
    annotation_custom(mylabel) +
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    geom_hline(yintercept = lod_cut[[mut]], linetype = "dashed", linewidth = 0.5) # add the LOD line
  
  plots[[i]] <- p
  
  # save each panel
  ggsave(
    filename = file.path(out_dir, paste0("SNV_", mut, "_panel.png")),
    plot = p, width = 9, height = 4.5, dpi = 300
  )
}

# optional: combined row of all three
combined <- cowplot::plot_grid(plotlist = plots, ncol = 1, labels = "AUTO", align = "v")
ggsave(file.path(out_dir, "SNV_all_mutations.png"), combined, width = 10, height = 12, dpi = 300)








