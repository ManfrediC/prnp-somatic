# --- pooled FA plot (one point per participant x assay) ---

library(tidyverse)
library(cowplot)
library(grid)
library(ggthemes)
library(ggpubr)

# paths
data_path <- "results/ddPCR/SNV_pooled_participant.xlsx"
out_dir   <- "manuscript/figures/ddpcr_fractional_abundance_pooled"

pooled_plot_data <- readxl::read_excel(data_path)

# mutations and LoD
mutation.list <- c("D178N", "E200K", "P102L")
lod_cut <- c(D178N = 0.056, E200K = 0.067, P102L = 0.13)

# initialise plot list
plots_pooled <- list()

# -------------------------------
# loop to create plots
# -------------------------------

for (mut in mutation.list) {
  
  # subset pooled data for this mutation
  df <- pooled_plot_data %>%
    filter(mutation == mut) %>%
    drop_na(fractional_abundance)
  
  # optional: remove heterozygous participant if applicable (example from before)
  # df <- df %>% filter(mutation != "E200K" | participant != "CJD30")
  
  # deterministic participant order: CJD1..CJDn then Control1..Controlm (keep only present)
  levels_order <- c(paste0("CJD", 1:200), paste0("Control", 1:200))
  levels_order <- intersect(levels_order, unique(df$participant))
  df <- df %>% mutate(participant = factor(participant, levels = levels_order))
  
  # corner label
  mylabel <- grobTree(textGrob(mut, x = 0.9, y = 0.85, hjust = 0, gp = gpar(fontsize = 15)))
  
  # plot: single colour; shape encodes Above LoB as before; no dodge needed
  p <- ggplot(
    df,
    aes(x = participant, y = fractional_abundance)
  ) +
    geom_point(aes(shape = detected_above_lob), size = 2, na.rm = TRUE) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, size = 0.3) +
    coord_cartesian(ylim = c(0, 0.25)) +
    xlab("\nParticipant") +
    ylab(paste0("Fractional abundance of ", mut, " allele\n")) +
    theme_bw() +
    # Above LoB legend identical to your regional plot
    scale_shape_manual(
      name   = "Above LoB",
      values = c(`FALSE` = 1, `TRUE` = 16),
      labels = c("No", "Yes")
    ) +
    annotation_custom(mylabel) +
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    geom_hline(yintercept = lod_cut[[mut]], linetype = "dashed", linewidth = 0.5)
  
  plots_pooled[[mut]] <- p
  
  ggsave(
    filename = file.path(out_dir, paste0("SNV_pooled_", mut, "_panel.pdf")),
    plot = p, width = 9, height = 4.5, dpi = 300
  )
}

# optional: combine the three pooled panels vertically
combined_pooled <- cowplot::plot_grid(plotlist = plots_pooled, ncol = 1, labels = "AUTO", align = "v")
ggsave(file.path(out_dir, "SNV_pooled_all_mutations.pdf"),
       combined_pooled, width = 10, height = 12, dpi = 300)
