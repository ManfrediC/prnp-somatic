library(tidyverse)
library(ggthemes)
library(cowplot)
library(magick)
library(ggpubr)
library(ggrepel)

#---------------------------------------------------------------------
# Setup
#---------------------------------------------------------------------
setwd("manuscript/figures/ddpcr_lod/e200k_old")

mutation.list <- c("D178N", "E200K", "P102L")
# mutation.list <- c("D178N")
# mut <- "D178N"

#---------------------------------------------------------------------
# Helper: run the full LoD workflow for one mutation
#---------------------------------------------------------------------
make_lod_plots <- function(mut) {
  # Sample IDs (left-to-right) and display labels
  sample.list <- c(paste0(mut, "_NTC"),
                   paste0(mut, "_WT"),
                   paste0(mut, "_0.05"),
                   paste0(mut, "_0.1"),
                   paste0(mut, "_0.33"),
                   paste0(mut, "_1"))
  
  sample.labs <- c("NTC",
                   "Non-mutant",
                   paste0(mut, " 0.05%"),
                   paste0(mut, " 0.1%"),
                   paste0(mut, " 0.33%"),
                   paste0(mut, " 1%"))
  names(sample.labs) <- sample.list
  
  # Import
  infile <- paste0(mut, "_LOD.csv")
  data <- read.csv(infile, fileEncoding = "UTF-8-BOM") %>% as_tibble()
  data$Sample <- factor(data$Sample, levels = sample.list)
  
  #-------------------------------------------------------------------
  # LoB: pooled p0 from the two blanks on the mutation channel
  #-------------------------------------------------------------------
  blank_data <- data %>%
    filter(Sample %in% c(paste0(mut, "_NTC"), paste0(mut, "_WT")),
           Target == mut)
  
  p0 <- sum(blank_data$Positives) / sum(blank_data$Accepted.Droplets)
  
  # Per-blank LoB as fractional abundance; take max (most conservative)
  lobfa_tbl_local <- blank_data %>%
    rowwise() %>%
    mutate(
      LoB_count = qbinom(0.95, size = Accepted.Droplets, prob = p0),
      LoBFA     = LoB_count / Accepted.Droplets
    ) %>%
    ungroup() %>%
    select(Well, Sample, Accepted.Droplets, Positives, LoB_count, LoBFA)
  
  lob_fa_local <- min(lobfa_tbl_local$LoBFA)
  print(lob_fa_local)
  print(lobfa_tbl_local)
  
  #-------------------------------------------------------------------
  # Concentration plot (same style; precompute text placement columns)
  #-------------------------------------------------------------------
  data <- data %>%
    mutate(
      is_ntc    = grepl("_NTC$", Sample),
      txt_hjust = ifelse(Target == mut, -0.4, 1.3),
      txt_vjust = ifelse(
        Target == mut, -1,
        ifelse(is_ntc, -1, 2)
      )
    )

  
  # # Replace zeros with small pseudo-value (0.00001%)
  # eps <- 1e-5
  # data.conc <- data %>%
  #   mutate(
  #     Conc = ifelse(Conc == 0, eps, Conc),
  #     TotalConfMin = ifelse(TotalConfMin == 0, eps, TotalConfMin),
  #     TotalConfMax = ifelse(TotalConfMax == 0, eps, TotalConfMax),
  #     PoissonConfMin = ifelse(PoissonConfMin == 0, eps, PoissonConfMin),
  #     PoissonConfMax = ifelse(PoissonConfMax == 0, eps, PoissonConfMax)
  #   )
  
  # plot
  concentration.plot <- ggplot(data, aes(x = Target, y = Conc, colour = Target)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = PoissonConfMin, ymax = PoissonConfMax), width = 0.5) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    annotation_logticks(
      sides = "l", outside = TRUE,
      short = unit(.5, "mm"), mid = unit(0.75, "mm"), long = unit(1.2, "mm")
    ) +
    coord_cartesian(clip = "off") +
    facet_grid(~Sample, labeller = labeller(Sample = sample.labs)) +
    theme_bw() +
    theme(panel.spacing = unit(0.75, "lines"),
          panel.grid.minor = element_blank()) +
    xlab("\nMosaic Samples") + ylab("Concentration (cpm)\n") +
    geom_text(aes(label = signif(Conc, digits = 3),
                  hjust = txt_hjust,
                  vjust = txt_vjust),
              size = 6 * 0.352777778,
              show.legend = FALSE) +
    theme(text = element_text(size = 8),
          strip.text.x = element_text(size = 6))
  
  #-------------------------------------------------------------------
  # Fraction plot (mutation channel only, exclude NTC as it's NA) + LoB line
  #-------------------------------------------------------------------
  frac.subset <- data %>%
    filter(Target == mut, Sample != paste0(mut, "_NTC")) %>%
    mutate(Sample = factor(Sample, levels = sample.list))
  
  # plot
  frac.plot <- ggplot(frac.subset, aes(x = Sample, y = FractionalAbundance)) +
    geom_point() +
    
    # error bars with Poisson CI
    geom_errorbar(aes(ymin = PoissonFractionalAbundanceMin,
                      ymax = PoissonFractionalAbundanceMax),
                  width = 0.3) +
    
    # base theme
    theme_bw(base_size = 8, base_family = "Helvetica") +   # global font
    theme(
      panel.spacing = unit(1, "lines"),
      panel.grid.minor = element_blank(),
      
      # axis text and titles
      axis.text.x = element_text(size = 8, vjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 8, face = "plain"),
      axis.title.y = element_text(size = 8, face = "plain"),
      
      # facet strips
      strip.text.x = element_text(size = 6),
      
      # plot margin
      plot.margin = grid::unit(c(1, 1, 1, 1), "mm")
    ) +
    
    # axis scales and custom labels
    scale_x_discrete(labels = sample.labs) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 12)) +
    
    # axis titles
    xlab("\nMosaic Samples") +
    ylab("Abundance of mutant allele (%)\n") +
    
    # FA values near points
    geom_text(aes(label = signif(FractionalAbundance, digits = 2)),
              size = 8 * 0.352777778,
              hjust = -0.75,
              vjust = -0.75) +
    
    # LoB line
    geom_hline(yintercept = lob_fa_local * 100,
               linetype = "dashed", linewidth = 0.6, colour = "#0072B2") +
    
    # LoB annotation
    annotate("text",
             x = Inf, y = lob_fa_local * 100,
             label = "LoB",
             hjust = 1.1, vjust = -0.3,
             size = 8 * 0.352777778,
             colour = "#0072B2")
  
  #-------------------------------------------------------------------
  # Export
  #-------------------------------------------------------------------
  ggsave(filename = paste0(mut, "_LoD_concentration.pdf"),
         plot = concentration.plot, device = "pdf",
         width = 14, height = 7, units = "cm")
  
  ggsave(filename = paste0(mut, "_fraction.pdf"),
         plot = frac.plot,
         width = 14, height = 7, units = "cm")
  
  invisible(list(concentration = concentration.plot,
                 fraction = frac.plot,
                 lob_fa = lob_fa_local,
                 p0 = p0))
}

#---------------------------------------------------------------------
# Run for all mutations
#---------------------------------------------------------------------
results <- lapply(mutation.list, make_lod_plots)

