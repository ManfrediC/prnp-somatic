
# this is the publication-ready version of the E200K LoD graph and analysis

library(tidyverse)
library(ggthemes)
library(cowplot)
library(magick)
library(ggpubr)
library(ggrepel)

#---------------------------------------------------------------------
# format
#---------------------------------------------------------------------

setwd("manuscript/legacy/figures/ddPCR_LoD/E200K_old")

mutation <- "E200K"

#list of samples (in desired order, left to right)
sample.list <- c("E200K_NTC", "E200K_WT", "E200K_0.05", "E200K_0.1", "E200K_0.33", "E200K_1")


#I changed the concentration column title from "Conc(copies/uL)" to "Conc" 
#and removed the space from "Fractional Abundance" to avoid problems

#import
data <- read.csv(paste0("Manfredi_",mutation,"_LOD.csv"), fileEncoding="UTF-8-BOM") %>%
  as_tibble()

#new labels
sample.labs <- c("NTC", "Non-mutant", paste0(mutation," 0.05%"), paste0(mutation," 0.1%"), 
                 paste0(mutation," 0.33%"), paste0(mutation," 1%"))
names(sample.labs) <- sample.list

#set order of samples or facets
data$Sample <- factor(data$Sample, levels = sample.list)

#---------------------------------------------------------------------
# Limit of Blank (LoB)
#---------------------------------------------------------------------

# blanks for E200K channel
blank_data <- data %>%
  filter(Sample %in% c("E200K_NTC", "E200K_WT"),
         Target == "E200K")

# pooled false-positive rate
p0 <- sum(blank_data$Positives) / sum(blank_data$Accepted.Droplets)

# calculate LoBFA per blank well
lobfa_tbl <- blank_data %>%
  rowwise() %>%
  mutate(
    LoB_count = qbinom(0.95, size = Accepted.Droplets, prob = p0),
    LoBFA     = LoB_count / Accepted.Droplets
  ) %>%
  ungroup()

p0

#---------------------------------------------------------------------
# local plate-specific LoBFA for E200K
#---------------------------------------------------------------------

# calculate LoBFA for each blank well under local p0
lobfa_tbl_local <- blank_data %>%
  rowwise() %>%
  mutate(
    LoB_count = qbinom(0.95, size = Accepted.Droplets, prob = p0),
    LoBFA     = LoB_count / Accepted.Droplets
  ) %>%
  ungroup() %>%
  select(Well, Sample, Accepted.Droplets, Positives, LoB_count, LoBFA)

# final LoBFA threshold = most conservative across blanks
lob_fa_local <- max(lobfa_tbl_local$LoBFA)

lob_fa_local

#---------------------------------------------------------------------
#concentration plot - 2025 version
#---------------------------------------------------------------------

concentration.plot

concentration.plot <- ggplot(data, mapping = aes(x = Target, y = Conc, color = Target)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = PoissonConfMin, ymax = PoissonConfMax), width = 0.5) + 
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + 
  annotation_logticks(
    sides = "l", outside = TRUE, short = unit(.5,"mm"), mid = unit(0.75,"mm"), long = unit(1.2,"mm")
  )  +
  coord_cartesian(clip = "off") +
  facet_grid(
    ~Sample,
    labeller = labeller(Sample = sample.labs)
  ) +
  theme_bw() +
  theme(panel.spacing = unit(0.75, "lines")) +
  theme(panel.grid.minor = element_blank()) +
  xlab("\nMosaic Samples") + ylab("Concentration (cpm)\n") + 
  geom_text(aes(label = signif(Conc, digits = 3)),
            size = 6*0.352777778, #0.35 is pt to mm conversion factor, I want font size 6
            hjust = if_else(data$Target == mutation, -0.4, 1.3), 
            vjust = if_else(data$Target == mutation, -1,
                            if_else(data$Sample == "NTC", -1, 2)),
            show.legend = FALSE)+
  theme(text = element_text(size=8)) + #overall text size 
  theme(strip.text.x = element_text(size = 6)) #facet label text size 

concentration.plot

#---------------------------------------------------------------------
# fraction plot - 2025 version
#---------------------------------------------------------------------

#the fraction is only given in the Target = mutation line
frac.subset <- subset(data, data$Target == mutation)
frac.subset$Sample <- factor(frac.subset$Sample, levels = sample.list)


frac.plot <- ggplot(frac.subset, 
                    mapping = aes(x = Sample, y = FractionalAbundance)) +
  geom_point() + 
  geom_errorbar(aes(ymin = PoissonFractionalAbundanceMin, ymax = PoissonFractionalAbundanceMax), width = 0.3) + 
  theme_bw() +
  theme(panel.spacing = unit(1, "lines")) +
  scale_x_discrete(labels = sample.labs) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 12)) +
  xlab("\nMosaic Samples") + ylab("Abundance of mutant allele (%)\n") + 
  geom_text(aes(label = signif(FractionalAbundance, digits = 2)), 
            size = 6*0.352777778, #0.35 is pt to mm conversion factor, I want font size 6 
            hjust= -1) +
  theme(text = element_text(size=8)) + #overall text size 
  theme(strip.text.x = element_text(size = 6)) + #facet label text size
  theme(plot.margin=grid::unit(c(1,0,0,0), "mm")) + #remove white space, but leave 1unit on top to avoid cutting off
  
  # add horizontal line for LoB (convert fraction -> percent)
  geom_hline(yintercept = lob_fa_local * 100, 
             linetype = "dashed", linewidth = 0.5, colour = "#0072B2")

frac.plot

#---------------------------------------------------------------------
# export plots
#---------------------------------------------------------------------

setwd("manuscript/figures/ddpcr_lod")

# concentration
ggsave(
  filename = "E200K_LoD_concentration.pdf",
  plot = concentration.plot,
  device = "pdf",
  width = 14, height = 7, units = "cm"
  )

# fraction
ggsave(
  filename = paste0(mutation,"_fraction.pdf"), 
  plot = frac.plot,
  width = 14, height = 7, units = "cm") #width of the conc. plot w/o legend in Inkscape

