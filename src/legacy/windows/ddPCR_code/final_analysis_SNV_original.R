#final analysis of SNV data

library(tidyverse)
library(ggthemes)
library(readxl)
library(cowplot)
library(grid)
library(ggpubr)
library(viridis)

#old
#setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/ddPCR/Final analysis/ddPCR-files/CSV files pool")

#snv.data <- read_excel("SNV_final_data_08-02-2021.xlsx")

setwd("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/ddPCR/Analysis/SNV graphs")

snv.data <- read_excel("SNV_data_final.xlsx")

mutation.list <- c("D178N", "E200K", "P102L")

#number of CJDs and controls
no.cjd <- 29
no.control <- 8

#loop through mutations
for (i in mutation.list) {

mutation <- i

##########
#format data
##########

#sample names
sample.names <- c("NTC", "WT", "MUT", paste0("CJD",1:no.cjd), paste0("Control",1:no.control))

sample.names

#shorter column names

oldnames = c("new_name", "brain_region",
             paste0(mutation,"_FractionalAbundance"),
             paste0(mutation,"_PoissonFractionalAbundanceMax"),
             paste0(mutation,"_PoissonFractionalAbundanceMin"))

newnames = c("patient","region", "abundance", "poissonmax", "poissonmin")

#subset data for single mutation
snv.subset <- snv.data[,colnames(snv.data) %in% oldnames] 

#rename columns
snv.subset <- snv.subset %>% 
  rename_at(vars(oldnames), ~ newnames)

#add levels to set order of samples in plot
snv.subset$patient <- factor(snv.subset$patient, levels = sample.names)

#subset without MUT and rows containing NA values
nomut <- subset(snv.subset, !(snv.subset$patient %in% c("WT","MUT","NTC"))) %>% 
  drop_na(abundance)

#colour-blind palette (with 6 colours) http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

#create label with mutation name
mylabel <- grobTree(textGrob(mutation, x=0.9,  y=0.85, hjust=0,
                             gp=gpar(fontsize=15)))


################
#create overview plot, assign to variable
################
assign(
  paste(mutation,"plot", sep = ""),

ggplot(nomut, mapping = aes(x = patient, y = abundance, colour = region)) + 
  geom_point(position = position_dodge(width = 0.6), na.rm = TRUE) +
  geom_errorbar(aes(ymin = poissonmax, ymax = poissonmin), width = 0, size = 0.3, position = position_dodge(width = 0.6)) +
  geom_hline(yintercept= if_else(mutation == "P102L", 0.1, 0.05), #0.1 cutoff for P102L, 0.05 for D178N and E200K
             linetype="dashed", color = "black", size = 0.5) +
  ylim(0,0.25)+
  xlab("\nPatient") + ylab(paste0("Frequency of ",mutation," allele (%)\n")) +
  theme_bw() +
  scale_colour_manual(name="Brain region",
                       labels = c("basal ganglia", "cerebellum", "frontal cortex", "hippocampus", "midbrain", "thalamus"), 
                       breaks = c("bg", "cb", "fr", "hc", "sn", "th"),
                       values = cbbPalette) + 
  annotation_custom(mylabel) +
  theme(text = element_text(size = 7)) +
  theme(legend.position=if_else(mutation == "P102L", "bottom", "none")) #only add legend under P102L

#end assign
)




################
#boxplots
################

#add patient category to dataframe

nomut$category <- if_else(nomut$patient %in% paste0("CJD",1:no.cjd), 
                          "CJD", "control")

#plot
assign(
  paste(mutation,"boxplot", sep = ""),


ggplot(nomut, mapping = aes(x = category, y = abundance, fill = category)) +
  stat_boxplot(geom ='errorbar', width = 0.2) + #add horizontal bars to boxplot whiskers
  geom_boxplot(outlier.shape = NA) + #as outliers are seen in jitter
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  scale_color_viridis(discrete = TRUE, alpha=0.6) + 
  ylab(paste0("frequency of ",mutation," allele (%)\n")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank()) +
  stat_compare_means(label.x = 1.3, 
                     label.y = 0.17)

)





#########################





################
#average MAF boxplots + csv
################

#select relevant rows and convert to wide format

widen <- select(nomut, patient, region, abundance) %>%
  pivot_wider(,names_from = region,
            values_from = abundance
            )

widen$category <- if_else(widen$patient %in% paste0("CJD",1:no.cjd), 
                          "CJD", "control")

widen$means <- rowMeans(widen[,2:7], na.rm = TRUE)

assign(
  paste(mutation,"means_boxplot", sep = ""),


ggplot(widen, mapping = aes(x = category, y = means, fill = category)) +
  stat_boxplot(geom ='errorbar', width = 0.2) + #add horizontal bars to boxplot whiskers
  geom_boxplot(outlier.shape = NA) + #as outliers are seen in jitter
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  scale_color_viridis(discrete = TRUE, alpha=0.6) + 
  ylab(paste0("frequency of ",mutation," allele (%)\n")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank()) +
  stat_compare_means(label.x = 1.3, 
                     label.y = 0.06)

)



write.table(widen, file = paste0("samples_mean_maf_",mutation,".csv"), 
            sep = ",", row.names = FALSE)

#end of for loop
}

################
#histogram
################

cjddata <- subset(nomut, category == "CJD")
ctrldata <- subset(nomut, category == "control")

ggplot(cjddata, aes(x = abundance)) +
  geom_histogram(binwidth = 0.002)

ggplot(ctrldata, aes(x = abundance)) +
  geom_histogram(binwidth = 0.0005)

###########
#save MAF plots
################


ggsave(
  filename = "D178N_ddPCRfrac.pdf", 
  plot = D178Nplot,
  width = 60, height = 22, units = "cm")

ggsave(
  filename = "E200K_ddPCRfrac.pdf", 
  plot = E200Kplot,
  width = 60, height = 22, units = "cm")

ggsave(
  filename = "P102L_ddPCRfrac.pdf", 
  plot = P102Lplot,
  width = 60, height = 22, units = "cm")

#make grid
grid <- plot_grid(
  D178Nplot, E200Kplot, P102Lplot,
  ncol = 1,
  nrow = 3, 
  rel_heights = c(1, 1, 1.428), 
  align= "hv",
  axis = "lr")

#height of bottom row is larger due to legend

grid

#############
#save boxplots
#############

ggsave(
  filename = "D178Nboxplot.pdf", 
  plot = D178Nboxplot,
  width = 10, height = 10, units = "cm")

ggsave(
  filename = "E200Kboxplot.pdf", 
  plot = E200Kboxplot,
  width = 10, height = 10, units = "cm")

ggsave(
  filename = "P102Lboxplot.pdf", 
  plot = P102Lboxplot,
  width = 10, height = 10, units = "cm")


boxgrid <- plot_grid(
  D178Nboxplot, E200Kboxplot, P102Lboxplot,
  labels = c("a","b","c"),
  ncol = 3,
  nrow = 1, 
  align= "hv",
  axis = "lr")

boxgrid


#############
#save means boxplots and csv
#############

ggsave(
  filename = "D178Nboxplot_means.pdf", 
  plot = D178Nmeans_boxplot,
  width = 10, height = 10, units = "cm")

ggsave(
  filename = "E200Kboxplot_means.pdf", 
  plot = E200Kmeans_boxplot,
  width = 10, height = 10, units = "cm")

ggsave(
  filename = "P102Lboxplot_means.pdf", 
  plot = P102Lmeans_boxplot,
  width = 10, height = 10, units = "cm")


boxgridmeans <- plot_grid(
  D178Nmeans_boxplot, E200Kmeans_boxplot, P102Lmeans_boxplot,
  labels = c("a","b","c"),
  ncol = 3,
  nrow = 1, 
  align= "hv",
  axis = "lr")

boxgridmeans


#############
#save grid plots in both directories
#############

directories = c("C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/ddPCR/Final analysis",
                "C:/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/Thesis/images")


ggsave(
  filename = "ddPCRgrid.pdf", 
  plot = grid,
  width = 26, height = 16, units = "cm")

ggsave(
  filename = "ddPCRboxplots.pdf", 
  plot = boxgrid,
  width = 26, height = 10, units = "cm")

ggsave(
  filename = "ddPCRboxplots_means.pdf", 
  plot = boxgridmeans,
  width = 26, height = 10, units = "cm")




#clean-up
#rm(list = ls())

