library(GenomicRanges)
library(trackViewer)
library(grid)

# Define SNVs
positions <- c(4691920, 4691920, 4691920, 4694249)
samples <- c("CJD2", "CJD6", "CJD23", "CJD23")
AAFs <- c(0.55, 0.52, 0.6, 0.82)
REFs <- c("G", "G", "G", "T")
ALTs <- c("A", "A", "A", "C")

# Create GRanges for SNVs
sample.gr <- GRanges(
  seqnames = "chr20",
  ranges = IRanges(start = positions, width = 1, names = paste0(samples, "_", positions)),
  strand = "*"
)

# Add metadata
mcols(sample.gr)$sample <- samples
mcols(sample.gr)$AAF <- AAFs
mcols(sample.gr)$REF <- REFs
mcols(sample.gr)$ALT <- ALTs

# Define PRNP exons
exon_starts <- c(4686456, 4699211)
exon_ends <- c(4686512, 4701588)
features <- GRanges(
  seqnames = "chr20",
  ranges = IRanges(start = exon_starts, end = exon_ends, names = c("Exon1", "Exon2")),
  strand = "+"
)

# exons: colours and height of box
features$fill <- c("#CC79A7", "#51C6E6")
features$height <- 0.07

# colours for samples
sample_colours <- c(
  "CJD2" = "#0072B2",   # blue
  "CJD6" = "#FF8833",   # orange
  "CJD23" = "#009E73"   # blue-green
)

# Assign colours based on sample
sample.gr$color <- sample_colours[sample.gr$sample]

# labels
sample.gr$label.parameter.label <- paste0(
  start(sample.gr), " ", sample.gr$REF, ">", sample.gr$ALT)

# reduce alpha for 4691920 (uncertain)
sample.gr$alpha <- 1
idx <- start(sample.gr) == 4691920
sample.gr$alpha[idx] <- 0.3

# scale length of stalks by factor
scaleFactor <- 10
sample.gr$score <- AAFs * scaleFactor

# build a y-axis: ticks every 1 stem-unit, labels as true AAF
ticks <- seq(0, scaleFactor, by = 1)
names(ticks) <- sprintf("%.2f", ticks / scaleFactor)

# ---- NEW: define x-axis window and ticks ----
view <- GRanges("chr20", IRanges(start = 4686000L, end = 4702000L))

xticks <- seq(4686000L, 4702000L, by = 2000L)
names(xticks) <- format(xticks, big.mark = "'", scientific = FALSE)
# --------------------------------------------

# plot on screen
grid::grid.newpage()
lolliplot(
  sample.gr,
  features,
  ranges = view,     # <- sets x-axis span
  xaxis  = xticks,   # <- 2 kb ticks
  yaxis  = ticks,
  ylab   = "",
  newpage = TRUE
)

# export to PDF
setwd("manuscript/figures/snv_lollipop")

svg("SNV_lollipop.svg", width = 9, height = 4.5, bg = "transparent")
lolliplot(
  sample.gr,
  features,
  ranges = view,
  xaxis  = xticks,
  yaxis  = ticks,
  ylab   = "",
  newpage = TRUE
)
dev.off()


