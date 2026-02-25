# this script builds a BED file with padded


#!/usr/bin/env Rscript
library(GenomicFeatures)
library(GenomeInfoDb)
library(rtracklayer)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "resources")
output_dir <- file.path(repo_root, "resources", "junctions")

gtf <- file.path(input_dir, "Homo_sapiens.GRCh38.110.gtf.gz")
bed <- file.path(output_dir, "PRNP.pad1kb.hg38.bed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------
# Derive PRNP gene_id from the GTF (should be ENSG00000171867)
# ------------------------------------

# Import only gene entries from the GTF
gtf_gr <- import(gtf)

mask <- gtf_gr$type == "gene" &
  !is.na(gtf_gr$gene_name) &
  gtf_gr$gene_name == "PRNP"

prnp_gene_rows <- gtf_gr[mask]
prnp_gene_rows
length(prnp_gene_rows)

if (length(prnp_gene_rows) != 1L) {
  stop("PRNP gene not found uniquely in GTF (found ", length(prnp_gene_rows), " entries).")
}

prnp_gene_id <- prnp_gene_rows$gene_id
prnp_gene_id

# ------------------------------------
# match PRNP
# ------------------------------------

# Build transcript database and harmonise chromosome naming
txdb <- makeTxDbFromGFF(gtf)
seqlevelsStyle(txdb) <- "UCSC"   # match hg38.fa


# match PRNP gene with the gene ID
g_all <- genes(txdb)
g <- g_all[prnp_gene_id]

# Pad the gene region by 1000bp
pad <- 1000
gr  <- GenomicRanges::trim(
  GenomicRanges::resize(g, width = width(g) + 2*pad, fix = "center")
)

# export to BED
export(gr, con = bed, format = "bed")
cat("Wrote:", bed, "\n")
