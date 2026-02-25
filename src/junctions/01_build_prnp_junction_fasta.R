#!/usr/bin/env Rscript

# 01_build_prnp_junction_fasta.R
# Build PRNP exon–exon junction reference sequences from Ensembl GTF + GRCh38 FASTA

library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(Biostrings)
library(GenomeInfoDb)

## ---- databases, main parameters ----
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
data_dir <- file.path(repo_root, "resources")

gtf_file <- file.path(data_dir, "Homo_sapiens.GRCh38.110.gtf.gz")
fasta_file <- file.path(data_dir, "hg38.fa")

# PRNP transcripts of interest
prnp_tx <- c("ENST00000379440",
             "ENST00000430350",
             "ENST00000457586")

# bases to take from each side of the junction
k <- 75L # as read length is about 150

# output
out_dir <- file.path(repo_root, "resources", "junctions")
out_fasta <- file.path(out_dir, "prnp_junctions.fa")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# represent the GTF chromosomes in UCSC style (we need to match Ensembl style to GATK/UCSC style)
txdb <- makeTxDbFromGFF(gtf_file)

## Make seqnames in the TxDb match UCSC-style (chr1, chr2, ...), as the hg38.fa uses chr-prefixed names.
seqlevelsStyle(txdb) <- "UCSC"

# Extract exons per transcript
ex_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)

# filter for PRNP (keep = Logical vector selecting our PRNP transcripts)
keep <- names(ex_by_tx) %in% prnp_tx
if (!any(keep)) stop("None of the requested PRNP transcript IDs were found in the GTF.")

# Reduce to only PRNP transcripts of interest
ex_by_tx <- ex_by_tx[keep]

# open the hg.38 FASTA
fa <- FaFile(fasta_file)
open(fa)

# Prepare containers for results
junction_seqs  <- DNAStringSet() # will store all junction sequences
junction_names <- character()    # will store the FASTA headers for each sequence

# Loop over transcripts and exon pairs
for (tx_id in names(ex_by_tx)) {
  exons <- sort(ex_by_tx[[tx_id]]) # Retrieves exons for the transcript and sorts them by genomic coordinate
  if (length(exons) < 2L) next     # If the transcript has fewer than 2 exons, there can be no exon–exon junctions, so next.
  
  for (i in seq_len(length(exons) - 1L)) {
    upstream   <- exons[i]
    downstream <- exons[i + 1L]
    
    # Take k bases from each side of the junction
    up_k <- resize(upstream,   width = min(width(upstream),   k), fix = "end")   
    dn_k <- resize(downstream, width = min(width(downstream), k), fix = "start")
    
    # Fetch sequences from the FASTA
    up_seq <- getSeq(fa, up_k)
    dn_seq <- getSeq(fa, dn_k)
    
    # Construct the junction sequence by concatenating the two DNA strings end-to-end
    junc_seq <- xscat(up_seq, dn_seq)
    
    # Build name for junction
    junc_name <- paste("PRNP", tx_id,
                       paste0("exon", i, "_exon", i + 1L),
                       sep = "|")
    
    # store in the containers
    junction_seqs  <- append(junction_seqs, junc_seq)
    junction_names <- c(junction_names, junc_name)
  }
}

# close FASTA
close(fa)

# sanity check (did it actually work)
if (length(junction_seqs) == 0L) {
  stop("No junction sequences were generated (check PRNP transcript IDs and GTF).")
}

# Assign names and write FASTA
names(junction_seqs) <- junction_names
writeXStringSet(junction_seqs, filepath = out_fasta)

message("Wrote ", length(junction_seqs), " junction sequences to: ", out_fasta)
