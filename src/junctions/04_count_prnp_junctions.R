#!/usr/bin/env Rscript

library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomeInfoDb)
library(Biostrings)
library(dplyr)


## ---- configure (adjust to your layout) ----
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
legacy_dir <- file.path(repo_root, "prnp-junctions")

default_gtf <- file.path(repo_root, "resources", "Homo_sapiens.GRCh38.110.gtf.gz")
if (!file.exists(default_gtf)) {
  default_gtf <- file.path(legacy_dir, "data", "Homo_sapiens.GRCh38.110.gtf.gz")
}

GTF      <- Sys.getenv("PRNP_JUNCTION_GTF", unset = default_gtf)
JUNC_DIR <- Sys.getenv("PRNP_JUNCTION_ALIGN_DIR", unset = file.path(repo_root, "results", "junctions", "junction_align"))
if (!file.exists(GTF)) stop("GTF not found: ", GTF)
if (!dir.exists(JUNC_DIR)) stop("Junction alignment directory not found: ", JUNC_DIR)

# Analysis parameters (must match the junction FASTA builder)
K          <- 75L  # number of bases per side used when building junction FASTA
MIN_OVER   <- 10L  # minimal overhang required on each side of the junction
MIN_MAPQ   <- 20L  # minimal mapping quality
MIN_FRAGS  <- 5L   # minimal number of junction fragments in a sample to consider it potentially positive (maybe adapt)
MIN_STARTS <- 3L   # minimal number of distinct start positions to avoid one PCR clone dominating
PRNP_TX    <- c("ENST00000379440","ENST00000430350","ENST00000457586")  # transcript IDs for PRNP splicing isoforms

## ---- build TxDb and junction lookup (join position per contig) ----
txdb <- makeTxDbFromGFF(GTF) # Reads the GTF and builds a TxDb object: a structured store of genes, transcripts, exons, etc
seqlevelsStyle(txdb) <- "UCSC"

ex_by_tx <- exonsBy(txdb, by="tx", use.names=TRUE)[PRNP_TX]
ex_by_tx <- ex_by_tx[!vapply(ex_by_tx, is.null, logical(1))]

# join_tbl = junction metadata per contig: contig name (rname) → where the junction is (join_bp)
join_tbl <- bind_rows(lapply(names(ex_by_tx), function(tx) { # Loop over each transcript ID (tx)
  exs <- GenomicRanges::sort(ex_by_tx[[tx]])                 # Get the exons for that transcript and sort them in genomic order 
  if (length(exs) < 2L) return(NULL)                         # If fewer than two exons, there is no exon–exon junction, skip
  do.call(rbind, lapply(seq_len(length(exs)-1L), function(i) {
    up <- min(width(exs[i]),      K)                         # minimum length of exons is K
    dn <- min(width(exs[i + 1L]), K)
    data.frame(
      rname   = paste("PRNP", tx, sprintf("exon%d_exon%d", i, i+1L), sep="|"), # rname = name of the junction contig
      join_bp = as.integer(up),                              # position (1-based) within the contig where the exon–exon boundary lies
      contigL = as.integer(up + dn),                         # total length of contig (up + dn)
      stringsAsFactors = FALSE
    )
  }))
}))
if (nrow(join_tbl) == 0L) stop("No PRNP junctions derived from the GTF.")

## ---- enumerate junction BAMs ----
bam_files <- list.files(JUNC_DIR, pattern="\\.PRNP\\.toJunc\\.bam$", full.names=TRUE)
if (!length(bam_files)) stop("No *.PRNP.toJunc.bam found in ", JUNC_DIR)

per_junction <- list() # detailed counts per junction, per sample
per_sample   <- list() # aggregated summary per sample

for (bam_junc in bam_files) {
  sample <- sub("\\.PRNP\\.toJunc\\.bam$", "", basename(bam_junc))    # sample name (e.g. CJD1)
  bam_qname <- file.path(JUNC_DIR, paste0(sample, ".PRNP.qname.bam")) # used to count total PRNP fragments
  
  ## denominator: unique fragments in PRNP window (primary, non-dup, non-supp)
  y <- scanBam(bam_qname, param=ScanBamParam(what=c("qname","flag")))[[1]]
  denom <- length(unique(y$qname[
    !bitwAnd(y$flag, 0x100) & !bitwAnd(y$flag, 0x400) & !bitwAnd(y$flag, 0x800)
  ]))
  
  ## junction alignments: filter and compute overhang support
  x <- scanBam(bam_junc, param=ScanBamParam(what=c("qname","flag","rname","pos","cigar","mapq")))[[1]]
  keep <- !bitwAnd(x$flag, 0x100) & !bitwAnd(x$flag, 0x400) & !bitwAnd(x$flag, 0x4) &
    !bitwAnd(x$flag, 0x800) & x$mapq >= MIN_MAPQ & x$rname %in% join_tbl$rname
  
  if (!any(keep)) {
    pj <- join_tbl %>% transmute(sample, rname, n_frags=0L, n_reads=0L, n_starts=0L,
                                 plus_frags=0L, minus_frags=0L, denom_frags=denom,
                                 JPHK=0)
    per_junction[[sample]] <- pj
    ps <- data.frame(
      sample=sample, prnp_total_frags=0L, prnp_total_reads=0L, prnp_total_starts=0L,
      prnp_plus_frags=0L, prnp_minus_frags=0L, denom_frags=denom,
      prnp_JPHK_total=0, n_junctions_pos=0L,
      passes_min_frags=FALSE, passes_strand=FALSE, passes_dispersion=FALSE,
      passes_all=FALSE, stringsAsFactors=FALSE
    )
    per_sample[[sample]] <- ps
    next
  }
  
  df <- data.frame(
    qname = x$qname[keep],
    flag  = x$flag[keep],
    rname = as.character(x$rname[keep]),
    pos   = x$pos[keep],
    cigar = x$cigar[keep],
    mapq  = x$mapq[keep],
    stringsAsFactors = FALSE
  )
  df$end    <- df$pos + cigarWidthAlongReferenceSpace(df$cigar) - 1L
  df$strand <- ifelse(bitwAnd(df$flag, 16) != 0, "-", "+")
  
  df <- df %>% left_join(join_tbl, by="rname")
  df$support <- (df$pos <= (df$join_bp - MIN_OVER + 1L)) &
    (df$end >= (df$join_bp + MIN_OVER))
  
  supp <- df %>% filter(support) %>%
    group_by(rname) %>%
    summarise(
      n_reads     = n(),
      n_frags     = n_distinct(qname),
      n_starts    = n_distinct(pos),
      plus_frags  = n_distinct(qname[strand=="+"]),
      minus_frags = n_distinct(qname[strand=="-"]),
      .groups = "drop"
    )
  
  pj <- join_tbl %>% select(rname) %>% left_join(supp, by="rname") %>%
    mutate(
      n_reads     = coalesce(n_reads,     0L),
      n_frags     = coalesce(n_frags,     0L),
      n_starts    = coalesce(n_starts,    0L),
      plus_frags  = coalesce(plus_frags,  0L),
      minus_frags = coalesce(minus_frags, 0L),
      sample      = sample,
      denom_frags = denom,
      JPHK        = 1e5 * n_frags / pmax(denom, 1L)
    ) %>%
    select(sample, rname, n_frags, n_reads, n_starts,
           plus_frags, minus_frags, denom_frags, JPHK)
  
  per_junction[[sample]] <- pj
  
  # aggregate per-sample
  ps <- pj %>%
    summarise(
      prnp_total_frags  = sum(n_frags),
      prnp_total_reads  = sum(n_reads),
      prnp_total_starts = sum(n_starts),
      prnp_plus_frags   = sum(plus_frags),
      prnp_minus_frags  = sum(minus_frags),
      denom_frags       = first(denom_frags),
      prnp_JPHK_total   = 1e5 * sum(n_frags) / pmax(first(denom_frags), 1L),
      n_junctions_pos   = sum(n_frags > 0L)
    )
  
  ps$sample            <- sample
  ps$passes_min_frags  <- ps$prnp_total_frags >= MIN_FRAGS
  ps$passes_strand     <- (ps$prnp_plus_frags > 0L & ps$prnp_minus_frags > 0L)
  ps$passes_dispersion <- ps$prnp_total_starts >= MIN_STARTS
  ps$passes_all        <- with(ps, passes_min_frags & passes_strand & passes_dispersion)
  
  per_sample[[sample]] <- ps %>% select(sample, everything())
}

per_junction <- bind_rows(per_junction)
per_sample   <- bind_rows(per_sample)

## ---- write outputs ----
OUTDIR <- Sys.getenv("PRNP_JUNCTION_COUNT_DIR", unset = file.path(repo_root, "results", "junctions", "junction_counts"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

outfile1 <- file.path(OUTDIR, "prnp_junction_counts.tsv")
outfile2 <- file.path(OUTDIR, "prnp_junction_summary.tsv")

write.table(per_junction, outfile1, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(per_sample,   outfile2, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote:\n  ", outfile1, "\n  ", outfile2, "\n")
