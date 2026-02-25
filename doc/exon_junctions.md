# Exon-exon junctions

## Construction of a PRNP exon–exon junction reference

To investigate whether PRNP cDNA insertions or retrocopies were present in the genomes of CJD and control brains, suggested by a report of APP exon-exon junctions in Alzheimer’s disease [44], we searched for DNA fragments spanning exon–exon junctions within our targeted PRNP sequencing libraries.

We first constructed a set of synthetic PRNP exon–exon junction reference sequences. Exon coordinates were extracted for three two-exon spliced PRNP transcripts of interest (ENST00000379440, ENST00000430350 and ENST00000457586). A synthetic junction sequence was generated for every exon pair (exon 1 → exon 2),  defined as up to 75 bp from the 3′ end of the upstream exon and up to 75 bp from the 5′ end of the downstream exon (truncated if the exon was shorter than 75 bp), by extracting sequences from the GRChg38 reference FASTA (GATK resource bundle) and end to end concatenation. This yielded, for each transcript, a single PRNP exon1–exon2 junction contig with a known junction coordinate at the boundary between the two exon derived segments. The resulting junction sequences were written to FASTA and indexed for alignment with BWA MEM.

## Extraction of PRNP region fragments from targeted sequencing BAMs

To isolate fragments likely to contain any PRNP derived cDNA insertions, we defined a minimal PRNP locus window consisting of PRNP padded with 1kb upstream and downstream. All fragments whose alignments overlapped the padded PRNP region were extracted with samtools from the final duplicate marked, base quality recalibrated BAM file used in the somatic variant analysis. Reads marked as duplicates were excluded and BAMs were converted to paired FASTQ files with samtools, retaining read pairs where both mates had passed the PRNP window and duplicate filters. Singletons and orphan reads were discarded from this analysis.

## Alignment to the PRNP junction reference

For each sample, the PRNP window read pairs were realigned to the synthetic PRNP exon-exon junction reference FASTA using BWA MEM with default settings. The resulting alignments were coordinate sorted and indexed with samtools.

## Quantification of junction spanning fragments and diagnostic filtering

To quantify evidence for genomic PRNP exon–exon junctions, we analysed the junction aligned BAMs in R using Rsamtools and GenomicAlignments. For each junction contig we knew the precise junction coordinate (the length of the upstream exon derived segment). For each sample, we read all alignments to the junction contigs and applied the following filters: (i) Excluded unmapped (0x4 flag), secondary (0x100), supplementary (0x800) and duplicate (0x400) alignments; (ii) retained alignments with mapping quality ≥ 20; (iii) retained reads whose aligned span crossed the exon-exon junction with a minimum overhang of 10 bp upstream and downstream of the junction coordinate.
