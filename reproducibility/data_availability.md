# Data Availability

Sequencing FASTQ files for the targeted PRNP sequencing analysis have been
submitted to NCBI SRA.

- BioProject: `PRJNA1484136`
- BioProject URL: <http://www.ncbi.nlm.nih.gov/bioproject/1484136>
- SRA submission ID: `SUB16298068`
- BioSample accessions: `SAMN61220231` through `SAMN61220265`
- BioSample/SRA release: 2027-07-30, or with release of linked data, whichever
  is first

The submitted SRA package contains 35 samples and 70 paired FASTQ files. The
per-sample mapping between project sample IDs, BioSample accessions, SRA
metadata, and FASTQ filenames is recorded in `sra_accessions.tsv`.

SRA experiment and run accessions are not yet assigned in the local records.
When they become available, update `sra_accessions.tsv` rather than replacing
the project sample IDs or FASTQ filenames; those names are used by the analysis
configuration and provenance manifests.

Processed analysis outputs and expected reproducibility artefacts are listed in
`final_outputs_manifest.tsv`, with checksums in `final_outputs.sha256`.
