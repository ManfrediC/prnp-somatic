# Data Availability

Sequencing FASTQ files for the targeted PRNP sequencing analysis have been
submitted to NCBI SRA and processed.

- BioProject: `PRJNA1484136`
- BioProject URL: <https://www.ncbi.nlm.nih.gov/bioproject/1484136>
- SRA submission ID: `SUB16298068`
- SRA study accession: `SRP714079`
- BioSample accessions: `SAMN61220231` through `SAMN61220265`
- BioSample/SRA release: 2027-07-30 or upon publication, whichever is first

The submitted SRA package contains 35 samples and 70 paired FASTQ files. The
per-sample mapping between project sample IDs, BioSample accessions, SRA
metadata, and FASTQ filenames is recorded in `sra_accessions.tsv`.

All 35 SRA run accessions are recorded in `sra_accessions.tsv`. SRX experiment
accessions were not included in the processed portal export available locally;
the table records that limitation explicitly rather than implying that the SRR
accessions are still pending.

Raw and project-derived ddPCR data have been uploaded to an unpublished Zenodo
draft with reserved DOI `10.5281/zenodo.21227118`. The draft contains 18 files,
is embargoed until 2027-07-30, and should be released earlier if the preprint or
manuscript is published first. The DOI becomes registered when the Zenodo
record is published.

Processed analysis outputs and expected reproducibility artefacts are listed in
`final_outputs_manifest.tsv`, with checksums in `final_outputs.sha256`.
