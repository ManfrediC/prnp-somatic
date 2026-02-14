# Toolchain (recommended baseline)

This repository was developed and validated with the tool versions listed below.
Where possible, use the same **major.minor** versions to minimise behavioural differences.

To see the versions on your machine:
- `make versions`

To write a snapshot into the repository:
- `make toolchain_lock` (writes `doc/tool_versions.lock.txt`)

## Core system tools
- bash: 5.1.16
- GNU Make: 4.3
- rsync: 3.2.7
- grep: 3.7
- sed: 4.8
- gawk: 5.1.0

## Alignment / HTS tools
- bwa: 0.7.17-r1188
- samtools: 1.13 (htslib 1.13+ds)
- bcftools: 1.13 (htslib 1.13+ds)
- bgzip/tabix: 1.19.1
- bedtools: 2.30.0

## Java / GATK
- OpenJDK: 17.0.18
- GATK: 4.5.0.0 (HTSJDK 4.1.0; Picard 3.1.1)

## Python
- Python: 3.10.12
- pip: 22.0.2

## QC helpers (optional)
- FastQC: 0.11.9
- MultiQC: 1.27

## R (optional)
- R: 4.4.1
