# Toolchain (recommended baseline)

This repository was developed and validated using the **Conda** environment defined in `env/environment.yml`.
The exact resolved package set used for analyses is recorded in `env/environment.lock.yml`.

To inspect versions on your machine (after activating Conda):
- `conda activate prnp-somatic`
- `make versions`

To write a snapshot into the repository:
- `make toolchain_lock` (writes `doc/tool_versions.lock.txt`)

## Conda-managed tools (baseline for reproducibility)
The versions below are the expected baseline when using `env/environment.yml`
(and are fully resolved/pinned in `env/environment.lock.yml`).

- bwa: 0.7.17-r1188
- samtools: 1.20 (htslib 1.20)
- bcftools: 1.20 (htslib 1.20)
- bgzip/tabix: 1.20 (htslib 1.20)
- bedtools: 2.31.1
- FastQC: 0.12.1
- MultiQC: 1.33
- OpenJDK: 17.0.18
- GATK: 4.5.0.0 (HTSJDK 4.1.0; Picard 3.1.1)
- Python: 3.10.14
- pip: 26.0.1
- R: 4.4.1

## System tools (not pinned; varies by OS)
These are commonly used utilities. Exact versions may differ across machines/CI images.

- bash: 5.1.16
- GNU Make: 4.3
- rsync: 3.2.7
- grep: 3.7
- sed: 4.8
- gawk: 5.1.0
