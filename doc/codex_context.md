# Codex context: PRNP somatic sequencing pipeline (publication reproducibility)

## Purpose of this repository
This repository provides a reproducible, publication-aligned implementation of a targeted deep sequencing workflow to detect somatic variation in the PRNP locus in human brain tissue.

The goal of the repository is reproducibility of the computational workflow described in the manuscript (inputs ? processing steps ? outputs), not the preservation of any particular biological interpretation.

---

## Non-negotiable invariants (publication pipeline)
Codex must treat the following as immutable unless explicitly instructed otherwise:

### Pipeline behaviour
- Do not change the functional behaviour, processing order, parameters, thresholds, or output formats of the publication pipeline.
- Any change that could alter results must be treated as a “new analysis version” rather than a refactor.

### Scientific thresholds / filters (sequencing)
- Somatic variant filtering logic and thresholds are fixed for the publication pipeline.
- The sequencing limit of detection (LoD) is fixed at **VAF = 0.0081 (0.81%)**, derived from spike-in experiments, and must not be modified for the publication pipeline.

### Resources and authoritative inputs
- **All files in `resources/` are immutable** (never modify in-place).
- **All authoritative files in `authoritative_files/` are immutable** (never modify in-place).
  - If updates are required, create a new versioned file and update references explicitly.

### Default stance on ambiguity
If you (Codex) are unsure whether a change will alter results or pipeline functionality in a relevant way, **STOP and ASK** before editing anything beyond trivial documentation.

---

## Publication version identity
- The repository should maintain a clearly identifiable snapshot representing the exact workflow used for the manuscript (e.g. a tag or branch such as `paper-v1`).
- Anything that changes behaviour should either:
  - be rejected for `paper-v1`, or
  - be implemented on a new version line (e.g. `paper-v2` / `analysis-next`) with explicit documentation of differences.

---

## Pipeline overview (sequencing workflow)
The variant detection workflow follows this required processing order:

1. **Tumour-only Mutect2 on non-diseased control samples** (no Panel of Normals).
2. **Verify** that control samples do not contain exonic PRNP mutations (as required before PoN generation).
3. **Build Panel of Normals (PoN)** using verified controls.
4. **Spike-in analysis** using the PoN to establish the empirical sequencing LoD.
5. **CJD sample processing** using the PoN and the established LoD.

Key note: ddPCR is part of the broader project but is intentionally **out of scope here for now**.

---

## Definition of “true somatic variant” (publication pipeline)
Somatic variants are considered only if they satisfy the publication pipeline constraints, including (but not limited to):

- Mutect2 PASS (post-filtering)
- Adequate depth and read support (including strand balance constraints)
- Base and mapping quality thresholds
- Population database frequency exclusion criteria
- VAF significance criteria (including VAF ? 0.5 where relevant)
- **VAF > 0.0081** (sequencing LoD)

Do not alter these criteria without explicit instruction.

---

## Environment contract
- Conda is the baseline environment mechanism.
- `env/environment.yml` defines the intended environment.
- `env/environment.lock.yml` defines the exact resolved package set for publication reproducibility.

Do not introduce undocumented toolchain dependencies.

---

## Directory contract (what lives where)

### Top-level directories
- `authoritative_files/`
  - Canonical, publication-aligned scripts and TSVs for QC and manifest validation.
  - Examples: `validate_manifest.sh`, `manifest.tsv`, `manifest_qc.tsv`, schema files, and `compute_sequencing_metrics.py`.
  - Treat as immutable for the publication pipeline unless a bug is identified and a controlled update is made.

- `bin/`
  - Small utilities used across the repository (e.g. inventory tooling).

- `config/`
  - User-editable configuration and examples.
  - `preprocessing.env` / `preprocessing.env.example` define pipeline environment variables.
  - `preprocessing_samples.tsv` defines sample metadata / input mapping.
  - `config.local.yaml` is local and should not be treated as canonical.

- `data/`
  - Documentation or placeholders for non-public / external datasets.
  - Should not contain large raw files unless explicitly intended.

- `doc/`
  - Static documentation artefacts (e.g. inventory tables, toolchain docs).

- `docs/`
  - Human-facing documentation intended to guide contributors and tools (including this file).

- `env/`
  - Conda environment specifications and lockfiles.

- `fastq/`
  - Raw FASTQ inputs, organised by batches (`CJD_16_samples`, `CJD_8_samples`, `first_CJD_seq`, `sequencing_of_dilutions`).
  - Includes vendor output and run metadata (e.g. stats reports, md5sums, dataset TSVs).

- `mk/`
  - Makefile includes for Conda and toolchain version reporting.

- `resources/`
  - Reference genomes, indexes, known-sites VCFs, intervals, BEDs, GTFs, capture targets, checksums.
  - Immutable: never modify in place.

- `results/`
  - Final or summarised outputs (QC outputs, final BAMs, run logs, run metadata).
  - Contents are derived; may be regenerated.

- `runs/`
  - Intermediate pipeline outputs (per-batch, per-sample directories).
  - Includes trimmed FASTQs, BAM intermediates, recal tables, logs, and RUN_META files.

- `src/`
  - Source code (including legacy material).
  - `src/pipelines/` contains the current maintained pipeline scripts:
    - `preflight_preprocessing.sh`
    - `preprocessing.sh`
    - `mutect2_controls_no_pon.sh`
  - `src/legacy/` contains historical scripts preserved for provenance; do not refactor. These are the scripts which we will edit and adapt in a new for in ´/src´ to make the code publication-ready

- `tests/`
  - Automated checks. Prefer smoke tests and invariant checks that are fast and deterministic.

---

## Expected outputs (high-level)
Examples of “success outputs” include:
- `results/final_bam/*.bam` and `.bai` (finalised BAMs)
- QC run directories under `results/` containing:
  - `logs/` (stderr/stdout logs)
  - `outputs/` (TSVs such as sequencing metrics)
  - `run_meta.txt` (run metadata)
  - in `src/`, our aim is to create a set of code that runs the entire analysis in a reproducible manner, using relative paths
  - the goal is that a reviewer can download the repo from GitHub and reproduce the analysis. This currently is not possible using the `src/legacy/` scripts

---

## How Codex should work in this repo (behavioural rules)

### Default workflow
1. Start by reading and summarising the relevant files (no edits).
2. Propose a short plan (files to change, risks, tests).
3. Make the smallest possible diff.
4. Show `git diff` and run lightweight checks.

### Do not do this without explicit permission
- Changing GATK/Picard/samtools command-line parameters used in the manuscript pipeline
- Changing thresholds or filtering logic
- Reorganising directory layout or renaming outputs used by downstream steps
- Modifying `resources/` or `authoritative_files/` in place

### If a bug is suspected
- Stop and ask.
- Provide:
  - a minimal reproduction (what command, what input, what output)
  - the expected behaviour
  - a proposed minimal fix
  - a test or validation check to prevent regression

---

## Safe / high-value tasks for Codex
Codex is encouraged to help with:
- Adding `--dry-run` modes and consistent CLI help text
- Adding smoke tests (e.g. manifest validation, “make -n” invariants, schema checks)
- Improving logging and “fail-fast” input validation
- Documentation updates that reflect the actual pipeline behaviour
- Refactors that do not change outputs (only when protected by tests)

---

## Quick start (for Codex sessions)
When starting work in this repository:
- Read this file first.
- Identify whether the task touches “publication pipeline invariants”.
- If it might, stop and ask before changing anything.