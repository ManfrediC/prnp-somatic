SHELL := /bin/bash
.ONESHELL:
REQUIRED_CONDA_ENV ?= prnp-somatic
CONDA_BIN ?= conda
PYTHON := $(shell if command -v python3 >/dev/null 2>&1; then echo python3; elif command -v python >/dev/null 2>&1; then echo python; fi)
R_SCRIPT ?= Rscript

ifeq ($(strip $(PYTHON)),)
  $(error Neither python3 nor python is available on PATH. Install a Python interpreter and retry.)
endif

.SHELLFLAGS := -eu -o pipefail -c

.PHONY: help versions toolchain_lock check_conda \
	ddpcr snv repeats junctions all check \
	dna_quality \
	qc_metrics clean_qc print_qc_paths verify_resources preprocessing_preflight preprocessing_dry preprocessing_run

# -------------------------------------------------------------------
# Paths inside the repo (adjust only if you move directories)
# -------------------------------------------------------------------
METRICS_SCRIPT  := src/sequencing_qc/compute_sequencing_metrics.py
DNA_QUALITY_SCRIPT := src/dna_quality/build_dna_quality_evidence_table.R
SAMPLE_MANIFEST_TSV := config/sample_manifest.tsv
SEQUENCING_METRICS_SCHEMA_TSV := src/sequencing_qc/sequencing_metrics_per_sample.schema.tsv

# -------------------------------------------------------------------
# Outputs
# -------------------------------------------------------------------
RESULTS_DIR       ?= results
SEQUENCING_QC_DIR ?= $(RESULTS_DIR)/sequencing_qc
QC_DIR            := $(SEQUENCING_QC_DIR)
DNA_QUALITY_OUT   := $(RESULTS_DIR)/dna_quality/sample_quality_evidence_table.tsv

QC_METRICS_TSV  := $(QC_DIR)/sequencing_metrics_per_sample.tsv
QC_METRICS_ERR  := $(QC_DIR)/compute_sequencing_metrics.stderr.log

# -------------------------------------------------------------------
# Toolchain reporting helpers
# -------------------------------------------------------------------
versions:
	@set -eu
	@echo "== Tool versions =="
	@echo "captured_utc: $$(date -u +%FT%TZ)"
	@echo "git_commit: $$(git rev-parse --short HEAD 2>/dev/null || echo NA)"
	@echo ""
	@echo "-- OS --"
	@uname -a || true
	@echo ""
	@echo "-- Core --"
	@command -v bash >/dev/null 2>&1 && bash --version 2>/dev/null | sed -n '1p' || echo "bash: NOT FOUND"
	@command -v make >/dev/null 2>&1 && make --version 2>/dev/null | sed -n '1p' || echo "make: NOT FOUND"
	@command -v rsync >/dev/null 2>&1 && rsync --version 2>/dev/null | sed -n '1p' || echo "rsync: NOT FOUND"
	@command -v grep >/dev/null 2>&1 && grep --version 2>/dev/null | sed -n '1p' || echo "grep: NOT FOUND"
	@command -v sed  >/dev/null 2>&1 && sed  --version 2>/dev/null | sed -n '1p' || echo "sed: NOT FOUND"
	@command -v awk  >/dev/null 2>&1 && awk  --version 2>/dev/null | sed -n '1p' || echo "awk: NOT FOUND"
	@echo ""
	@echo "-- HTS / alignment --"
	@command -v bwa >/dev/null 2>&1 && { bwa 2>&1 | awk -F': ' '/^Version:/{print "bwa " $$2; found=1} END{if(!found) print "bwa (version not detected)"}' || true; } || echo "bwa: NOT FOUND"
	@command -v samtools >/dev/null 2>&1 && samtools --version 2>/dev/null | sed -n '1,2p' || echo "samtools: NOT FOUND"
	@command -v bcftools >/dev/null 2>&1 && bcftools --version 2>/dev/null | sed -n '1,2p' || echo "bcftools: NOT FOUND"
	@command -v bgzip >/dev/null 2>&1 && bgzip --version 2>/dev/null | sed -n '1p' || echo "bgzip: NOT FOUND"
	@command -v tabix >/dev/null 2>&1 && tabix --version 2>/dev/null | sed -n '1p' || echo "tabix: NOT FOUND"
	@command -v bedtools >/dev/null 2>&1 && bedtools --version 2>/dev/null | sed -n '1p' || echo "bedtools: NOT FOUND"
	@echo ""
	@echo "-- Java / GATK --"
	@command -v java >/dev/null 2>&1 && java -version 2>&1 | sed -n '1,3p' || echo "java: NOT FOUND"
	@command -v gatk >/dev/null 2>&1 && gatk --version 2>/dev/null || echo "gatk: NOT FOUND"
	@echo ""
	@echo "-- Python --"
	@command -v python3 >/dev/null 2>&1 && python3 --version 2>/dev/null || echo "python3: NOT FOUND"
	@command -v pip3 >/dev/null 2>&1 && pip3 --version 2>/dev/null || echo "pip3: NOT FOUND"
	@echo ""
	@echo "-- QC helpers (optional) --"
	@command -v fastqc >/dev/null 2>&1 && fastqc --version 2>/dev/null || echo "fastqc: NOT FOUND"
	@command -v multiqc >/dev/null 2>&1 && multiqc --version 2>/dev/null || echo "multiqc: NOT FOUND"
	@echo ""
	@echo "-- R (optional) --"
	@command -v R >/dev/null 2>&1 && R --version 2>/dev/null | sed -n '1,2p' || echo "R: NOT FOUND"

toolchain_lock:
	@set -eu
	@mkdir -p reproducibility
	# Persist a single, greppable snapshot of the current toolchain.
	@$(MAKE) -s versions > reproducibility/tool_versions.lock.txt
	@echo "Wrote: reproducibility/tool_versions.lock.txt"

# -------------------------------------------------------------------
# Conda guardrail (non-base env required; optional exact name check)
# -------------------------------------------------------------------
check_conda:
	@set -eu
	@if [ -z "$$CONDA_PREFIX" ] || [ "$$CONDA_DEFAULT_ENV" = "base" ]; then \
	  echo "ERROR: a non-base Conda environment is required."; \
	  echo "Create one from env specs in env/ (see env/README.md)."; \
	  echo "Then activate it before running this target."; \
	  exit 1; \
	fi
	@if [ -n "$${REQUIRED_CONDA_ENV:-}" ] && [ "$$CONDA_DEFAULT_ENV" != "$$REQUIRED_CONDA_ENV" ]; then \
	  echo "ERROR: active Conda environment is '$$CONDA_DEFAULT_ENV', expected '$$REQUIRED_CONDA_ENV'."; \
	  echo "Run: conda activate $$REQUIRED_CONDA_ENV"; \
	  exit 1; \
		fi

# -------------------------------------------------------------------
# Main workflow entrypoints
# -------------------------------------------------------------------
help:
	@echo "Targets (run one step at a time):"
	@echo "  make ddpcr                     Run ddPCR workflow (requires env: prnp-somatic-ddpcr)"
	@echo "  make snv                       Run SNV Stage-12 wrapper (requires env: prnp-somatic)"
	@echo "  make repeats                   Run PRNP ORR repeat workflow (requires env: prnp-repeats)"
	@echo "  make repeats_manual_controls   Run manual PRNP ORR mosaic review on controls"
	@echo "  make repeats_manual_cjd        Run manual PRNP ORR mosaic review on CJD samples"
	@echo "  make repeats_manual_filter_cjd Filter CJD manual review summary against controls"
	@echo "  make junctions                 Run junction workflow (requires env: prnp-junctions)"
	@echo "  make dna_quality               Build DNA-quality evidence table"
	@echo "  make all                       Run ddpcr + snv + junctions via conda run"
	@echo "  make check                     Run resource/output integrity checks"
	@echo "  make versions                  Show key tool versions (fast)"
	@echo "  make qc_metrics                 Run compute_sequencing_metrics.py -> TSV"
	@echo "  make clean_qc                   Remove results/sequencing_qc/"
	@echo ""
	@echo "Common options:"
	@echo "  RESULTS_DIR=results              Base results directory"
	@echo "  SEQUENCING_QC_DIR=results/sequencing_qc"
	@echo ""
	@echo "Examples:"
	@echo "  make versions"
	@echo "  make ddpcr"
	@echo "  make junctions"
	@echo "  make dna_quality"
	@echo "  make qc_metrics"

ddpcr: REQUIRED_CONDA_ENV=prnp-somatic-ddpcr
ddpcr: check_conda
	@bash src/ddpcr/run_ddpcr.sh

snv: REQUIRED_CONDA_ENV=prnp-somatic
snv: check_conda
	@bash src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh

repeats: REQUIRED_CONDA_ENV=prnp-repeats
repeats: check_conda
	@bash src/repeats/01_run_prnp_orr.sh

repeats_manual_controls: REQUIRED_CONDA_ENV=prnp-repeats
repeats_manual_controls: check_conda
	@$(PYTHON) src/repeats/07_run_manual_mosaic_prnp_orr_cohort.py --cohort controls

repeats_manual_cjd: REQUIRED_CONDA_ENV=prnp-repeats
repeats_manual_cjd: check_conda
	@$(PYTHON) src/repeats/07_run_manual_mosaic_prnp_orr_cohort.py --cohort cjd

repeats_manual_filter_cjd: REQUIRED_CONDA_ENV=prnp-repeats
repeats_manual_filter_cjd: check_conda
	@$(PYTHON) src/repeats/08_filter_manual_mosaic_prnp_orr_cohort.py \
		--input-summary results/repeats/manual_cohort/cjd/cohort_summary.tsv \
		--background-summary results/repeats/manual_cohort/controls/cohort_summary.tsv \
		--require-background-exceedance \
		--output-prefix results/repeats/manual_cohort/cjd/filtered/default

junctions: REQUIRED_CONDA_ENV=prnp-junctions
junctions: check_conda
	@TMPDIR=/tmp TEMP=/tmp TMP=/tmp bash src/junctions/run_junctions.sh

all:
	# Run each workflow in its expected environment without requiring manual activation.
	@command -v "$(CONDA_BIN)" >/dev/null 2>&1 || { echo "ERROR: conda not found in PATH."; exit 1; }
	@echo "== [1/3] ddPCR (env: prnp-somatic-ddpcr) =="
	@"$(CONDA_BIN)" run -n prnp-somatic-ddpcr bash src/ddpcr/run_ddpcr.sh
	@echo "== [2/3] SNV Stage-12 wrapper (env: prnp-somatic) =="
	@"$(CONDA_BIN)" run -n prnp-somatic bash src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh
	@echo "== [3/3] Junctions (env: prnp-junctions) =="
	@"$(CONDA_BIN)" run -n prnp-junctions bash -lc 'TMPDIR=/tmp TEMP=/tmp TMP=/tmp bash src/junctions/run_junctions.sh'

check: REQUIRED_CONDA_ENV=prnp-somatic
check: check_conda verify_resources
	# Validate tracked final-output checksums after resource integrity checks.
	@bash reproducibility/verify_output_checksums.sh --mode check

# -------------------------------------------------------------------
# QC helpers
# -------------------------------------------------------------------
qc_metrics: check_conda
	mkdir -p "$(QC_DIR)"
	echo "== compute_sequencing_metrics =="
	echo "TSV: $(QC_METRICS_TSV)"
	echo "STDERR: $(QC_METRICS_ERR)"
	PRNP_MANIFEST="$(SAMPLE_MANIFEST_TSV)" \
	PRNP_SCHEMA="$(SEQUENCING_METRICS_SCHEMA_TSV)" \
	$(PYTHON) "$(METRICS_SCRIPT)" \
	  > "$(QC_METRICS_TSV)" \
	  2> >(tee "$(QC_METRICS_ERR)" >&2)
	test -s "$(QC_METRICS_TSV)"
	awk -F'\t' 'NR==1{n=NF; next} NF!=n{print "Column mismatch at line " NR ": " NF " vs " n; exit 1}' "$(QC_METRICS_TSV)"
	@echo "OK: wrote $$(( $$(wc -l < "$(QC_METRICS_TSV)") - 1 )) rows (excluding header)"

dna_quality:
	mkdir -p "$$(dirname "$(DNA_QUALITY_OUT)")"
	echo "== build DNA-quality evidence table =="
	echo "TSV: $(DNA_QUALITY_OUT)"
	$(R_SCRIPT) "$(DNA_QUALITY_SCRIPT)" \
	  --output "$(DNA_QUALITY_OUT)"
	test -s "$(DNA_QUALITY_OUT)"
	@echo "OK: wrote $(DNA_QUALITY_OUT)"

clean_qc:
	# Remove the canonical sequencing QC output folder.
	rm -rf "$(QC_DIR)"
	@echo "Removed: $(QC_DIR)"

print_qc_paths:
	@echo "METRICS_SCRIPT=$(METRICS_SCRIPT)"
	@echo "SAMPLE_MANIFEST_TSV=$(SAMPLE_MANIFEST_TSV)"
	@echo "SEQUENCING_METRICS_SCHEMA_TSV=$(SEQUENCING_METRICS_SCHEMA_TSV)"
	@echo "QC_DIR=$(QC_DIR)"
	@echo "QC_METRICS_TSV=$(QC_METRICS_TSV)"
	@echo "QC_METRICS_ERR=$(QC_METRICS_ERR)"
 
verify_resources:
	cd resources && sha256sum -c SHA256SUMS.txt

# -------------------------------------------------------------------
# Preprocessing wrappers (same scripts used in src/pipelines)
# -------------------------------------------------------------------
preprocessing_preflight: check_conda
	@src/pipelines/preflight_preprocessing.sh

preprocessing_dry: check_conda
	@DRY_RUN=1 src/pipelines/preprocessing.sh

preprocessing_run: check_conda
	@DRY_RUN=0 src/pipelines/preprocessing.sh
