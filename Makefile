SHELL := /bin/bash
.ONESHELL:
REQUIRED_CONDA_ENV ?= prnp-somatic

.SHELLFLAGS := -eu -o pipefail -c

.PHONY: help versions toolchain_lock check_conda qc_validate qc_metrics clean_qc print_qc_paths verify_resources preprocessing_preflight preprocessing_dry preprocessing_run

# -------------------------------------------------------------------
# Paths inside the repo (adjust only if you move directories)
# -------------------------------------------------------------------
AUTH_DIR ?= authoritative_files
VALIDATE_SCRIPT := $(AUTH_DIR)/validate_manifest.sh
METRICS_SCRIPT  := $(AUTH_DIR)/compute_sequencing_metrics.py
MANIFEST_TSV    := $(AUTH_DIR)/manifest.tsv
MANIFEST_QC_TSV := $(AUTH_DIR)/manifest_qc.tsv

# -------------------------------------------------------------------
# Outputs (choose a run label without editing the Makefile)
#   make qc_metrics QC_RUN=2026-02-14_test
# -------------------------------------------------------------------
RESULTS_DIR ?= results
QC_RUN      ?= latest
QC_DIR      := $(RESULTS_DIR)/qc/$(QC_RUN)

QC_VALIDATE_LOG := $(QC_DIR)/validate_manifest.log
QC_METRICS_TSV  := $(QC_DIR)/sequencing_metrics_per_sample.tsv
QC_METRICS_ERR  := $(QC_DIR)/compute_sequencing_metrics.stderr.log

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
	@mkdir -p doc
	@$(MAKE) -s versions > doc/tool_versions.lock.txt
	@echo "Wrote: doc/tool_versions.lock.txt"

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

help:
	@echo "Targets (run one step at a time):"
	@echo "  make versions                  Show key tool versions (fast)"
	@echo "  make qc_validate                Run validate_manifest.sh"
	@echo "  make qc_metrics                 Run compute_sequencing_metrics.py -> TSV"
	@echo "  make clean_qc                   Remove results/qc/<QC_RUN>/"
	@echo ""
	@echo "Common options:"
	@echo "  QC_RUN=latest                   Output subfolder label (default: latest)"
	@echo "  RESULTS_DIR=results              Base results directory"
	@echo ""
	@echo "Examples:"
	@echo "  make versions"
	@echo "  make qc_metrics QC_RUN=2026-02-14_test"

qc_validate: check_conda
	mkdir -p "$(QC_DIR)"
	echo "== validate_manifest ==" | tee "$(QC_VALIDATE_LOG)"
	bash "$(VALIDATE_SCRIPT)" "$(MANIFEST_TSV)" "$(MANIFEST_QC_TSV)" 2>&1 | tee -a "$(QC_VALIDATE_LOG)"
	@echo "Log: $(QC_VALIDATE_LOG)"

qc_metrics: qc_validate
	mkdir -p "$(QC_DIR)"
	echo "== compute_sequencing_metrics =="
	echo "TSV: $(QC_METRICS_TSV)"
	echo "STDERR: $(QC_METRICS_ERR)"
	python "$(METRICS_SCRIPT)" \
	  > "$(QC_METRICS_TSV)" \
	  2> >(tee "$(QC_METRICS_ERR)" >&2)
	test -s "$(QC_METRICS_TSV)"
	awk -F'\t' 'NR==1{n=NF; next} NF!=n{print "Column mismatch at line " NR ": " NF " vs " n; exit 1}' "$(QC_METRICS_TSV)"
	@echo "OK: wrote $$(( $$(wc -l < "$(QC_METRICS_TSV)") - 1 )) rows (excluding header)"

clean_qc:
	rm -rf "$(QC_DIR)"
	@echo "Removed: $(QC_DIR)"

print_qc_paths:
	@echo "AUTH_DIR=$(AUTH_DIR)"
	@echo "VALIDATE_SCRIPT=$(VALIDATE_SCRIPT)"
	@echo "METRICS_SCRIPT=$(METRICS_SCRIPT)"
	@echo "MANIFEST_TSV=$(MANIFEST_TSV)"
	@echo "MANIFEST_QC_TSV=$(MANIFEST_QC_TSV)"
 
verify_resources:
	cd resources && sha256sum -c SHA256SUMS.txt

preprocessing_preflight: check_conda
	@src/pipelines/preflight_preprocessing.sh

preprocessing_dry: check_conda
	@DRY_RUN=1 src/pipelines/preprocessing.sh

preprocessing_run: check_conda
	@DRY_RUN=0 src/pipelines/preprocessing.sh
