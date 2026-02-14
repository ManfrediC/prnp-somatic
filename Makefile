SHELL := /bin/bash
.ONESHELL:

# Modular targets
include mk/versions.mk

.SHELLFLAGS := -eu -o pipefail -c

.PHONY: help versions qc_validate qc_metrics clean_qc

# -------------------------------------------------------------------
# Paths inside the repo (adjust only if you move directories)
# -------------------------------------------------------------------
LEGACY_AUTH_DIR ?= src/legacy/ubuntu/authoritative_files
VALIDATE_SCRIPT := $(LEGACY_AUTH_DIR)/validate_manifest.sh
METRICS_SCRIPT  := $(LEGACY_AUTH_DIR)/compute_sequencing_metrics.py

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

qc_validate:
	mkdir -p "$(QC_DIR)"
	echo "== validate_manifest ==" | tee "$(QC_VALIDATE_LOG)"
	bash -lc "cd '$(LEGACY_AUTH_DIR)' && bash './$$(basename "$(VALIDATE_SCRIPT)")'" \
	  2>&1 | tee -a "$(QC_VALIDATE_LOG)"
	@echo "Log: $(QC_VALIDATE_LOG)"

qc_metrics: qc_validate
	mkdir -p "$(QC_DIR)"
	echo "== compute_sequencing_metrics ==" 
	echo "TSV: $(QC_METRICS_TSV)"
	echo "STDERR: $(QC_METRICS_ERR)"
	# stdout is the TSV; stderr is kept as a log (and also shown on screen)
	bash -lc "cd '$(LEGACY_AUTH_DIR)' && python3 './$$(basename "$(METRICS_SCRIPT)")'" \
	  > "$(QC_METRICS_TSV)" \
	  2> >(tee "$(QC_METRICS_ERR)" >&2)
	# sanity check: consistent column counts
	awk -F'\t' 'NR==1{n=NF; next} NF!=n{print "Column mismatch at line " NR ": " NF " vs " n; exit 1}' "$(QC_METRICS_TSV)"
	@echo "OK: wrote $$(($(wc -l < "$(QC_METRICS_TSV)") - 1)) rows (excluding header)"

clean_qc:
	rm -rf "$(QC_DIR)"
	@echo "Removed: $(QC_DIR)"
