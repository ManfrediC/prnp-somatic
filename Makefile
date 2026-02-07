SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c

.PHONY: help qc_metrics

RUN_ID ?= $(shell date -u +%F_%H%M%S)
RUN_DIR ?= results/$(RUN_ID)_qc_metrics
LOG_DIR := $(RUN_DIR)/logs
OUT_DIR := $(RUN_DIR)/outputs

LEGACY_AUTH_DIR := src/legacy/ubuntu/authoritative_files
VALIDATE_SCRIPT := $(LEGACY_AUTH_DIR)/validate_manifest.sh
METRICS_SCRIPT  := $(LEGACY_AUTH_DIR)/compute_sequencing_metrics.py

METRICS_TSV := $(OUT_DIR)/sequencing_metrics_per_sample.tsv
METRICS_STDERR_LOG := $(LOG_DIR)/compute_sequencing_metrics.stderr.log

help:
	@echo "Targets:"
	@echo "  make qc_metrics            Run manifest validation + sequencing metrics"
	@echo ""
	@echo "Options:"
	@echo "  RUN_ID=...                 Override run id (default: UTC timestamp)"
	@echo "  RUN_DIR=...                Override output directory"

qc_metrics:
	mkdir -p "$(LOG_DIR)" "$(OUT_DIR)"

	{
	  echo "run_id: $(RUN_ID)"
	  echo "run_dir: $(RUN_DIR)"
	  echo "date_utc: $$(date -u +%FT%TZ)"
	  echo "git_commit: $$(git rev-parse HEAD 2>/dev/null || echo NA)"
	} | tee "$(RUN_DIR)/run_meta.txt"

	echo ""
	echo "== Step 1: validate manifest =="
	bash -lc "cd '$(LEGACY_AUTH_DIR)' && bash './validate_manifest.sh'" \
	  2>&1 | tee "$(LOG_DIR)/validate_manifest.log"

	echo ""
	echo "== Step 2: compute sequencing metrics =="
	echo "Writing TSV to: $(METRICS_TSV)"
	# stdout -> TSV file; stderr -> console + stderr log
	bash -lc "cd '$(LEGACY_AUTH_DIR)' && PIPELINE_RUN_ID='$(RUN_ID)' python3 './compute_sequencing_metrics.py'" \
	  > "$(METRICS_TSV)" \
	  2> >(tee "$(METRICS_STDERR_LOG)" >&2)

	# Quick integrity check: consistent column count across rows
	awk -F'\t' 'NR==1{n=NF; next} NF!=n{print "Column mismatch at line " NR ": " NF " vs " n; exit 1}' "$(METRICS_TSV)"
	echo "OK: TSV written with consistent columns."
	echo "Rows: $$(($(wc -l < "$(METRICS_TSV)") - 1)) (excluding header)"

	echo ""
	echo "Done."
	echo "Outputs: $(OUT_DIR)"
	echo "Logs:    $(LOG_DIR)"
