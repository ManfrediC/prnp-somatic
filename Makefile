# Use bash and make each recipe fail loudly on errors.
SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c

.PHONY: help qc_metrics

# ---- Run metadata ----
# Default: UTC timestamp so runs are unambiguous across machines/timezones.
RUN_ID ?= $(shell date -u +%F_%H%M%S)
RUN_DIR ?= results/$(RUN_ID)_qc_metrics
LOG_DIR := $(RUN_DIR)/logs
OUT_DIR := $(RUN_DIR)/outputs

# ---- Current "authoritative" scripts (as-is) ----
LEGACY_AUTH_DIR := src/legacy/ubuntu/authoritative_files
VALIDATE_SCRIPT := $(LEGACY_AUTH_DIR)/validate_manifest.sh
METRICS_SCRIPT  := $(LEGACY_AUTH_DIR)/compute_sequencing_metrics.py

help:
	@echo "Targets:"
	@echo "  make qc_metrics            Run manifest validation + sequencing metrics"
	@echo ""
	@echo "Options:"
	@echo "  RUN_ID=...                 Override run id (default: UTC timestamp)"
	@echo "  RUN_DIR=...                Override output directory"
	@echo ""
	@echo "Example:"
	@echo "  make qc_metrics RUN_ID=2026-02-07_qc"

qc_metrics:
	mkdir -p "$(LOG_DIR)" "$(OUT_DIR)"

	# Record minimal provenance (useful for reviewers and for you).
	{
	  echo "run_id: $(RUN_ID)"
	  echo "run_dir: $(RUN_DIR)"
	  echo "date_utc: $$(date -u +%FT%TZ)"
	  echo "git_commit: $$(git rev-parse HEAD 2>/dev/null || echo NA)"
	} | tee "$(RUN_DIR)/run_meta.txt"

	echo ""
	echo "== Step 1: validate manifest =="
	# Run from the script's directory to avoid relative-path surprises.
	bash -lc "cd '$(LEGACY_AUTH_DIR)' && bash './validate_manifest.sh'" \
	  2>&1 | tee "$(LOG_DIR)/validate_manifest.log"

	echo ""
	echo "== Step 2: compute sequencing metrics =="
	bash -lc "cd '$(LEGACY_AUTH_DIR)' && python3 './compute_sequencing_metrics.py'" \
	  2>&1 | tee "$(LOG_DIR)/compute_sequencing_metrics.log"

	echo ""
	echo "Done. Logs in: $(LOG_DIR)"
	echo "Outputs folder (for copied artefacts later): $(OUT_DIR)"
