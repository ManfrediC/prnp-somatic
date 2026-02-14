.PHONY: check_conda

check_conda:
	@set -eu
	@if [ -z "$$CONDA_PREFIX" ] || [ "$$CONDA_DEFAULT_ENV" = "base" ]; then \
	  echo "ERROR: Conda environment 'prnp-somatic' is not active."; \
	  echo "Run: conda env create -f env/environment.yml"; \
	  echo "Then: conda activate prnp-somatic"; \
	  exit 1; \
	fi
