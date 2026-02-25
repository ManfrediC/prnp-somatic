# mk modules

This directory contains modular `make` fragments included by the repository root `Makefile`.

## Files

- `versions.mk`
  - Defines target: `versions`
  - Prints versions for OS/core tools, HTS tools, Java/GATK, Python, optional QC helpers, and R.

- `toolchain.mk`
  - Defines target: `toolchain_lock`
  - Writes a reproducible snapshot of `make versions` output to `doc/tool_versions.lock.txt`.

- `conda.mk`
  - Defines target: `check_conda`
  - Fails early unless a non-`base` conda environment is active.
  - Used by root Makefile targets that require a configured analysis environment.

## How these are used

The root `Makefile` includes these modules:

```make
include mk/versions.mk
include mk/toolchain.mk
include mk/conda.mk
```

That keeps shared target logic in one place while the root `Makefile` focuses on workflow/QC orchestration.
