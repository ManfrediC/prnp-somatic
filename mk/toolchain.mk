.PHONY: toolchain_lock

toolchain_lock:
	@set -eu
	@mkdir -p doc
	@$(MAKE) -s versions > doc/tool_versions.lock.txt
	@echo "Wrote: doc/tool_versions.lock.txt"
