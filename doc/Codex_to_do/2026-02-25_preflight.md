You are working inside the prnp-somatic repository.

Preflight only. Do NOT execute any write operations or workflow commands yet.

Objective:
Prepare to check whether the R utilities and other required tools for the junction workflow are installed, then execute the junction workflow, while also performing a repo hygiene/documentation pass across workflows (see below):

Task scope (ignore original numbering and treat this as a unified task list):
- Check whether the R utilities and other required tools for the junction part of the workflow have been installed, then execute the junction workflow.
- Ensure every workflow has one entrypoint (e.g., run_ddpcr.sh, run_junctions.sh) and uses only relative paths. Create missing entrypoint scripts and check scripts already created.
- Define expected output files per workflow in README.
- Create a list of final outputs and a small checksums script for verification.
- Keep raw large data ignored, but document exactly where it should be placed.
- Keep generated index/runtime files ignored (e.g., BWA sidecars, run outputs).
- Standardise each workflow README with: inputs, command, outputs.
- In /add/prnp-somatic/src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R, comment out the step that applies the AAF filter (temporary manual-check need), then run the script to regenerate its final output.
- Create a list of packages and tools with current versions required to run everything (ignore TeX for now). Record provenance of required reference artefacts in the same reference document. If provenance is unclear, leave placeholders for manual completion by me.

Important:
- Use only relative paths in any new scripts or README commands (repo-root-relative execution assumptions are fine).
- Do not install anything yet. First assess what is missing.
- Do not use destructive git or rm commands.
- Do not push commits.
- Do not move or delete raw data.

Preflight deliverables (mandatory, in this order, saved in /add/prnp-somatic/doc/Codex_to_do):
1) A concise step-by-step execution plan (grouped into phases)
2) A permission manifest including:
   - Files and directories you expect to read
   - Files and directories you expect to create/modify
   - Shell commands you expect to run
   - Whether network access is needed (and exactly why)
   - Any potentially risky commands and safer alternatives
3) A dependency-check plan for the junction workflow:
   - Which tools/packages you will check
   - How you will detect installed vs missing
   - What you will do if something is missing (stop and report, unless explicitly authorised later)
4) A run plan for the junction workflow and the QC regeneration step, including expected output locations
5) A short list of assumptions/unknowns that could block execution

Stop after the preflight manifest and wait for approval.