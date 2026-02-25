You are working inside the prnp-somatic repository.

Execution phase. You may now perform write operations and run workflow commands, but only within the approved scope from the preflight plan/manifest (/add/prnp-somatic/doc/Codex_to_do/2026-02-25_preflight_plan.md)

Objective:
Execute the approved plan to (1) assess and, if possible, run the junction workflow, and (2) complete the repo hygiene/documentation tasks across workflows, including the temporary AAF-filter change and QC regeneration step.

Execution prerequisites:
- Follow the previously approved preflight plan and permission manifest.
- If the current repo state differs materially from the preflight assumptions, pause and report the mismatch before making broad changes.
- Use only repo-root-relative paths in any new scripts and README command examples.

Permissions and safety constraints:
- You may modify files only inside this repository workspace.
- No network access unless it was explicitly justified in the approved manifest and separately authorised.
- Do not install missing tools/packages automatically.
- If required dependencies for junction workflow execution are missing, stop and report exactly what is missing and how to install/check it manually.
- Do not use destructive git commands or recursive/force deletions (e.g., git reset --hard, git clean -fdx, rm -rf).
- Do not push commits.
- Do not move or delete raw data.
- Preserve existing repository structure and style where possible.

Tasks to execute:
1) Check whether the R utilities and other required tools/packages for the junction workflow are installed.
2) If dependencies are present, execute the junction workflow and report success/failure plus output locations.
3) Ensure each workflow has one entrypoint and uses only relative paths:
   - sanity-check existing entrypoint scripts
   - create missing entrypoint scripts where needed
4) Standardise workflow README files to include:
   - Inputs
   - Command
   - Outputs
5) Define expected output files per workflow in the relevant README(s).
6) Create verification assets:
   - a final outputs list
   - a small checksum verification script
7) Ensure ignore rules/documentation reflect:
   - raw large data remain ignored, with exact placement documented
   - generated index/runtime files remain ignored (e.g., BWA sidecars, run outputs)
8) Edit src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R:
   - comment out the step that applies the AAF filter (temporary manual-check need)
   - rerun the script to regenerate its final output
9) Create a single reference document listing:
   - required tools/packages and current versions (ignore TeX)
   - provenance of required reference artefacts
   - placeholders where provenance is unclear

Documentation/output location for this Codex task:
- Save a copy of the execution plan summary, permission manifest used, and final execution report under:
  - doc/Codex_to_do/
- If the directory does not exist, create it.
- Use clear, deterministic filenames (repo-root-relative), for example:
  - doc/Codex_to_do/preflight_plan_approved.md
  - doc/Codex_to_do/permission_manifest_approved.md
  - doc/Codex_to_do/execution_report.md
  (If different names were proposed/approved in preflight, use those.)

Operational rules:
- Work in small, verifiable steps.
- After each major phase, provide a brief progress update (what changed, what was checked, what outputs were generated).
- For command failures, try at most 2 sensible fixes, then stop and report the blocker.
- Do not broaden scope beyond the approved task list without explicitly flagging it first.
- If you encounter ambiguity about what counts as a “workflow”, use the classification from preflight; if still ambiguous, report and proceed conservatively.
- Prefer minimal edits over broad refactors.
- do not alter analytical steps in the scripts, we are only working on making what already exists more reproducible

Definition of done:
- Junction workflow dependency check completed; workflow executed if dependencies are present; outputs reported.
- Entrypoints audited and missing ones created as needed; relative-path usage enforced in entrypoints/README command examples.
- Workflow README files standardised with Inputs / Command / Outputs.
- Expected outputs documented per workflow.
- Final outputs list and checksum verification script created.
- Ignore/documentation updated for raw data placement and generated index/runtime files.
- AAF filter application step in src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R commented out temporarily; script rerun; regenerated final output path reported.
- Reference document created with tools/packages versions + artefact provenance/placeholders.
- Final execution report written under doc/Codex_to_do/.

Final deliverable format (mandatory):
1) Summary of completed tasks
2) Files created/modified (grouped by directory)
3) Commands executed (high-level list)
4) Outputs generated/regenerated (with paths)
5) Missing dependencies or blockers (if any)
6) Suggested manual checks (especially related to the temporarily disabled AAF filter)
7) Notes on any placeholders left for provenance