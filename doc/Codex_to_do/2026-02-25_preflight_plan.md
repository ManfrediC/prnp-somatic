# 2026-02-25 Preflight Plan

## 1) Step-by-step execution plan (grouped into phases)

### Phase A: Baseline and scope lock
1. Confirm current repo state and target workflows:
   - `src/ddPCR`
   - `src/junctions`
   - `src/pipelines`
2. Confirm the exact files to touch for this task:
   - workflow entrypoints and README standardization
   - checksum/final-output listing utilities
   - temporary AAF-filter bypass in `src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R`
3. Confirm no destructive operations and no moving/deleting raw data.

### Workflow classification (explicit scope)
In scope workflow families:
- `ddPCR` reproducibility packaging (`src/ddPCR`)
- `junctions` dependency check + run + output verification (`src/junctions`)
- `pipelines` Stage-12 CJD+dilutions QC regeneration with temporary AAF bypass (`src/pipelines/12_*`)
- cross-workflow reproducibility docs/checksum utilities

Out of scope workflow families:
- full preprocessing and full Mutect2 end-to-end reruns (`src/pipelines/preprocessing.sh`, stages 1-11 reruns)
- manuscript rendering/build (R/LaTeX compile for figures/tables)
- legacy script migration outside the three families above (`src/legacy/**`, `prnp-junctions/**` edits)
- TeX toolchain/environment setup

### Phase B: Dependency and tool preflight (no install)
1. Check command-line tools required for junction workflow.
2. Check required R packages for junction scripts.
3. Check required input files/dirs exist for junction run and Stage-12 QC regeneration.
4. If anything is missing, stop and report only (do not install unless explicitly authorized later).

### Phase C: Reproducibility/hygiene edits (minimal changes)
1. Ensure one clear entrypoint per workflow family:
   - ddPCR: add/confirm `src/ddPCR/run_ddpcr.sh`
   - junctions: confirm `src/junctions/run_junctions.sh`
   - pipelines: add/confirm `src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh` as a single publication-path wrapper for Stage-12 QC regeneration
2. Ensure entrypoints use repo-root-relative paths only.
3. Standardize README sections per workflow: `Inputs`, `Command`, `Outputs`.
4. Define expected final outputs per workflow in README(s) and/or a central reproducibility doc.
5. Add a lightweight checksum verifier and final-output manifest script (read-only verification).
6. Confirm `.gitignore` coverage for raw large data and generated runtime/index files.

### Phase D: Requested analysis edit + execution
1. In `src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R`, temporarily bypass/comment out the AAF-filter application block (minimal, clearly marked temporary change).
2. Run the junction workflow via the single entrypoint.
3. Run the Stage-12 QC regeneration command for CJD+dilutions.
4. Capture output inventory and checksums.

### Phase E: Verification and closeout
1. Verify produced files against expected outputs list.
2. Compare key regenerated outputs to legacy counterparts where available.
3. Record tool/package versions and reference artifact provenance into a single reference document (with placeholders where provenance is unknown).
4. Summarize residual gaps and approval gates before final publication-ready pass.

### Exact proposed new files (repo-root-relative)
1. `src/ddPCR/run_ddpcr.sh`
2. `src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh`
3. `bin/verify_output_checksums.sh`
4. `doc/reproducibility/final_outputs_manifest.tsv`
5. `doc/reproducibility/tooling_and_reference_provenance.md`

### Estimated edit scope
- Estimated files touched: ~11 to ~15 files total
- Estimated new files: 5 files (listed above)
- Estimated existing files modified: ~6 to ~10 files
- Estimated directories touched: ~6 directories
  - `src/ddPCR`
  - `src/junctions`
  - `src/pipelines`
  - `doc/reproducibility`
  - `bin`
  - repo root (for `.gitignore` if needed)

## 2) Permission manifest

### Files/directories expected to read
- `README.md`
- `src/ddPCR/*`
- `src/junctions/*`
- `src/pipelines/*`
- `config/junctions.env*`
- `config/preprocessing.env*`
- `ddPCR/*`
- `resources/*`
- `results/final_bam/*`
- `runs/mutect2_cjd_dilutions_with_pon/*`
- `results/mutect2_cjd_dilutions_with_pon/*`
- `prnp-junctions/*` (legacy comparison context)
- `doc/*` and `JOURNAL.md`

### Files/directories expected to create/modify
- `src/ddPCR/run_ddpcr.sh` (if missing)
- `src/pipelines/run_cjd_dilutions_variant_qc_with_pon.sh` (new wrapper)
- `src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R` (temporary AAF bypass)
- `src/ddPCR/README.md`
- `src/junctions/README.md`
- `src/pipelines/README.md`
- `bin/verify_output_checksums.sh`
- `doc/reproducibility/final_outputs_manifest.tsv`
- `doc/reproducibility/tooling_and_reference_provenance.md`
- `.gitignore` (if ignore coverage gaps are found)
- `doc/Codex_to_do/*` (planning/manifest artifacts)

### Ignore-rule files to inspect/update explicitly
Inspect:
- `.gitignore` (primary repo-wide ignore policy)
- `ddPCR/.gitignore` (raw ddPCR CSV policy)
- `prnp-junctions/.gitignore` (legacy junction policy; inspect for comparison only)
- `manuscript/.gitignore` (ensure manuscript-local ignores do not conflict)
- `manuscript/legacy/.gitignore` (legacy manuscript ignores)

Potentially update (only if needed for this task):
- `.gitignore`
- `ddPCR/.gitignore`

### Shell commands expected to run
- Discovery/read-only:
  - `ls`, `find`, `rg --files`, `rg -n`, `sed -n`, `cat`
- Dependency checks:
  - `command -v Rscript bwa samtools bedtools`
  - `Rscript -e 'requireNamespace(...)'`
  - `Rscript -e 'sessionInfo()'`
- Junction run:
  - `bash src/junctions/run_junctions.sh`
- Stage-12 regeneration:
  - `bash src/pipelines/12_cjd_dilutions_variant_qc_with_pon.sh`
- Output verification:
  - `find results ...`, `sha256sum ...`, `diff`/`cmp` on selected artifacts

### Whether network access is needed (and why)
- Preflight dependency detection: **No** network required.
- If packages/tools are missing and install is later authorized: **Yes**, for:
  - Conda channel package retrieval
  - CRAN/Bioconductor package retrieval (if conda packages unavailable)

### Potentially risky commands and safer alternatives
- Risk: deleting/overwriting outputs with `rm` or in-place cleanup.
  - Safer alternative: write to timestamped run/output directories and compare.
- Risk: destructive git commands (`git reset --hard`, checkout force).
  - Safer alternative: explicit targeted edits with patching and review.
- Risk: broad `find ... -delete`.
  - Safer alternative: no deletion in this task; keep additive outputs and compare via manifest/checksums.

## 3) Dependency-check plan for junction workflow

### Tools/packages to check
- CLI tools:
  - `Rscript`
  - `bwa`
  - `samtools`
  - `bedtools`
- R packages used directly by `src/junctions/*.R`:
  - `GenomicFeatures`
  - `GenomicRanges`
  - `GenomeInfoDb`
  - `Rsamtools`
  - `Biostrings`
  - `rtracklayer`
  - `GenomicAlignments`
  - `dplyr`
  - `txdbmaker` (required by recent BioC when `makeTxDbFromGFF()` is called)
- File prerequisites:
  - `resources/Homo_sapiens.GRCh38.110.gtf.gz`
  - `resources/hg38.fa`
  - `results/final_bam/*.bam`
  - `config/junctions.env` or defaults

### How installed vs missing will be detected
- CLI tools:
  - `command -v <tool>` and non-zero exit indicates missing.
- R packages:
  - `Rscript -e 'pkgs <- c(...); print(pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)])'`
- Inputs:
  - `test -f` / `test -d` and glob count checks for BAM inputs.

### If missing
- Stop immediately and report missing components with exact names.
- Do not install automatically.
- Wait for explicit authorization before any package/tool installation.

## 4) Run plan (junction workflow + QC regeneration step)

### Junction workflow run plan
1. Confirm prerequisites pass (tools, R packages, required files, BAM presence).
2. Execute:
   - `bash src/junctions/run_junctions.sh`
3. Expected outputs:
   - `resources/junctions/prnp_junctions.fa`
   - `resources/junctions/PRNP.pad1kb.hg38.bed`
   - `results/junctions/junction_align/*.PRNP.qname.bam`
   - `results/junctions/junction_align/*.PRNP.R1.fq.gz`
   - `results/junctions/junction_align/*.PRNP.R2.fq.gz`
   - `results/junctions/junction_align/*.PRNP.toJunc.bam`
   - `results/junctions/junction_align/*.PRNP.toJunc.bam.bai`
   - `results/junctions/junction_counts/prnp_junction_counts.tsv`
   - `results/junctions/junction_counts/prnp_junction_summary.tsv`

### Stage-12 QC regeneration run plan
1. Apply minimal temporary edit in:
   - `src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R`
   - Bypass/comment out AAF-filter application block (retain traceable marker comment).
2. Verify upstream inputs exist:
   - `runs/mutect2_cjd_dilutions_with_pon/<group>/annot_with_gnomad/*.vcf.gz`
   - `runs/mutect2_cjd_dilutions_with_pon/<group>/readcount_qc/metrics/*.metrics`
   - `resources/annotations/manual_population_freq.tsv`
3. Execute:
   - `bash src/pipelines/12_cjd_dilutions_variant_qc_with_pon.sh`
4. Expected outputs (per `cjd` and `dilutions`):
   - `results/mutect2_cjd_dilutions_with_pon/variant_qc/<group>/summary_combined_variants.tsv`
   - `results/mutect2_cjd_dilutions_with_pon/variant_qc/<group>/filtered_variants.tsv`
   - `results/mutect2_cjd_dilutions_with_pon/variant_qc/<group>/filtered_prnp_variants.tsv`
   - `results/mutect2_cjd_dilutions_with_pon/variant_qc/<group>/filtered_out_variants.tsv`
   - `results/mutect2_cjd_dilutions_with_pon/variant_qc/<group>/filter_counts.tsv`
   - `results/mutect2_cjd_dilutions_with_pon/variant_qc/<group>/run_settings.tsv`
   - `results/mutect2_cjd_dilutions_with_pon/variant_qc/<group>/final_withPoN_variants.tsv`

## 5) Assumptions/unknowns that could block execution

1. Ambiguity in “every workflow” scope:
   - assumed to mean at least `ddPCR`, `junctions`, and the main `pipelines` publication path.
2. Stage-12 regeneration depends on prior pipeline outputs existing under `runs/mutect2_cjd_dilutions_with_pon/`; if missing, regeneration cannot proceed.
3. Junction run depends on `results/final_bam/*.bam`; if absent/incomplete, junction outputs will be partial or empty.
4. R/Bioconductor version skew may require `txdbmaker` even when core packages are installed.
5. Some reference provenance fields may be incomplete in current docs and will require manual confirmation from you.
6. If dependencies are missing and installation is required, network policy/permissions may become a blocker until explicitly approved.

---

Preflight manifest complete. Stop here and wait for approval before execution.
