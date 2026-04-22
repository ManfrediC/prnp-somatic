# PRNP ORR Repeat Workflow

This directory contains the repeat-aware PRNP octapeptide repeat region (ORR)
workflow for screening OPRI and OPRD candidates from the existing short-read
BAMs in `results/final_bam`.

## Inputs

- `results/final_bam/CJD*.bam`
- `results/final_bam/Ctrl*.bam`
- `resources/prnp_orr.hg38.bed`
- `resources/repeats/prnp_orr.expansionhunter.json`
- `resources/repeats/prnp_orr.gangstr.bed` (optional orthogonal caller target)
- `resources/chr2_chr4_chr20.fasta`

## Main entrypoint

From repository root:

```bash
conda activate prnp-repeats
bash src/repeats/run_prnp_orr.sh
```

If `REViewer` is installed in its own environment, set
`REVIEWER_CONDA_ENV=prnp-reviewer` in `config/repeats.env` and the workflow
will invoke it through `conda run`.

If `GangSTR` is installed locally, set `RUN_GANGSTR=1` in `config/repeats.env`.
The workflow will run GangSTR in `--targeted --nonuniform` mode against the
mutable PRNP ORR repeat block and summarize its VCF output as an orthogonal
local evidence layer.

Optional local overrides:

```bash
cp config/repeats.env.example config/repeats.env
bash src/repeats/run_prnp_orr.sh
```

For a clean rerun that archives the current live outputs first:

```bash
cp config/repeats.env.example config/repeats.env
printf 'ARCHIVE_EXISTING_RUN=1\nARCHIVE_RUN_LABEL="pre_rerun_2026-04-18"\n' >> config/repeats.env
bash src/repeats/run_prnp_orr.sh
```

## Outputs

- `results/repeats/sample_manifest.tsv`
- `results/repeats/sample_calls.tsv`
- `results/repeats/sample_review.tsv`
- `results/repeats/candidate_calls.tsv`
- `results/repeats/cohort_summary.tsv`
- `results/repeats/subclonal_read_support.tsv`
- `results/repeats/gangstr_calls.tsv`
- `results/repeats/somatic_screen.tsv`
- `results/repeats/run_settings.tsv`
- `results/repeats/raw/expansionhunter/`
- `results/repeats/raw/gangstr/` (when enabled)
- `results/repeats/review/reviewer/`
- `results/repeats/logs/`
- `results/repeats/old_runs/`

For cross-machine work, the repository keeps the text-based repeat artefacts in
git:

- `results/repeats/raw/expansionhunter/**/*.json`
- `results/repeats/raw/expansionhunter/**/*.vcf`
- `results/repeats/raw/gangstr/**/*.vcf`
- `results/repeats/raw/gangstr/**/*.samplestats.tab`
- `results/repeats/review/reviewer/**/*.metrics.tsv`
- `results/repeats/review/reviewer/**/*.phasing.tsv`
- `results/repeats/logs/*.log`
- `results/repeats/manual/*.tsv`
- `results/repeats/manual_cohort/*/samples/*.tsv`

The heavy artefacts stay local:

- `results/repeats/raw/**/` realigned BAMs and indexes
- `results/repeats/review/reviewer/**/*.svg`
- `results/repeats/old_runs/`

## Manual Mosaic Review

For one-sample manual review of possible low-level somatic OPRI/OPRD evidence
from the existing alignment, use:

```bash
conda activate prnp-repeats
python src/repeats/manual_mosaic_prnp_orr.py \
  --bam results/final_bam/CJD1.bam \
  --reference-fasta resources/chr2_chr4_chr20.fasta \
  --output-prefix results/repeats/manual/CJD1 \
  --sample-calls-tsv results/repeats/sample_calls.tsv \
  --gangstr-calls-tsv results/repeats/gangstr_calls.tsv
```

This writes:

- `results/repeats/manual/<sample>.reads.tsv`
- `results/repeats/manual/<sample>.summary.tsv`

The helper stays conservative:

- only primary, mapped reads are inspected
- two-sided anchor support is required for high-confidence block-level calls
- one-sided reads are retained as lower-confidence manual-review evidence
- the synthetic panel is a local rescoring step over extracted PRNP ORR
  sequences, not a full BAM remap
- by default, only exact two-sided reads with non-reference length or nearby
  indel/soft-clip evidence are remapped to synthetic contigs; set
  `--synthetic-remap-mode all_two_sided_exact` for a broader comparison pass

Panel definitions live in `resources/repeats/prnp_orr_manual_panel.tsv`.
Where a browsable primary source gave an exact human block-level architecture,
the panel marks it as `published_exact`; where only the published copy number
was available, the panel uses an explicit `representative_copy_number` model.

For control-first cohort calibration:

```bash
conda activate prnp-repeats
python src/repeats/run_manual_mosaic_prnp_orr_cohort.py --cohort controls
python src/repeats/run_manual_mosaic_prnp_orr_cohort.py --cohort cjd
```

This writes:

- `results/repeats/manual_cohort/controls/cohort_summary.tsv`
- `results/repeats/manual_cohort/controls/cohort_overview.tsv`
- `results/repeats/manual_cohort/controls/samples/<sample>.reads.tsv`
- `results/repeats/manual_cohort/controls/samples/<sample>.summary.tsv`
- `results/repeats/manual_cohort/cjd/cohort_summary.tsv`
- `results/repeats/manual_cohort/cjd/cohort_overview.tsv`
- `results/repeats/manual_cohort/cjd/samples/<sample>.reads.tsv`
- `results/repeats/manual_cohort/cjd/samples/<sample>.summary.tsv`

The cohort summary adds sortable per-sample metrics such as:

- exact two-sided nonreference-read count
- plus/minus strand split and unique start-site count
- top recurring indel pattern among exact nonreference reads
- synthetic-remap nonreference count
- manual review priority (`background`, `low`, `medium`, `high`)

To apply a transparent post-processing filter to the CJD cohort after the
control background is available:

```bash
conda activate prnp-repeats
python src/repeats/filter_manual_mosaic_prnp_orr_cohort.py \
  --input-summary results/repeats/manual_cohort/cjd/cohort_summary.tsv \
  --background-summary results/repeats/manual_cohort/controls/cohort_summary.tsv \
  --require-background-exceedance \
  --output-prefix results/repeats/manual_cohort/cjd/filtered/default
```

This keeps the raw cohort summary unchanged and writes:

- `results/repeats/manual_cohort/cjd/filtered/default.annotated.tsv`
- `results/repeats/manual_cohort/cjd/filtered/default.candidates.tsv`
- `results/repeats/manual_cohort/cjd/filtered/default.overview.tsv`

The default filter is intentionally explicit and auditable:

- minimum signal from exact or synthetic nonreference reads
- bidirectional support for exact nonreference reads
- minimum unique start-site count for exact nonreference reads
- cap on one-sided indel/soft-clip burden
- optional requirement that a CJD sample exceed control maxima on selected metrics

Convenience wrappers are available from repository root:

```bash
make repeats_manual_controls
make repeats_manual_cjd
make repeats_manual_filter_cjd
```

## Guardrails

- Existing live outputs are never mixed into a new run.
- Set `ARCHIVE_EXISTING_RUN=1` to move the previous live run into `results/repeats/old_runs/` before starting again.
- The workflow verifies that the manifest, per-sample caller outputs, reviewer outputs, and summary tables all cover the expected 32-sample cohort before reporting success.

## Interpretation

- `ExpansionHunter` is the primary caller.
- `REViewer` images are generated for manual assessment.
- `subclonal_read_support.tsv` ranks weak non-reference read tails from the
  ExpansionHunter JSON histograms.
- `GangSTR` is an optional orthogonal local STR caller for the same locus, not
  a validated low-VAF somatic caller.
- Any non-reference or uncertain result should be considered a candidate only
  until orthogonal validation has been performed.
