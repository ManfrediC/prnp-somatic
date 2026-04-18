# PRNP ORR Repeat Workflow

This directory contains the repeat-aware PRNP octapeptide repeat region (ORR)
workflow for screening OPRI and OPRD candidates from the existing short-read
BAMs in `results/final_bam`.

## Inputs

- `results/final_bam/CJD*.bam`
- `results/final_bam/Ctrl*.bam`
- `resources/prnp_orr.hg38.bed`
- `resources/repeats/prnp_orr.expansionhunter.json`
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
- `results/repeats/run_settings.tsv`
- `results/repeats/raw/expansionhunter/`
- `results/repeats/review/reviewer/`
- `results/repeats/logs/`
- `results/repeats/old_runs/`

## Guardrails

- Existing live outputs are never mixed into a new run.
- Set `ARCHIVE_EXISTING_RUN=1` to move the previous live run into `results/repeats/old_runs/` before starting again.
- The workflow verifies that the manifest, per-sample caller outputs, reviewer outputs, and summary tables all cover the expected 32-sample cohort before reporting success.

## Interpretation

- `ExpansionHunter` is the primary caller.
- `REViewer` images are generated for manual assessment.
- Any non-reference or uncertain result should be considered a candidate only
  until orthogonal validation has been performed.
