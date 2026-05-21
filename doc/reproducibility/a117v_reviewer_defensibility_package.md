# A117V Reviewer Defensibility Package

This packet serves to clarify how display names for the A177V-Spike-in samples were chosen.

## Core Position

The remediation is a label and interpretation correction, not a data change.

- Raw sample IDs, FASTQ paths, BAM paths, read counts, VAFs, QC thresholds, and
  variant filtering rules were not changed.
- The former manuscript-facing labels `A117V 0.5% spike-in` and
  `A117V 1% spike-in` were replaced with preparation/observation-compatible
  labels:
  - `NA99A1_undil` / `NA99A1_undil_D02` -> `A117V low`
  - `NA995A05_undil` / `NA995A05_undil_G01` -> `A117V high`
- The conservative LoD anchor remains `A117V low`, observed at 35/4315 reads
  = 0.81% VAF.
- `A117V high` is no longer presented as a nominal 0.5% LoD sample. It is
  described as the higher-input sample by preparation calculation, observed at
  92/6961 reads = 1.32% VAF.

## Evidence Summary

1. The original preparation workbook indicates that the sample corresponding to
   `NA995A05_undil_G01` used 5% heterozygous A117V-source input. That implies
   an expected mutant VAF of roughly 2.4-2.5%, not a true 0.5% mutant VAF.
2. The sequencing submission names and sequencing metadata remain internally
   consistent. There is no evidence that the downstream analysis swapped BAMs
   or sample identities.
3. Legacy and current outputs consistently show a higher A117V signal in
   `NA995A05_undil` than in `NA99A1_undil`.
4. The remediation makes the manuscript more conservative: the lower observed
   VAF sample, `A117V low` at 0.81%, remains the LoD anchor.

Detailed private audit trail: `doc/reproducibility/a117v_spike_in_source_audit.md`.

## What Changed

Current registries now carry a `display_label` column for the A117V spike-ins.
Current scripts read or propagate this display label while keeping raw IDs for
file lookup and joins. Manuscript-facing tables and prose use `A117V high` and
`A117V low`. The canonical DOCX was edited with tracked changes authored by
`Manfredi`.

Detailed private implementation log:
`doc/plans/a117v_remediation_implemented_changes.md`.

## What Did Not Change

- No legacy files were edited.
- No raw external files were modified.
- No FASTQ/BAM/sample-ID provenance fields were renamed.
- No numeric A117V read counts, VAFs, depths, strand counts, or quality metrics
  were changed.
- No filtering thresholds were relaxed to preserve the A117V result.

## Suggested Reviewer Response

If asked directly:

> During final source-file reconciliation, we found that the historical label
> of one A117V spike-in sample did not accurately describe its preparation
> calculation. We therefore preserved the raw sequencing identifiers as
> provenance, but changed the manuscript-facing labels to `A117V low` and
> `A117V high`. The correction did not alter any reads, counts, VAFs, QC
> thresholds, or variant calls. It made the LoD interpretation more
> conservative, because the lower observed spike-in, `A117V low` at 0.81% VAF,
> remains the LoD anchor.

If misconduct or post hoc sanitisation is implied:

> The correction is explicitly auditable: the original raw identifiers are
> preserved, the previous labels are retained in the private audit trail, and
> the numeric data are unchanged. The manuscript no longer claims detection of
> a true nominal 0.5% A117V spike-in. Source files and checksum evidence can be
> supplied confidentially if the editor requires verification.

## Checksum Strategy

Purpose: provide tamper-evident anchoring without making private source files
public.

Recommended public item:

- Commit a small checksum manifest to GitHub, for example
  `doc/reproducibility/a117v_source_checksums.tsv`.

Recommended columns:

```text
artifact_id	source_class	sha256	bytes	mtime_utc	redacted_path	notes
```

Artifacts to hash:

- Original preparation workbook used for A117V spike-in calculation.
- Original sequencing submission workbook/metadata for the A117V spike-ins.
- Current registry files after remediation:
  - `authoritative_files/manifest.tsv`
  - `config/preprocessing_samples.tsv`
- Current manuscript-facing A117V table source.
- Canonical DOCX after tracked changes.

How this helps:

- A GitHub commit anchors the hash values at a public timestamp.
- If reviewers request confidential verification, the private files can be
  provided and independently hashed.
- Matching hashes show that the files supplied later are the same files whose
  hashes were committed.

Important limitation:

- Checksums prove file identity after the hash is published. They do not, by
  themselves, prove that no error or alteration occurred before the checksum
  commit. Present them as tamper-evidence and audit support, not as absolute
  proof of historical intent.

Lightweight command pattern:

```bash
sha256sum path/to/artifact > doc/reproducibility/a117v_source_checksums.tsv
```

For a final public anchor, commit only the checksum manifest, not the private
source files or this reviewer package.
