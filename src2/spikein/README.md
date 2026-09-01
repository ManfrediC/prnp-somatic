# Spike-in marker recovery

This directory follows the reviewed plan in `doc/spikein_plan.md`. New code is
adapted from the existing numbered pipeline stages, with explicit sample roles
and outputs confined to `results2/spikein/`. Existing inputs are read-only.

## Current step

The sample manifest, source genotyping, candidate selection, pure-source
read-count collection, parsing and marker finalisation are implemented and have
run. The fixed set contains four other SNPs plus the separate A117V control.
Stage 6 raw mixture counting, stage 7 read-level recovery, stage 8 Mutect2
calling and stage 9 filtering have run for the high and low mixtures.
Source genotyping completed successfully on 2026-08-31;
all four VCFs and their indexes passed the script's readability and sample checks.
The full console log is `results2/spikein/logs/source_genotyping_console.log`;
VCFs, `run.log` and `run_settings.tsv` are under `results2/spikein/discovery/`.
Candidate selection found eight other SNPs plus A117V. Direct counts are available
for all nine sites in both pure samples; marker finalisation produced the fixed set.

The WT background rule is now ALT VAF at most 0.001, with no separate maximum
ALT-read count. All other thresholds remain unchanged, including minimum ALT
support in donor and mixtures. Stage 2 uses AD-based AF; direct-count review
uses ALT divided by summed A/C/G/T counts. WT ALT counts remain in the tables.

The stage-2 rerun is under `results2/spikein/discovery/candidates_wt_vaf_only/`.
Previous outputs, including `candidates_wt_alt4/`, retain their original settings.
Both candidate tables are byte-identical: eight provisional other SNPs plus A117V.
The updated `pure_source_qc_review.tsv` reuses stage 3's raw counts: four other
SNPs pass, four still fail base quality, and A117V passes separately. All reviewed
measurements match the previous review and raw counts; only WT AF was reapplied,
with other QC decisions retained. No eligibility decision changed.

This is a review, not a final informative marker set. Only stage 2 was rerun for
this threshold change; no BAM recount, later-stage implementation or mixture
analysis was performed. Log: `results2/spikein/logs/wt_vaf_only_rerun.log`.
The parser added in the side conversation applies no thresholds and is unaffected.

The manifest records the four user-confirmed roles. `A100_1to2` is the pure
heterozygous donor despite its name. Fractions are provenance metadata only;
the caller does not use them. High uses 10/209 and low remains unknown.

## First stage

Run through WSL from the repository root after `conda activate prnp-spikein`.
The default is a read-only preflight and dry run:

```bash
bash src2/spikein/1_sources_haplotypecaller.sh
```

Preflight checks all four roles, distinct files and sample IDs, BAM header/EOF
integrity with `samtools quickcheck`, coordinate-sort metadata, contig agreement
and available BAM index readability. The index check discards a single-position
query result; it does not measure marker recovery. The stage reads the existing
BAMs and BAIs directly. Each index must use the canonical `<sample>.bam.bai`
name beside its BAM; alternative `<sample>.bai` names are not accepted.
The stage stops if an index is missing or unreadable and never creates indexes.
Quickcheck is not a complete scan for internal BAM
corruption, and an index-readability check does not prove that an index is
current for every locus.

The printed calling commands use only pure donor and pure WT. Each source
produces a diploid HaplotypeCaller gVCF, followed by CombineGVCFs and
GenotypeGVCFs. Calling uses the full authoritative capture interval list and
the existing subset FASTA. HC per-alignment-start downsampling is disabled;
native PairHMM uses one thread for reproducible processing, and the two sources
run sequentially. A117V identity and marker eligibility are checked in the
later pure-sample counting/selection stage, not by this script.

When the user authorises actual source genotyping:

```bash
DRY_RUN=0 bash src2/spikein/1_sources_haplotypecaller.sh
```

Settings are deliberately limited:

| Setting | Default | Purpose |
|---|---|---|
| `DRY_RUN` | `1` | `1` checks and prints; `0` executes |
| `SAMPLES_TSV` | `src2/spikein/samples.tsv` | Explicit five-column manifest |
| `OUT_ROOT` | `results2/spikein/discovery` | Fresh stage directory below `results2/spikein` |
| `JAVA_MEM_GB` | `8` | Heap for one GATK process at a time |

The script does not source `config/preprocessing.env` or accept overrides of the
authoritative reference/capture paths. It refuses an existing output directory,
including during dry run; use a fresh directory for a later repeat run. There
is no automatic resume or overwrite. Output paths must not contain whitespace
because the installed GATK launcher splits Java options on whitespace.
No BAM copies, reheadering or input-side
index creation are performed. GATK reads the original BAM with an explicit
index path.

An actual run writes two source gVCFs, `combined.g.vcf.gz`, `joint.vcf.gz` and
their indexes. `run_settings.tsv` records input/source hashes and caller
settings; `run.log` records commands, versions, role metadata and Git state.
Temporary files remain in the stage's `work/` directory. Output validation
requires readable VCFs/indexes and exactly the expected pure sample names; an
empty joint variant set is allowed.

## Verification of this implementation step

Run the focused preflight/dry-run checks in WSL:

```bash
python -B src2/spikein/test_marker_recovery.py
```

The checks use a guard that fails if a variant caller is invoked. They cover
pure-only calls, a dry run without writes, missing/duplicate/invalid manifest
roles, reused BAMs, mismatched sample names, missing inputs/indexes and unsafe
or existing output paths. One empty synthetic BAM exercises the missing-index
failure; no project BAM is copied or changed. Temporary test fixtures remain
under `results2/spikein/work` and are removed when each check finishes.

The first real-input preflight reported older BAI timestamps for `A100_1to2`
and `NA995A05_undil`. With user authorisation, samtools rebuilt comparison
copies under `results2/spikein/work/index_check/`. Both copies are byte-for-byte
identical to their original BAIs. On the user's subsequent instruction, the
originals were moved to the repository's `legacy/bai/`, and the rebuilt BAIs
were moved into the same canonical `results/final_bam/*.bam.bai` locations.
The script reads these canonical files; it never searches the archive.
This was a one-off authorised replacement, not part of the pipeline.
The original comparison remains in `results2/spikein/logs/stage1_index_validation.tsv`;
the moves and rechecked hashes are recorded in
`results2/spikein/logs/bai_relocation.tsv`. BAM contents were not changed.
Actual calling and the script's output validation have now completed. Marker
eligibility and the direct-read A117V identity check remain for the next stage.

The previously proposed **VAF >= 0.0081** comparator has been superseded. The
spike-in analysis now derives its empirical LoD from the minimum technically
recovered marker VAF across the two mixtures, including A117V. This first stage
does not apply a VAF threshold or select informative markers.

## Candidate selection

`2_define_markers.py` follows the existing Python CSV/TSV conventions. It reads
the completed joint VCF through the installed `bcftools query` and requires
exactly the two pure sample names in the manifest. `bcftools norm` checks every
original REF allele against the reference; its output is discarded. Selection
uses the original, unsplit VCF, already restricted to capture targets by stage 1.
Stage 2 does not repeat the capture checks. It reads only the A117V allele
definition from the historical A117V CSV. No BAMs are read.

Run in the WSL environment:

```bash
python -B src2/spikein/2_define_markers.py
```

The default output directory is `results2/spikein/discovery/candidates/`.
Existing output directories are refused. `--output-dir` accepts a fresh directory
under `results2/spikein/` for a repeat run; input paths are fixed for this project.

Thresholds are named constants at the top of the script and follow the plan.
Selection requires an original biallelic SNV, adequate source DP/GQ, donor allele
balance consistent with its genotype, and confident WT reference genotype with
minimal ALT support. AF is calculated from AD. Population frequency and somatic
filters are not used. The unfiltered VCF's `FILTER=.` is recorded, not rejected.

The three outputs are:

- `candidate_markers.tsv`: eligible SNPs plus the separately labelled A117V control.
- `candidate_qc.tsv`: original VCF records, genotype fields and exclusion reasons.
- `run_settings.tsv`: thresholds, sample IDs, Git commit, command line and the
  stage-1 settings path/hash. Detailed provenance remains with stage 1.

The next counting script will read chromosome and position directly from
`candidate_markers.tsv`. Stage 2 no longer writes a separate site file. The
previous run's site file and detailed settings are retained as historical outputs.

`candidate_qc_pass` covers the genotype and allele checks, not technical
read QC. Every candidate has `validation_status=pending_pure_read_counts`.
A117V is always included for counting, even if its exact VCF record is absent;
unavailable control genotype fields remain NA. Its direct-read identity check
must pass before finalisation. Only A117V has a supplied rsID in the current data.

The run reconciled all 30 joint records and selected eight other SNPs plus A117V.
The three analytical files were byte-identical on repeat execution. No test
suite or hypothetical data fixtures were added. Logs:
`results2/spikein/logs/candidate_selection.log` and
`results2/spikein/logs/candidate_selection_validation.log`.
The shorter implementation also produced identical analytical files; its check
is recorded in `results2/spikein/logs/candidate_bcftools_validation.log`.
The latest simplification retained both analytical tables byte-for-byte and
produced only the three outputs listed above. Verification:
`results2/spikein/logs/candidate_minimal_validation.log`.

## Pure-source read counting

`3_fixed_site_readcount_qc.sh` adapts the existing stages 3 and 10. It reads the
nine sites from `results2/spikein/discovery/candidates_revised/candidate_markers.tsv`
and selects only `pure_donor` and `pure_wt` from the manifest. No mixture is read.

Run from the WSL `prnp-spikein` environment. The default prints commands without
creating files; set `DRY_RUN=0` for an authorised counting run:

```bash
bash src2/spikein/3_fixed_site_readcount_qc.sh
```

The current script runs bam-readcount directly on each original BAM.
bam-readcount excludes duplicate, unmapped, secondary and QC-failed reads in its
own pileup logic, so no filtered BAM or new index is needed. A shared `sites.tsv`
uses bam-readcount's 1-based inclusive coordinates. Explicit `-q 0 -b 0` retains
the existing counting convention; the depth cap is 10,000,000.

Raw counts go into `readcounts/<sample>.txt`; commands, input table hashes and
tool messages go into `run.log`. Tool versions are recorded with environment
setup in `doc/spikein_plan.md`. `OUT_ROOT` can name a fresh directory below
`results2/spikein/`; existing output directories are refused.

Each requested coordinate must be returned exactly once. Zero-depth rows are
accepted, but missing rows and reported depths at or above the cap stop the run.
The authorised run completed successfully on 2026-08-31 in 95 seconds using the
earlier duplicate-filtered-BAM implementation. It has not been rerun with the
direct-BAM script. Each pure
sample returned all nine sites exactly once. Reported depth was 4,486–7,473 in
donor and 3,848–6,544 in WT, well below the cap. Reference bases match the candidate
table after case normalisation; bam-readcount preserves lower-case FASTA bases.
The console log is `results2/spikein/logs/pure_readcount_console.log`.

Repeated warnings concern missing per-read `SM` tags, used only for optional
single-end mapping quality. Counts, ordinary mapping quality, base quality and
strand counts are unaffected. Omit that optional metric or mark it unavailable
in the later parser; see the source check in `doc/spikein_plan.md`.

The separate parser is implemented below. Technical QC and A117V identity
validation remain for the next step. No final marker set or mixture results
have been produced.

## Read-count conversion

`4_readcount_to_tsv.py` adapts the existing parser to the two pure samples in
the manifest. Run in WSL from the repository root:

```bash
python -B src2/spikein/4_readcount_to_tsv.py
```

It reads stage 3's raw files and writes one `<sample>_metrics.tsv` per sample
under `results2/spikein/readcount_qc/pure/metrics/`. An existing output directory
is refused; `--output-dir` can select a fresh directory under `results2/spikein`.

Each allele row contains `CHROM`, `POS`, `REF`, reported `DEPTH`, `ACGT_DEPTH`,
`BASE`, `COUNT`, `MEAN_BQ`, `MEAN_MQ`, `FWD` and `REV`. `ACGT_DEPTH` is the sum of
A/C/G/T observations for later SNV VAF calculation; it excludes N and indel
events. Reference and allele bases are upper case. Zero-count alleles are
retained with mean qualities recorded as `NA`.

Quality labels follow the
[bam-readcount 1.0.1 output definition](https://github.com/genome/bam-readcount/blob/v1.0.1/README.md#output):
the raw mapping-quality field precedes base quality. The unavailable optional
single-end mapping quality and unused metrics are omitted. This stage performs
no marker selection, VAF filtering or technical QC and does not read BAMs.

The conversion ran in WSL and produced 55 donor and 54 WT allele rows, covering
all nine sites in each sample. Every exported field was checked against the
raw counts. A repeat invocation refused to overwrite the tables. The run and
verification log is `results2/spikein/logs/readcount_conversion.log`.

## Marker finalisation

`5_marker_recovery.py` reads only the latest
provisional candidates and the two pure-source metrics tables. It independently
applies the agreed direct-count criteria, including WT AF at most 0.001 with no
WT ALT-count cap. The stage-2 genotype decision remains a prerequisite.

A117V must be the single labelled control, be heterozygous in the donor and pass
all pure-source checks. The completed run wrote `marker_qc.tsv`, its passing
subset `informative_markers.tsv`, and minimal settings under
`results2/spikein/markers/`. The informative table was hashed before any mixture
was read.

The audit contains nine provisional markers. Four other SNPs plus A117V form the
fixed informative set; four other SNPs fail mean base-quality criteria. The
informative table SHA-256 is
`56d5f63410fb338e610d14f200a9c9ea05507b11a3d6a1f92ac3d69c24557a8f`.
No mixture, caller or recovery-status logic ran. Log:
`results2/spikein/logs/marker_finalisation.log`.

## Mixture read counting

`6_mixtures_fixed_site_readcount_qc.sh` is a minimal adaptation of stage 3. It
verifies the frozen marker-table hash and checks the manifest against the hash
recorded by source genotyping. It then selects only the manifest's high and low
BAMs. bam-readcount reads those BAMs directly with explicit `-q 0 -b 0` and a
depth cap of 10,000,000; its own pileup logic excludes duplicate-flagged reads.

The default is a write-free dry run:

```bash
bash src2/spikein/6_mixtures_fixed_site_readcount_qc.sh
```

The authorised run completed on 2026-09-01 in 51.12 seconds and wrote `sites.tsv`
and raw counts under `results2/spikein/readcount_qc/mixtures/`. It checked
coordinates, reference bases and the depth cap. It did not parse counts or
assign any recovery status. Bash syntax, ShellCheck and the real-input dry run
passed; the dry-run log is
`results2/spikein/logs/mixture_readcount_direct_bam_dry_run.log`.

The script was then reorganised into small functions matching its existing
validation blocks. All comments and commands were retained. A saved pre-refactor
run and the canonical run produced byte-identical `sites.tsv` and raw count files
for both mixtures. Paired SHA-256 hashes are recorded in
`results2/spikein/logs/mixture_readcount_refactor_sha256.txt`. Ruff and Radon do
not analyse Bash; their direct invocation therefore supplies no valid code-quality
or complexity result for this stage.

## Mixture read recovery

`7_mixture_read_recovery.py` is drafted to evaluate the frozen informative-marker
set in the high and low mixture count tables. It applies the agreed depth,
support, strand, base-quality and mapping-quality thresholds and records expected
AF where the source fraction is established. After technical recovery is fixed,
it derives the empirical LoD as the minimum recovered marker VAF across both
mixtures, including A117V. The LoD is not applied back to the spike-in rows. The
script does not read BAMs or call variants.

The authorised run completed on 2026-09-01 and wrote
`mixture_read_recovery.tsv`, `empirical_lod.tsv` and `run_settings.tsv` under
`results2/spikein/read_recovery/`. All ten marker-by-mixture observations passed
the technical read-level criteria. The empirical LoD is 26/3,891 = 0.0066820869
(0.668%), supported by chr20:4693455 G>A in the low mixture. Log:
`results2/spikein/logs/mixture_read_recovery_console.log`.

## Mixture Mutect2 calling

`8_mixtures_mutect2_with_pon.sh` closely adapts the established stage 8. It
checks that stage 5's marker table is unchanged, then calls only the high and
low mixtures across the full capture panel. It uses the same panel of normals,
gnomAD germline resource and allele-frequency prior as the existing pipeline.
The final script sets `--initial-tumor-lod 0`, spike-in-only
`--max-population-af 1.0` and `--max-reads-per-alignment-start 0`. Emission LOD
3, activity-probability threshold 0.002 and all other caller parameters retain
their GATK defaults. Its default output is `results2/spikein/mutect2_final/`.

The default is a write-free dry run:

```bash
bash src2/spikein/8_mixtures_mutect2_with_pon.sh
```

The completed run wrote raw VCFs, statistics and orientation data under
`results2/spikein/mutect2/`. This stage does not run `FilterMutectCalls`, apply
somatic filters or change the stage-7 recovery result and empirical LoD. The
high and low VCFs contain 33 and 32 raw records, respectively. Mutect2 represents
A117V and the adjacent fixed marker as one `CA` to `TG` record in both mixtures;
later fixed-marker comparison must atomise this record before exact matching.

An earlier separate uncapped comparison under
`results2/spikein/mutect2_no_alignment_start_cap/` changed only the
alignment-start cap to 0. The high-mixture VCF increased from 33 to 63 raw
records and newly emitted chr20:4693455 G>A, giving three of five atomised
marker matches. The low-mixture
VCF retained the same 32 raw record keys and two of five atomised marker matches.
No filtering or further parameter tuning followed this comparison.

The capped `--max-population-af 1.0` run under
`results2/spikein/mutect2_max_population_af_1/` retained all 33 high and 32 low
canonical record keys. Marker emission was unchanged: only A117V and the
adjacent marker were present after atomisation in each mixture. Three non-marker
records had small caller-field changes. FilterMutectCalls was not run.

The combined run under
`results2/spikein/mutect2_max_population_af_1_no_alignment_start_cap/` produced
63 high and 32 low raw records. Its record keys and marker emissions match the
earlier uncapped run: three atomised markers in high and two in low. Raising
`max-population-af` did not add or remove a raw record when alignment-start
downsampling was disabled. FilterMutectCalls was not run.

Focused assembly-region diagnostics are under
`results2/spikein/mutect2_diagnostics/`. The two markers at chr20:4688888 and
chr20:4690735 are in inactive regions in both mixtures. chr20:4693455 is active
in high and inactive in low. `--genotype-germline-sites true` does not change
these states or marker emission. `--initial-tumor-lod 0` activates and emits
chr20:4693455 in low, but the first two markers remain inactive. Forced calling
shows strong later assembly evidence at all three missing sites, so their
ordinary non-emission occurs at the earlier activity decision. The complete
marker result is `mutect2_diagnostics/activity_profile_marker_summary.tsv`.

Matched full-panel runs then compared `--tumor-lod-to-emit 3` with 0 while
holding the other three settings fixed. Lowering the threshold increased the
high raw VCF from 63 to 67 records and the low raw VCF from 38 to 39 records.
It did not add a fixed marker. Both runs emitted chr20:4693455 G>A, A117V and
chr20:4699571 A>G in both mixtures; chr20:4688888 G>A and chr20:4690735 T>G
remained absent. The atomised comparison is
`results2/spikein/logs/mutect2_initial0_emit_lod_comparison.tsv`.

The emission-LOD 3 callset was then processed with the original project filter
settings. A117V and chr20:4699571 A>G passed in both mixtures. The emitted
chr20:4693455 G>A call was filtered for `clustered_events` in both mixtures.
The two inactive markers remained absent. The ten-row marker result is
`results2/spikein/logs/filtermutectcalls_initial0_emit3_marker_results.tsv`.

Rerunning FilterMutectCalls with only `--max-events-in-region 5` allowed
chr20:4693455 G>A to pass in both mixtures. A117V and chr20:4699571 A>G
continued to pass, while the two markers absent from the raw VCF remained
absent. PASS records increased from 6 to 13 in high and from 4 to 10 in low.
Of the added non-marker PASS alleles, one in high and two in low lacked an exact
allele match in the pure-source joint VCF. The concise comparison is
`results2/spikein/logs/filtermutectcalls_initial0_emit3_max_events_5_comparison.tsv`.

Full-panel caller runs compared activity-probability thresholds 0.002, 0.0002,
0.00002 and 0. All four produced the same 63 high and 38 low raw records, the
same atomised alleles and the same three emitted markers in each mixture. The
comparison is `results2/spikein/logs/mutect2_active_probability_comparison.tsv`.

## Mixture Mutect2 filtering

`9_mixtures_filtermutectcalls.sh` learns the established orientation model for
each final stage-8 call and runs FilterMutectCalls with
`--max-events-in-region 5`. All other filtering parameters retain their GATK
defaults. It writes through a temporary directory, validates both VCFs and
publishes the completed output as
`results2/spikein/filtermutectcalls/`.

The exploratory Mutect2 and FilterMutectCalls runs, logs and comparison script
are preserved under `legacy/spikein_mutect2_exploration_2026-09-01/`, with
their original repository-relative paths and SHA-256 hashes. Stage 9 no longer
recreates the historical four-run comparison.

The original and maximum-population-AF-1 runs produced 33/5 total/PASS records
in high and 32/4 in low. The two uncapped runs produced 63/6 in high and 32/3
in low. With the cap, the A117V/adjacent-marker record was filtered for
`strand_bias` in high and passed in low. Without the cap, that record passed in
both mixtures; chr20:4693455 G>A was also emitted in high but filtered for
`clustered_events`. The remaining fixed markers were not emitted.

After atomisation, no PASS non-marker allele lacked an exact pure-source
joint-VCF match in any run. Removing the cap increased high-mixture non-marker
alleles from 37 to 68 and alleles without a pure-source match from 15 to 45,
but all 45 were filtered. Raising maximum population AF did not change these
aggregate outcomes. The directory contains the eight complete filtered VCFs,
eight orientation models, one run log, settings and SHA-256 records.
