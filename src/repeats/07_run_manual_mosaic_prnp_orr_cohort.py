#!/usr/bin/env python3
"""Run the manual PRNP ORR mosaic reviewer across a selected BAM cohort."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import Counter
from pathlib import Path


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the manual PRNP ORR mosaic review helper across selected "
            "samples and aggregate one sortable cohort TSV."
        )
    )
    parser.add_argument(
        "--bam-dir",
        default=Path("results/final_bam"),
        type=Path,
        help="Directory containing input BAMs.",
    )
    parser.add_argument(
        "--reference-fasta",
        default=Path("resources/chr2_chr4_chr20.fasta"),
        type=Path,
        help="Reference FASTA used by the alignments.",
    )
    parser.add_argument(
        "--manual-script",
        default=Path("src/repeats/06_manual_mosaic_prnp_orr.py"),
        type=Path,
        help="One-sample manual review script to execute.",
    )
    parser.add_argument(
        "--output-root",
        default=Path("results/repeats/manual_cohort"),
        type=Path,
        help="Base output directory for cohort runs.",
    )
    parser.add_argument(
        "--cohort",
        choices=("controls", "cjd", "all"),
        default="controls",
        help="Named cohort selection when --sample-id is not used.",
    )
    parser.add_argument(
        "--cohort-label",
        default="",
        help="Output folder label. Defaults to the cohort name.",
    )
    parser.add_argument(
        "--sample-id",
        action="append",
        default=[],
        help="Sample ID to run. May be repeated. Overrides --cohort when set.",
    )
    parser.add_argument(
        "--sample-calls-tsv",
        default=Path("results/repeats/sample_calls.tsv"),
        type=Path,
        help="Optional sample_calls.tsv for annotations.",
    )
    parser.add_argument(
        "--gangstr-calls-tsv",
        default=Path("results/repeats/gangstr_calls.tsv"),
        type=Path,
        help="Optional gangstr_calls.tsv for annotations.",
    )
    parser.add_argument(
        "--synthetic-remap-mode",
        default="candidate_only",
        choices=("off", "candidate_only", "all_two_sided_exact"),
        help="Synthetic remap mode to pass to the one-sample helper.",
    )
    parser.add_argument(
        "--synthetic-flank-size",
        default=40,
        type=int,
        help="Synthetic flank size to pass to the one-sample helper.",
    )
    parser.add_argument(
        "--min-mapq",
        default=20,
        type=int,
        help="Minimum mapping quality for the one-sample helper.",
    )
    parser.add_argument(
        "--anchor-probe",
        default=12,
        type=int,
        help="Anchor probe size for the one-sample helper.",
    )
    parser.add_argument(
        "--min-anchor-bases",
        default=8,
        type=int,
        help="Minimum anchor bases for the one-sample helper.",
    )
    parser.add_argument(
        "--nearby-window",
        default=6,
        type=int,
        help="Nearby indel/soft-clip window for the one-sample helper.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun samples even if per-sample outputs already exist.",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Small TSV helpers and cohort discovery
# ---------------------------------------------------------------------------

def read_tsv(path: Path) -> list[dict[str, str]]:
    # The cohort wrapper only needs plain TSV parsing, so keep the helper small
    # and dependency-free.
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    # Create the target directory on demand so ad hoc cohort labels can be
    # written into a fresh output tree without extra setup.
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def cohort_label(args: argparse.Namespace) -> str:
    # Preserve an explicit output label whenever the caller provides one so
    # repeated tuning runs can coexist under the same results root.
    if args.cohort_label:
        return args.cohort_label
    if args.sample_id:
        return "custom"
    return args.cohort


def sample_group(sample_id: str) -> str:
    # Derive the group directly from the sample ID prefix to avoid depending on
    # an additional manifest for quick manual follow-up runs.
    if sample_id.startswith("Ctrl"):
        return "control"
    if sample_id.startswith("CJD"):
        return "cjd"
    return "other"


def discover_bams(args: argparse.Namespace) -> list[tuple[str, Path]]:
    # Use BAM basenames as the cohort source of truth so the wrapper can be run
    # against the current live BAM directory without maintaining a separate
    # manifest just for manual follow-up work.
    all_bams = {path.stem: path for path in sorted(args.bam_dir.glob("*.bam"))}
    if args.sample_id:
        missing = [sample_id for sample_id in args.sample_id if sample_id not in all_bams]
        if missing:
            raise FileNotFoundError(f"Requested BAMs not found for sample IDs: {', '.join(missing)}")
        return [(sample_id, all_bams[sample_id]) for sample_id in args.sample_id]

    def keep(sample_id: str) -> bool:
        if args.cohort == "controls":
            return sample_id.startswith("Ctrl")
        if args.cohort == "cjd":
            return sample_id.startswith("CJD")
        return sample_id.startswith("Ctrl") or sample_id.startswith("CJD")

    return [(sample_id, path) for sample_id, path in all_bams.items() if keep(sample_id)]


def encode_counter(counter: Counter[str]) -> str:
    # Keep cohort-level count summaries compact but still human-readable.
    if not counter:
        return ""
    return ",".join(f"{label}:{counter[label]}" for label in sorted(counter))


# ---------------------------------------------------------------------------
# Sample scoring and one-row-per-sample aggregation
# ---------------------------------------------------------------------------

def review_priority(
    synthetic_nonref: int,
    exact_nonref: int,
    plus_reads: int,
    minus_reads: int,
    unique_starts: int,
    one_sided_suspicious: int,
) -> str:
    # This is deliberately a triage label, not a biological call. It ranks
    # samples by how much manual attention they deserve before any background-
    # aware cohort filter is applied.
    if synthetic_nonref >= 2 and plus_reads > 0 and minus_reads > 0 and unique_starts >= 2:
        return "high"
    if exact_nonref >= 3 and plus_reads > 0 and minus_reads > 0 and unique_starts >= 2:
        return "medium"
    if synthetic_nonref >= 1 or exact_nonref >= 1 or one_sided_suspicious >= 5:
        return "low"
    return "background"


def summarize_sample(
    sample_id: str,
    reads_path: Path,
    summary_path: Path,
) -> dict[str, str]:
    # The per-sample cohort row is built from the detailed read table rather
    # than from the one-sample summary alone, because strand balance, unique
    # starts, and indel-pattern recurrence are central for mosaic triage.
    sample_summary = read_tsv(summary_path)[0]
    reads = read_tsv(reads_path)

    # Exact nonreference reads are the reads used for filter criteria (i), (ii),
    # and (iii): they span both ORR boundaries exactly and have a nonreference
    # ORR length signal in the original BAM alignment.
    exact_nonref_rows = [
        row
        for row in reads
        if row["anchor_tier"] == "two_sided_exact"
        and row["bam_signal"] in {"expansion_like", "contraction_like", "complex_length"}
    ]
    # Panel-scored reads are retained as an audit trail for direct ORR sequence
    # matching, but the final post-processing filter uses the synthetic count
    # below rather than this panel count.
    panel_nonref_rows = [
        row
        for row in reads
        if row["best_event_type"] in {"opri", "oprd"}
        and row["panel_confidence"] in {"high", "medium"}
    ]
    # Synthetic nonreference reads are the second route through criterion (i):
    # high/medium-confidence assignments to synthetic OPRI/OPRD alleles.
    synthetic_nonref_rows = [
        row
        for row in reads
        if row["synthetic_best_event_type"] in {"opri", "oprd"}
        and row["synthetic_remap_confidence"] in {"high", "medium"}
    ]
    synthetic_reference_rows = [
        row
        for row in reads
        if row["synthetic_best_event_type"] == "reference"
        and row["synthetic_remap_confidence"] in {"high", "medium"}
    ]

    # These exact-read summaries feed criteria (ii) and (iii): bidirectional
    # support and at least two unique exact-read start sites.
    plus_reads = sum(1 for row in exact_nonref_rows if row["strand"] == "+")
    minus_reads = sum(1 for row in exact_nonref_rows if row["strand"] == "-")
    unique_starts = len({row["reference_start_1based"] for row in exact_nonref_rows})

    # These counters are descriptive audit fields, helping reviewers see whether
    # a sample is driven by a repeated indel pattern or a specific length shift.
    unique_indel_ops = Counter(row["indel_ops"] for row in exact_nonref_rows if row["indel_ops"])
    delta_bp_counts = Counter(row["delta_bp_vs_reference"] for row in exact_nonref_rows if row["delta_bp_vs_reference"])
    synthetic_best_alleles = Counter(
        row["synthetic_best_allele_id"]
        for row in reads
        if row["synthetic_best_allele_id"]
    )
    synthetic_reasons = Counter(
        row["synthetic_remap_reason"]
        for row in reads
        if row["synthetic_remap_attempted"] == "yes"
    )
    panel_best_alleles = Counter(
        row["best_allele_id"]
        for row in reads
        if row["best_allele_id"] and row["best_event_type"] in {"opri", "oprd"}
    )

    total_exact_nonref = len(exact_nonref_rows)
    strand_balance = (
        f"{(min(plus_reads, minus_reads) / total_exact_nonref):.3f}"
        if total_exact_nonref
        else ""
    )

    # Keep both the raw counts and the condensed counters so the cohort TSV can
    # be sorted mechanically while still exposing the dominant artifact pattern.
    row = {
        "sample_id": sample_id,
        "group": sample_group(sample_id),
        "bam_path": sample_summary["bam_path"],
        "reads_tsv": str(reads_path),
        "summary_tsv": str(summary_path),
        "total_overlapping_primary_reads": sample_summary["total_overlapping_primary_reads"],
        "two_sided_exact_reads": sample_summary["two_sided_exact_reads"],
        "two_sided_window_reads": sample_summary["two_sided_window_reads"],
        "left_only_reads": sample_summary["left_only_reads"],
        "right_only_reads": sample_summary["right_only_reads"],
        "unanchored_reads": sample_summary["unanchored_reads"],
        "exact_nonreference_reads": str(total_exact_nonref),
        "exact_nonreference_plus_reads": str(plus_reads),
        "exact_nonreference_minus_reads": str(minus_reads),
        "exact_nonreference_unique_start_sites": str(unique_starts),
        "exact_nonreference_strand_balance": strand_balance,
        "exact_nonreference_delta_bp_counts": encode_counter(delta_bp_counts),
        "exact_nonreference_indel_ops": encode_counter(unique_indel_ops),
        "panel_high_or_medium_nonreference_reads": str(len(panel_nonref_rows)),
        "panel_nonreference_allele_counts": encode_counter(panel_best_alleles),
        "synthetic_remapped_reads": sample_summary["synthetic_remapped_reads"],
        "synthetic_high_or_medium_nonreference_reads": str(len(synthetic_nonref_rows)),
        "synthetic_reference_like_reads": str(len(synthetic_reference_rows)),
        "synthetic_best_allele_counts": encode_counter(synthetic_best_alleles),
        "synthetic_remap_reason_counts": encode_counter(synthetic_reasons),
        "synthetic_vs_panel_counts": sample_summary["synthetic_vs_panel_counts"],
        "one_sided_indel_or_softclip_reads": sample_summary["one_sided_indel_or_softclip_reads"],
        "manual_review_priority": review_priority(
            synthetic_nonref=len(synthetic_nonref_rows),
            exact_nonref=total_exact_nonref,
            plus_reads=plus_reads,
            minus_reads=minus_reads,
            unique_starts=unique_starts,
            one_sided_suspicious=int(sample_summary["one_sided_indel_or_softclip_reads"]),
        ),
        "eh_interpretation": sample_summary["eh_interpretation"],
        "eh_total_repeat_counts": sample_summary["eh_total_repeat_counts"],
        "gangstr_interpretation": sample_summary["gangstr_interpretation"],
        "gangstr_total_repeat_counts": sample_summary["gangstr_total_repeat_counts"],
    }
    return row


def build_overview(label: str, cohort: str, rows: list[dict[str, str]]) -> dict[str, str]:
    # The overview is intentionally compact: it answers "how noisy was this
    # cohort?" without opening the full sortable sample table.
    priority_counts = Counter(row["manual_review_priority"] for row in rows)
    synthetic_nonref_positive = sum(
        1 for row in rows if int(row["synthetic_high_or_medium_nonreference_reads"]) > 0
    )
    exact_nonref_positive = sum(
        1 for row in rows if int(row["exact_nonreference_reads"]) > 0
    )
    return {
        "cohort_label": label,
        "cohort": cohort,
        "sample_count": str(len(rows)),
        "priority_counts": encode_counter(priority_counts),
        "samples_with_exact_nonreference_reads": str(exact_nonref_positive),
        "samples_with_synthetic_high_or_medium_nonreference_reads": str(synthetic_nonref_positive),
    }


def main() -> int:
    args = parse_args()
    selected_bams = discover_bams(args)
    if not selected_bams:
        raise FileNotFoundError("No BAMs matched the requested cohort selection.")

    # Keep one subdirectory per cohort label so repeated exploratory runs do not
    # overwrite each other unless the caller explicitly reuses a label.
    label = cohort_label(args)
    cohort_root = args.output_root / label
    samples_root = cohort_root / "samples"
    samples_root.mkdir(parents=True, exist_ok=True)

    aggregated_rows: list[dict[str, str]] = []
    for sample_id, bam_path in selected_bams:
        output_prefix = samples_root / sample_id
        reads_path = output_prefix.with_suffix(".reads.tsv")
        summary_path = output_prefix.with_suffix(".summary.tsv")

        # Reuse prior per-sample outputs by default so threshold tweaks and
        # cohort regrouping do not force a full rerun of the one-sample helper.
        if args.force or not (reads_path.exists() and summary_path.exists()):
            command = [
                sys.executable,
                str(args.manual_script),
                "--bam",
                str(bam_path),
                "--sample-id",
                sample_id,
                "--reference-fasta",
                str(args.reference_fasta),
                "--output-prefix",
                str(output_prefix),
                "--synthetic-remap-mode",
                args.synthetic_remap_mode,
                "--synthetic-flank-size",
                str(args.synthetic_flank_size),
                "--min-mapq",
                str(args.min_mapq),
                "--anchor-probe",
                str(args.anchor_probe),
                "--min-anchor-bases",
                str(args.min_anchor_bases),
                "--nearby-window",
                str(args.nearby_window),
            ]
            if args.sample_calls_tsv.exists():
                command.extend(["--sample-calls-tsv", str(args.sample_calls_tsv)])
            if args.gangstr_calls_tsv.exists():
                command.extend(["--gangstr-calls-tsv", str(args.gangstr_calls_tsv)])
            subprocess.run(command, check=True)

        aggregated_rows.append(
            summarize_sample(
                sample_id=sample_id,
                reads_path=reads_path,
                summary_path=summary_path,
            )
        )

    # Keep the sortable cohort summary in manual-review order so the most
    # interesting samples surface first in spreadsheet-style inspection.
    aggregated_rows.sort(
        key=lambda row: (
            # High-priority rows float to the top, then stronger synthetic or
            # exact nonreference support, then sample ID for stable ordering.
            {"high": 0, "medium": 1, "low": 2, "background": 3}.get(
                row["manual_review_priority"], 4
            ),
            -int(row["synthetic_high_or_medium_nonreference_reads"]),
            -int(row["exact_nonreference_reads"]),
            row["sample_id"],
        )
    )

    cohort_summary_path = cohort_root / "cohort_summary.tsv"
    cohort_overview_path = cohort_root / "cohort_overview.tsv"

    # Keep the schemas explicit because the downstream filter consumes these
    # exact columns without additional metadata files.
    fieldnames = [
        "sample_id",
        "group",
        "bam_path",
        "reads_tsv",
        "summary_tsv",
        "total_overlapping_primary_reads",
        "two_sided_exact_reads",
        "two_sided_window_reads",
        "left_only_reads",
        "right_only_reads",
        "unanchored_reads",
        "exact_nonreference_reads",
        "exact_nonreference_plus_reads",
        "exact_nonreference_minus_reads",
        "exact_nonreference_unique_start_sites",
        "exact_nonreference_strand_balance",
        "exact_nonreference_delta_bp_counts",
        "exact_nonreference_indel_ops",
        "panel_high_or_medium_nonreference_reads",
        "panel_nonreference_allele_counts",
        "synthetic_remapped_reads",
        "synthetic_high_or_medium_nonreference_reads",
        "synthetic_reference_like_reads",
        "synthetic_best_allele_counts",
        "synthetic_remap_reason_counts",
        "synthetic_vs_panel_counts",
        "one_sided_indel_or_softclip_reads",
        "manual_review_priority",
        "eh_interpretation",
        "eh_total_repeat_counts",
        "gangstr_interpretation",
        "gangstr_total_repeat_counts",
    ]
    overview_fieldnames = [
        "cohort_label",
        "cohort",
        "sample_count",
        "priority_counts",
        "samples_with_exact_nonreference_reads",
        "samples_with_synthetic_high_or_medium_nonreference_reads",
    ]

    # Write both the sortable sample table and a one-row overview so downstream
    # filters can consume the detailed metrics while humans can skim cohort
    # noise levels quickly.
    write_tsv(cohort_summary_path, aggregated_rows, fieldnames)
    write_tsv(
        cohort_overview_path,
        [build_overview(label=label, cohort=args.cohort, rows=aggregated_rows)],
        overview_fieldnames,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
