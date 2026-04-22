#!/usr/bin/env python3
"""Apply transparent post-processing filters to manual PRNP ORR cohort summaries."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply explicit, auditable filters to a manual PRNP ORR cohort "
            "summary without changing the underlying raw read summaries."
        )
    )
    parser.add_argument(
        "--input-summary",
        required=True,
        type=Path,
        help="Input cohort_summary.tsv from 07_run_manual_mosaic_prnp_orr_cohort.py.",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        type=Path,
        help="Prefix for filtered outputs, e.g. results/repeats/manual_cohort/cjd/filtered/default",
    )
    parser.add_argument(
        "--background-summary",
        default=None,
        type=Path,
        help="Optional background cohort_summary.tsv, typically the controls.",
    )
    parser.add_argument(
        "--min-exact-nonreference-reads",
        default=2,
        type=int,
        help="Minimum exact two-sided nonreference reads required unless synthetic support is stronger.",
    )
    parser.add_argument(
        "--min-synthetic-nonreference-reads",
        default=2,
        type=int,
        help="Minimum synthetic high/medium-confidence nonreference reads required for a synthetic-driven pass.",
    )
    parser.add_argument(
        "--min-plus-strand-reads",
        default=1,
        type=int,
        help="Minimum plus-strand exact nonreference reads.",
    )
    parser.add_argument(
        "--min-minus-strand-reads",
        default=1,
        type=int,
        help="Minimum minus-strand exact nonreference reads.",
    )
    parser.add_argument(
        "--min-unique-start-sites",
        default=2,
        type=int,
        help="Minimum unique exact nonreference start sites.",
    )
    parser.add_argument(
        "--max-one-sided-suspicious-reads",
        default=25,
        type=int,
        help="Maximum one-sided indel/soft-clip reads tolerated before a sample is filtered out.",
    )
    parser.add_argument(
        "--require-background-exceedance",
        action="store_true",
        help="Require a sample to exceed the provided background maxima on at least one selected metric.",
    )
    parser.add_argument(
        "--background-metrics",
        default="exact_nonreference_reads,synthetic_high_or_medium_nonreference_reads,exact_nonreference_unique_start_sites",
        help=(
            "Comma-separated metrics that can satisfy the background-exceedance test. "
            "Supported: exact_nonreference_reads, synthetic_high_or_medium_nonreference_reads, "
            "exact_nonreference_unique_start_sites, one_sided_indel_or_softclip_reads"
        ),
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Small TSV helpers
# ---------------------------------------------------------------------------

def read_tsv(path: Path) -> list[dict[str, str]]:
    # The filter works entirely from TSV exports so it stays easy to rerun and
    # audit without reopening BAM-level evidence.
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    # Create the destination directory automatically so new filter labels can be
    # written into a clean output tree.
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def parse_int(row: dict[str, str], field: str) -> int:
    # Treat blank fields as zero so threshold comparisons stay explicit and do
    # not need repetitive missing-value handling.
    value = row.get(field, "")
    return int(value) if value else 0


def encode_counter(counter: Counter[str]) -> str:
    # Keep overview counters compact while preserving the underlying labels.
    if not counter:
        return ""
    return ",".join(f"{label}:{counter[label]}" for label in sorted(counter))


# ---------------------------------------------------------------------------
# Background calibration helpers
# ---------------------------------------------------------------------------

def load_background_maxima(path: Path | None, metrics: list[str]) -> dict[str, int]:
    # The control cohort is summarized only by per-metric maxima here. That is
    # deliberate: the filter stays easy to audit because every pass/fail reason
    # can be traced to a simple threshold comparison.
    maxima = {metric: 0 for metric in metrics}
    if path is None or not path.exists():
        return maxima
    rows = read_tsv(path)
    for metric in metrics:
        maxima[metric] = max(parse_int(row, metric) for row in rows) if rows else 0
    return maxima


def background_exceedance_flags(
    row: dict[str, str],
    maxima: dict[str, int],
    metrics: list[str],
) -> tuple[bool, str]:
    # Record exactly which metrics beat the control maxima so later review does
    # not have to reconstruct why a sample counted as background-exceeding.
    exceeded = [
        f"{metric}>{maxima[metric]}"
        for metric in metrics
        if parse_int(row, metric) > maxima[metric]
    ]
    return (bool(exceeded), ",".join(exceeded))


def annotate_row(
    row: dict[str, str],
    args: argparse.Namespace,
    background_maxima: dict[str, int],
    background_metrics: list[str],
) -> dict[str, str]:
    # Keep the filter logic fully explicit. Each rule gets its own pass/fail
    # column so the post-processed table remains a transparent overlay on top of
    # the raw cohort summary instead of replacing it.
    exact_nonref = parse_int(row, "exact_nonreference_reads")
    synthetic_nonref = parse_int(row, "synthetic_high_or_medium_nonreference_reads")
    plus_reads = parse_int(row, "exact_nonreference_plus_reads")
    minus_reads = parse_int(row, "exact_nonreference_minus_reads")
    unique_sites = parse_int(row, "exact_nonreference_unique_start_sites")
    one_sided_suspicious = parse_int(row, "one_sided_indel_or_softclip_reads")

    signal_by_exact = exact_nonref >= args.min_exact_nonreference_reads
    signal_by_synthetic = synthetic_nonref >= args.min_synthetic_nonreference_reads
    has_min_signal = signal_by_exact or signal_by_synthetic

    # Exact-read evidence is expected to be bidirectional and supported by more
    # than one start site. Synthetic-only signal can still pass the early gate,
    # but exact-read artefacts remain heavily constrained.
    passes_bidirectional = (
        exact_nonref == 0
        or (
            plus_reads >= args.min_plus_strand_reads
            and minus_reads >= args.min_minus_strand_reads
        )
    )
    passes_unique_sites = exact_nonref == 0 or unique_sites >= args.min_unique_start_sites
    passes_one_sided_cap = one_sided_suspicious <= args.max_one_sided_suspicious_reads

    exceeds_background, exceeded_background_metrics = background_exceedance_flags(
        row=row,
        maxima=background_maxima,
        metrics=background_metrics,
    )
    # Background exceedance is optional so the same script can be used either as
    # a within-cohort technical filter or as a control-aware case-vs-control
    # screen.
    passes_background = (
        not args.require_background_exceedance or exceeds_background
    )

    fail_reasons: list[str] = []
    if not has_min_signal:
        fail_reasons.append("below_min_signal")
    if not passes_bidirectional:
        fail_reasons.append("fails_bidirectional_support")
    if not passes_unique_sites:
        fail_reasons.append("fails_unique_start_sites")
    if not passes_one_sided_cap:
        fail_reasons.append("too_many_one_sided_suspicious_reads")
    if not passes_background:
        fail_reasons.append("does_not_exceed_background")

    annotated = {
        **row,
        "filter_signal_by_exact": "yes" if signal_by_exact else "no",
        "filter_signal_by_synthetic": "yes" if signal_by_synthetic else "no",
        "filter_has_min_signal": "yes" if has_min_signal else "no",
        "filter_passes_bidirectional": "yes" if passes_bidirectional else "no",
        "filter_passes_unique_sites": "yes" if passes_unique_sites else "no",
        "filter_passes_one_sided_cap": "yes" if passes_one_sided_cap else "no",
        "filter_background_required": "yes" if args.require_background_exceedance else "no",
        "filter_exceeds_background": "yes" if exceeds_background else "no",
        "filter_exceeded_background_metrics": exceeded_background_metrics,
        "filter_passes_background": "yes" if passes_background else "no",
        "filter_pass": "yes" if not fail_reasons else "no",
        "filter_fail_reasons": ",".join(fail_reasons) if fail_reasons else "none",
    }
    return annotated


def build_overview(rows: list[dict[str, str]]) -> dict[str, str]:
    # Summarize both pass counts and the dominant reasons for rejection. This
    # makes threshold tuning easier without reopening the full annotated table.
    pass_counts = Counter(row["filter_pass"] for row in rows)
    fail_reason_counts = Counter()
    for row in rows:
        reasons = row["filter_fail_reasons"]
        if reasons and reasons != "none":
            for reason in reasons.split(","):
                fail_reason_counts[reason] += 1
    return {
        "sample_count": str(len(rows)),
        "filter_pass_counts": encode_counter(pass_counts),
        "filter_fail_reason_counts": encode_counter(fail_reason_counts),
    }


def main() -> int:
    args = parse_args()
    rows = read_tsv(args.input_summary)

    # Parse the background metric list once up front so both threshold loading
    # and row annotation share the exact same metric set.
    background_metrics = [
        metric.strip()
        for metric in args.background_metrics.split(",")
        if metric.strip()
    ]
    background_maxima = load_background_maxima(args.background_summary, background_metrics)

    # Annotate every row, then keep a candidate-only subset as a convenience
    # view rather than baking selection into the raw summary table.
    annotated_rows = [
        annotate_row(
            row=row,
            args=args,
            background_maxima=background_maxima,
            background_metrics=background_metrics,
        )
        for row in rows
    ]
    annotated_rows.sort(
        key=lambda row: (
            # Keep surviving candidates first, then preserve the original manual
            # review priority ordering within the annotated table.
            row["filter_pass"] != "yes",
            row["manual_review_priority"],
            -parse_int(row, "synthetic_high_or_medium_nonreference_reads"),
            -parse_int(row, "exact_nonreference_reads"),
            row["sample_id"],
        )
    )

    pass_rows = [row for row in annotated_rows if row["filter_pass"] == "yes"]

    filtered_path = args.output_prefix.with_suffix(".annotated.tsv")
    candidates_path = args.output_prefix.with_suffix(".candidates.tsv")
    overview_path = args.output_prefix.with_suffix(".overview.tsv")

    fieldnames = list(annotated_rows[0].keys()) if annotated_rows else []
    # Emit the full annotated table plus a candidate-only subset so downstream
    # review can stay focused without losing the rejected rows or fail reasons.
    write_tsv(filtered_path, annotated_rows, fieldnames)
    write_tsv(candidates_path, pass_rows, fieldnames)
    write_tsv(overview_path, [build_overview(annotated_rows)], ["sample_count", "filter_pass_counts", "filter_fail_reason_counts"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
