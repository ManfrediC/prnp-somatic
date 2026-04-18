#!/usr/bin/env python3
"""Summarize low-level non-reference read support from ExpansionHunter JSONs."""

from __future__ import annotations

# The goal of this helper is not to make somatic calls. Instead, it extracts
# the per-repeat-length read histograms that ExpansionHunter stores in JSON
# outputs so we can inspect weak non-reference tails across the cohort.

import argparse
import csv
import json
import re
from pathlib import Path


COUNT_PATTERN = re.compile(r"\(([-0-9]+),\s*([0-9]+)\)")
REFERENCE_MIDDLE_REPEAT_COUNT = 2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect ExpansionHunter per-read support histograms for low-level "
            "non-reference PRNP ORR evidence."
        )
    )
    parser.add_argument("--sample-calls", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def parse_count_histogram(raw: str) -> dict[int, int]:
    # ExpansionHunter stores histograms as strings like "(1, 14), (2, 8503)".
    # Convert those into repeat_count -> read_count dictionaries.
    return {int(repeat_count): int(read_count) for repeat_count, read_count in COUNT_PATTERN.findall(raw)}


def encode_count_histogram(counts: dict[int, int]) -> str:
    if not counts:
        return ""
    return ",".join(f"{repeat_count}:{counts[repeat_count]}" for repeat_count in sorted(counts))


def sum_matching(counts: dict[int, int], predicate) -> int:
    return sum(read_count for repeat_count, read_count in counts.items() if predicate(repeat_count))


def max_supported_repeat(counts: dict[int, int], predicate) -> tuple[str, str]:
    supported = [(repeat_count, read_count) for repeat_count, read_count in counts.items() if predicate(repeat_count)]
    if not supported:
        return "", ""
    repeat_count, read_count = max(supported, key=lambda item: (item[1], item[0]))
    return str(repeat_count), str(read_count)


def read_sample_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_row(sample_row: dict[str, str]) -> dict[str, str]:
    sample_id = sample_row["sample_id"]
    json_path = Path(sample_row["json_path"])
    with json_path.open(encoding="utf-8") as handle:
        payload = json.load(handle)

    variant = payload["LocusResults"]["PRNP_ORR"]["Variants"]["PRNP_ORR_MIDDLE_R2_BLOCK"]

    spanning_counts = parse_count_histogram(variant["CountsOfSpanningReads"])
    flanking_counts = parse_count_histogram(variant["CountsOfFlankingReads"])
    inrepeat_counts = parse_count_histogram(variant["CountsOfInrepeatReads"])

    spanning_total = sum(spanning_counts.values())
    flanking_total = sum(flanking_counts.values())
    inrepeat_total = sum(inrepeat_counts.values())

    spanning_reference = spanning_counts.get(REFERENCE_MIDDLE_REPEAT_COUNT, 0)
    flanking_reference = flanking_counts.get(REFERENCE_MIDDLE_REPEAT_COUNT, 0)
    inrepeat_reference = inrepeat_counts.get(REFERENCE_MIDDLE_REPEAT_COUNT, 0)

    spanning_lt_ref = sum_matching(
        spanning_counts, lambda repeat_count: repeat_count < REFERENCE_MIDDLE_REPEAT_COUNT
    )
    spanning_gt_ref = sum_matching(
        spanning_counts, lambda repeat_count: repeat_count > REFERENCE_MIDDLE_REPEAT_COUNT
    )
    flanking_lt_ref = sum_matching(
        flanking_counts, lambda repeat_count: repeat_count < REFERENCE_MIDDLE_REPEAT_COUNT
    )
    flanking_gt_ref = sum_matching(
        flanking_counts, lambda repeat_count: repeat_count > REFERENCE_MIDDLE_REPEAT_COUNT
    )
    inrepeat_lt_ref = sum_matching(
        inrepeat_counts, lambda repeat_count: repeat_count < REFERENCE_MIDDLE_REPEAT_COUNT
    )
    inrepeat_gt_ref = sum_matching(
        inrepeat_counts, lambda repeat_count: repeat_count > REFERENCE_MIDDLE_REPEAT_COUNT
    )

    top_spanning_gt_ref_repeat, top_spanning_gt_ref_count = max_supported_repeat(
        spanning_counts, lambda repeat_count: repeat_count > REFERENCE_MIDDLE_REPEAT_COUNT
    )
    top_flanking_gt_ref_repeat, top_flanking_gt_ref_count = max_supported_repeat(
        flanking_counts, lambda repeat_count: repeat_count > REFERENCE_MIDDLE_REPEAT_COUNT
    )

    def fraction(numerator: int, denominator: int) -> str:
        return f"{(numerator / denominator):.8f}" if denominator else ""

    return {
        "sample_id": sample_id,
        "group": sample_row["group"],
        "vcf_interpretation": sample_row["interpretation"],
        "locus_coverage": sample_row["LC"],
        "spanning_counts": encode_count_histogram(spanning_counts),
        "spanning_total_reads": str(spanning_total),
        "spanning_reference_reads": str(spanning_reference),
        "spanning_contraction_like_reads": str(spanning_lt_ref),
        "spanning_expansion_like_reads": str(spanning_gt_ref),
        "spanning_contraction_like_fraction": fraction(spanning_lt_ref, spanning_total),
        "spanning_expansion_like_fraction": fraction(spanning_gt_ref, spanning_total),
        "top_spanning_expansion_like_repeat_count": top_spanning_gt_ref_repeat,
        "top_spanning_expansion_like_read_count": top_spanning_gt_ref_count,
        "flanking_counts": encode_count_histogram(flanking_counts),
        "flanking_total_reads": str(flanking_total),
        "flanking_reference_reads": str(flanking_reference),
        "flanking_contraction_like_reads": str(flanking_lt_ref),
        "flanking_expansion_like_reads": str(flanking_gt_ref),
        "flanking_contraction_like_fraction": fraction(flanking_lt_ref, flanking_total),
        "flanking_expansion_like_fraction": fraction(flanking_gt_ref, flanking_total),
        "top_flanking_expansion_like_repeat_count": top_flanking_gt_ref_repeat,
        "top_flanking_expansion_like_read_count": top_flanking_gt_ref_count,
        "inrepeat_counts": encode_count_histogram(inrepeat_counts),
        "inrepeat_total_reads": str(inrepeat_total),
        "inrepeat_reference_reads": str(inrepeat_reference),
        "inrepeat_contraction_like_reads": str(inrepeat_lt_ref),
        "inrepeat_expansion_like_reads": str(inrepeat_gt_ref),
    }


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def main() -> int:
    args = parse_args()
    sample_rows = read_sample_rows(args.sample_calls)

    # Sort by expansion-like spanning support first because those reads are the
    # most plausible starting point for reviewing weak OPRI candidates.
    output_rows = sorted(
        (build_row(sample_row) for sample_row in sample_rows),
        key=lambda row: (
            -int(row["spanning_expansion_like_reads"]),
            -float(row["spanning_expansion_like_fraction"] or 0.0),
            row["sample_id"],
        ),
    )

    fieldnames = [
        "sample_id",
        "group",
        "vcf_interpretation",
        "locus_coverage",
        "spanning_counts",
        "spanning_total_reads",
        "spanning_reference_reads",
        "spanning_contraction_like_reads",
        "spanning_expansion_like_reads",
        "spanning_contraction_like_fraction",
        "spanning_expansion_like_fraction",
        "top_spanning_expansion_like_repeat_count",
        "top_spanning_expansion_like_read_count",
        "flanking_counts",
        "flanking_total_reads",
        "flanking_reference_reads",
        "flanking_contraction_like_reads",
        "flanking_expansion_like_reads",
        "flanking_contraction_like_fraction",
        "flanking_expansion_like_fraction",
        "top_flanking_expansion_like_repeat_count",
        "top_flanking_expansion_like_read_count",
        "inrepeat_counts",
        "inrepeat_total_reads",
        "inrepeat_reference_reads",
        "inrepeat_contraction_like_reads",
        "inrepeat_expansion_like_reads",
    ]
    write_rows(args.output, output_rows, fieldnames)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
