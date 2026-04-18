#!/usr/bin/env python3
"""Merge repeat-screen evidence into a lightweight somatic review table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a merged PRNP ORR somatic screening table."
    )
    parser.add_argument("--results-root", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def main() -> int:
    args = parse_args()
    sample_calls_rows = {
        row["sample_id"]: row for row in read_tsv(args.results_root / "sample_calls.tsv")
    }
    sample_review_rows = {
        row["sample_id"]: row for row in read_tsv(args.results_root / "sample_review.tsv")
    }
    subclonal_rows = {
        row["sample_id"]: row
        for row in read_tsv(args.results_root / "subclonal_read_support.tsv")
    }
    gangstr_rows = {
        row["sample_id"]: row for row in read_tsv(args.results_root / "gangstr_calls.tsv")
    }

    output_rows: list[dict[str, str]] = []
    for sample_id, sample_call in sample_calls_rows.items():
        sample_review = sample_review_rows[sample_id]
        subclonal = subclonal_rows[sample_id]
        gangstr = gangstr_rows[sample_id]
        output_rows.append(
            {
                "sample_id": sample_id,
                "group": sample_call["group"],
                "eh_interpretation": sample_call["interpretation"],
                "eh_total_repeat_counts": sample_call["total_repeat_counts"],
                "eh_locus_coverage": sample_call["LC"],
                "review_artifact": sample_review["review_artifact"],
                "review_status": sample_review["review_status"],
                "spanning_contraction_like_reads": subclonal["spanning_contraction_like_reads"],
                "spanning_contraction_like_fraction": subclonal["spanning_contraction_like_fraction"],
                "spanning_expansion_like_reads": subclonal["spanning_expansion_like_reads"],
                "spanning_expansion_like_fraction": subclonal["spanning_expansion_like_fraction"],
                "flanking_contraction_like_reads": subclonal["flanking_contraction_like_reads"],
                "flanking_contraction_like_fraction": subclonal["flanking_contraction_like_fraction"],
                "flanking_expansion_like_reads": subclonal["flanking_expansion_like_reads"],
                "flanking_expansion_like_fraction": subclonal["flanking_expansion_like_fraction"],
                "gangstr_interpretation": gangstr["gangstr_interpretation"],
                "gangstr_vcf_status": gangstr["vcf_status"],
                "gangstr_filter": gangstr["filter"],
                "gangstr_total_repeat_counts": gangstr["total_repeat_counts"],
                "gangstr_enclosing_nonreference_reads": gangstr["enclosing_nonreference_reads"],
                "gangstr_enclosing_nonreference_fraction": gangstr["enclosing_nonreference_fraction"],
                "gangstr_top_enclosing_nonreference_repeat_count": gangstr["top_enclosing_nonreference_repeat_count"],
                "gangstr_top_enclosing_nonreference_read_count": gangstr["top_enclosing_nonreference_read_count"],
                "gangstr_flanking_nonreference_reads": gangstr["flanking_nonreference_reads"],
                "gangstr_flanking_nonreference_fraction": gangstr["flanking_nonreference_fraction"],
                "gangstr_top_flanking_nonreference_repeat_count": gangstr["top_flanking_nonreference_repeat_count"],
                "gangstr_top_flanking_nonreference_read_count": gangstr["top_flanking_nonreference_read_count"],
                "gangstr_review_required": gangstr["gangstr_review_required"],
                "gangstr_vcf_path": gangstr["vcf_path"],
                "manual_review_label": "not_reviewed",
                "manual_review_notes": "none",
            }
        )

    fieldnames = [
        "sample_id",
        "group",
        "eh_interpretation",
        "eh_total_repeat_counts",
        "eh_locus_coverage",
        "review_artifact",
        "review_status",
        "spanning_contraction_like_reads",
        "spanning_contraction_like_fraction",
        "spanning_expansion_like_reads",
        "spanning_expansion_like_fraction",
        "flanking_contraction_like_reads",
        "flanking_contraction_like_fraction",
        "flanking_expansion_like_reads",
        "flanking_expansion_like_fraction",
        "gangstr_interpretation",
        "gangstr_vcf_status",
        "gangstr_filter",
        "gangstr_total_repeat_counts",
        "gangstr_enclosing_nonreference_reads",
        "gangstr_enclosing_nonreference_fraction",
        "gangstr_top_enclosing_nonreference_repeat_count",
        "gangstr_top_enclosing_nonreference_read_count",
        "gangstr_flanking_nonreference_reads",
        "gangstr_flanking_nonreference_fraction",
        "gangstr_top_flanking_nonreference_repeat_count",
        "gangstr_top_flanking_nonreference_read_count",
        "gangstr_review_required",
        "gangstr_vcf_path",
        "manual_review_label",
        "manual_review_notes",
    ]
    write_tsv(args.results_root / "somatic_screen.tsv", output_rows, fieldnames)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
