#!/usr/bin/env python3
"""Merge repeat-screen evidence into a lightweight somatic review table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a merged PRNP ORR somatic screening table."
    )
    parser.add_argument("--results-root", required=True, type=Path)
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Small TSV helpers
# ---------------------------------------------------------------------------

def read_tsv(path: Path) -> list[dict[str, str]]:
    # All repeat workflow tables are plain TSVs, so a tiny shared reader keeps
    # the merge step dependency-light and easy to audit.
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    # Create the target directory on demand so this summary can be rebuilt into
    # a clean results tree after partial reruns.
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


# ---------------------------------------------------------------------------
# Main merge flow
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    # Key every source table by sample ID so the final somatic screen row can
    # place orthogonal evidence side by side without repeated joins later.
    sample_calls_rows = {
        row["sample_id"]: row for row in read_tsv(args.results_root / "sample_calls.tsv")
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
        # Pull together the high-level EH call, subclonal EH support, and
        # GangSTR signal into one compact review row per sample.
        subclonal = subclonal_rows[sample_id]
        gangstr = gangstr_rows[sample_id]
        output_rows.append(
            {
                "sample_id": sample_id,
                "group": sample_call["group"],
                "eh_interpretation": sample_call["interpretation"],
                "eh_total_repeat_counts": sample_call["total_repeat_counts"],
                "eh_locus_coverage": sample_call["LC"],
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

    # Keep the export schema explicit so manual review notes and downstream
    # figure/table scripts can rely on stable column order and names.
    fieldnames = [
        "sample_id",
        "group",
        "eh_interpretation",
        "eh_total_repeat_counts",
        "eh_locus_coverage",
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
