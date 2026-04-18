#!/usr/bin/env python3
"""Summarize per-sample GangSTR outputs for exploratory PRNP ORR screening."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


PAIR_PATTERN = re.compile(r"(-?\d+)\s*[:=,]\s*(-?\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize GangSTR VCF outputs for the PRNP ORR screen."
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--reference-total-repeats", required=True, type=int)
    parser.add_argument("--variable-repeat-offset", required=True, type=int)
    parser.add_argument("--gangstr-enabled", choices=("0", "1"), default="0")
    return parser.parse_args()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_info_field(field: str) -> dict[str, str]:
    info = {}
    for item in field.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
    return info


def split_allele_values(raw: str) -> list[str]:
    if raw in {"", "."}:
        return []
    return [part for part in re.split(r"[/|,]", raw) if part != ""]


def parse_repeat_counts(raw: str) -> list[int]:
    counts = []
    for token in split_allele_values(raw):
        try:
            counts.append(int(token))
        except ValueError:
            pass
    return counts


def parse_confidence_intervals(raw: str) -> list[tuple[int, int] | None]:
    intervals: list[tuple[int, int] | None] = []
    for token in split_allele_values(raw):
        match = re.fullmatch(r"(-?\d+)-(-?\d+)", token)
        if not match:
            intervals.append(None)
            continue
        intervals.append((int(match.group(1)), int(match.group(2))))
    return intervals


def parse_read_summary(raw: str) -> dict[int, int]:
    if raw in {"", "."}:
        return {}
    summary: dict[int, int] = {}
    for token in raw.split("|"):
        token = token.strip()
        if not token:
            continue
        match = PAIR_PATTERN.fullmatch(token)
        if not match:
            continue
        summary[int(match.group(1))] = int(match.group(2))
    return summary


def encode_histogram(counts: dict[int, int]) -> str:
    if not counts:
        return ""
    return ",".join(f"{repeat_count}:{counts[repeat_count]}" for repeat_count in sorted(counts))


def sum_matching(counts: dict[int, int], predicate) -> int:
    return sum(read_count for repeat_count, read_count in counts.items() if predicate(repeat_count))


def max_supported_repeat(counts: dict[int, int], predicate) -> tuple[str, str]:
    supported = [
        (repeat_count, read_count)
        for repeat_count, read_count in counts.items()
        if predicate(repeat_count)
    ]
    if not supported:
        return "", ""
    repeat_count, read_count = max(supported, key=lambda item: (item[1], item[0]))
    return str(repeat_count), str(read_count)


def parse_sample_vcf(vcf_path: Path) -> dict[str, str]:
    if not vcf_path.exists():
        return {"vcf_status": "missing"}

    header_sample = ""
    record_data: dict[str, str] = {"vcf_status": "ok"}
    with vcf_path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#CHROM"):
                cols = line.split("\t")
                if len(cols) > 9:
                    header_sample = cols[9]
                continue
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 10:
                record_data["vcf_status"] = "malformed"
                return record_data
            record_data.update(
                {
                    "sample_name": header_sample,
                    "chrom": cols[0],
                    "pos": cols[1],
                    "id": cols[2],
                    "ref_base": cols[3],
                    "alt": cols[4],
                    "filter": cols[6],
                    "info_raw": cols[7],
                    "format": cols[8],
                    "sample_field": cols[9],
                }
            )
            break
        else:
            record_data["vcf_status"] = "no_records"
            return record_data

    info = parse_info_field(record_data["info_raw"])
    record_data["info_END"] = info.get("END", "")
    record_data["info_PERIOD"] = info.get("PERIOD", "")
    record_data["info_RU"] = info.get("RU", "")
    record_data["info_REF"] = info.get("REF", "")

    keys = record_data["format"].split(":")
    values = record_data["sample_field"].split(":")
    fmt = dict(zip(keys, values))
    for key in ["GT", "DP", "Q", "REPCN", "REPCI", "RC", "ENCLREADS", "FLNKREADS", "ML"]:
        record_data[key] = fmt.get(key, "")
    return record_data


def classify_sample(
    record: dict[str, str],
    reference_total_repeats: int,
    variable_repeat_offset: int,
) -> dict[str, str]:
    repcn = parse_repeat_counts(record.get("REPCN", ""))
    repci = parse_confidence_intervals(record.get("REPCI", ""))
    total_counts = [count + variable_repeat_offset for count in repcn]

    interpretation = "disabled"
    if record.get("vcf_status") == "disabled":
        interpretation = "disabled"
    elif record.get("vcf_status") != "ok":
        interpretation = "failed_qc"
    elif len(repcn) != 2:
        interpretation = "failed_qc"
    else:
        has_wide_ci = False
        if len(repci) == len(repcn):
            for count, interval in zip(repcn, repci):
                if interval is None:
                    has_wide_ci = True
                    break
                if interval[0] != count or interval[1] != count:
                    has_wide_ci = True
        else:
            has_wide_ci = True

        if record.get("filter", "") not in {"PASS", "."}:
            interpretation = "uncertain"
        elif any(total > reference_total_repeats for total in total_counts):
            interpretation = "candidate_OPRI"
        elif any(total < reference_total_repeats for total in total_counts):
            interpretation = "candidate_OPRD"
        elif has_wide_ci:
            interpretation = "uncertain"
        else:
            interpretation = "reference"

    row = {
        "gangstr_interpretation": interpretation,
        "gangstr_review_required": "yes" if interpretation not in {"reference", "disabled"} else "no",
        "middle_repeat_counts": ",".join(str(value) for value in repcn),
        "total_repeat_counts": ",".join(str(value) for value in total_counts),
        "delta_vs_reference": ",".join(
            str(value - reference_total_repeats) for value in total_counts
        ),
    }
    for idx in range(2):
        middle = repcn[idx] if idx < len(repcn) else ""
        total = total_counts[idx] if idx < len(total_counts) else ""
        delta = (
            total_counts[idx] - reference_total_repeats
            if idx < len(total_counts)
            else ""
        )
        ci_raw = split_allele_values(record.get("REPCI", ""))
        ci = ci_raw[idx] if idx < len(ci_raw) else ""
        row[f"allele{idx + 1}_middle_repeat_count"] = str(middle)
        row[f"allele{idx + 1}_total_repeat_count"] = str(total)
        row[f"allele{idx + 1}_delta_vs_reference"] = str(delta)
        row[f"allele{idx + 1}_middle_repeat_ci"] = ci
    return row


def build_disabled_row(sample_id: str, group: str) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "group": group,
        "vcf_path": "",
        "vcf_status": "disabled",
        "chrom": "",
        "pos": "",
        "filter": "",
        "GT": "",
        "DP": "",
        "Q": "",
        "REPCN": "",
        "REPCI": "",
        "RC": "",
        "ENCLREADS": "",
        "FLNKREADS": "",
        "RU": "",
        "REF": "",
        "allele1_middle_repeat_count": "",
        "allele2_middle_repeat_count": "",
        "allele1_middle_repeat_ci": "",
        "allele2_middle_repeat_ci": "",
        "allele1_total_repeat_count": "",
        "allele2_total_repeat_count": "",
        "allele1_delta_vs_reference": "",
        "allele2_delta_vs_reference": "",
        "middle_repeat_counts": "",
        "total_repeat_counts": "",
        "delta_vs_reference": "",
        "gangstr_interpretation": "disabled",
        "gangstr_review_required": "no",
        "enclosing_counts": "",
        "enclosing_total_reads": "0",
        "enclosing_nonreference_reads": "0",
        "enclosing_nonreference_fraction": "",
        "top_enclosing_nonreference_repeat_count": "",
        "top_enclosing_nonreference_read_count": "",
        "flanking_counts": "",
        "flanking_total_reads": "0",
        "flanking_nonreference_reads": "0",
        "flanking_nonreference_fraction": "",
        "top_flanking_nonreference_repeat_count": "",
        "top_flanking_nonreference_read_count": "",
    }


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def main() -> int:
    args = parse_args()
    manifest_rows = read_manifest(args.manifest)
    raw_root = args.results_root / "raw" / "gangstr"
    gangstr_enabled = args.gangstr_enabled == "1"

    output_rows: list[dict[str, str]] = []
    for manifest_row in manifest_rows:
        sample_id = manifest_row["sample_id"]
        group = manifest_row["group"]
        if not gangstr_enabled:
            output_rows.append(build_disabled_row(sample_id, group))
            continue

        vcf_path = raw_root / sample_id / f"{sample_id}.vcf"
        record = parse_sample_vcf(vcf_path)
        classified = classify_sample(
            record,
            reference_total_repeats=args.reference_total_repeats,
            variable_repeat_offset=args.variable_repeat_offset,
        )

        enclosing_counts = parse_read_summary(record.get("ENCLREADS", ""))
        flanking_counts = parse_read_summary(record.get("FLNKREADS", ""))
        enclosing_total = sum(enclosing_counts.values())
        flanking_total = sum(flanking_counts.values())
        reference_middle_repeat_count = args.reference_total_repeats - args.variable_repeat_offset
        enclosing_nonreference = sum_matching(
            enclosing_counts, lambda repeat_count: repeat_count != reference_middle_repeat_count
        )
        flanking_nonreference = sum_matching(
            flanking_counts, lambda repeat_count: repeat_count != reference_middle_repeat_count
        )
        top_enclosing_repeat, top_enclosing_reads = max_supported_repeat(
            enclosing_counts, lambda repeat_count: repeat_count != reference_middle_repeat_count
        )
        top_flanking_repeat, top_flanking_reads = max_supported_repeat(
            flanking_counts, lambda repeat_count: repeat_count != reference_middle_repeat_count
        )

        row = {
            "sample_id": sample_id,
            "group": group,
            "vcf_path": str(vcf_path),
            "vcf_status": record.get("vcf_status", ""),
            "chrom": record.get("chrom", ""),
            "pos": record.get("pos", ""),
            "filter": record.get("filter", ""),
            "GT": record.get("GT", ""),
            "DP": record.get("DP", ""),
            "Q": record.get("Q", ""),
            "REPCN": record.get("REPCN", ""),
            "REPCI": record.get("REPCI", ""),
            "RC": record.get("RC", ""),
            "ENCLREADS": record.get("ENCLREADS", ""),
            "FLNKREADS": record.get("FLNKREADS", ""),
            "RU": record.get("info_RU", ""),
            "REF": record.get("info_REF", ""),
            "enclosing_counts": encode_histogram(enclosing_counts),
            "enclosing_total_reads": str(enclosing_total),
            "enclosing_nonreference_reads": str(enclosing_nonreference),
            "enclosing_nonreference_fraction": (
                f"{(enclosing_nonreference / enclosing_total):.8f}" if enclosing_total else ""
            ),
            "top_enclosing_nonreference_repeat_count": top_enclosing_repeat,
            "top_enclosing_nonreference_read_count": top_enclosing_reads,
            "flanking_counts": encode_histogram(flanking_counts),
            "flanking_total_reads": str(flanking_total),
            "flanking_nonreference_reads": str(flanking_nonreference),
            "flanking_nonreference_fraction": (
                f"{(flanking_nonreference / flanking_total):.8f}" if flanking_total else ""
            ),
            "top_flanking_nonreference_repeat_count": top_flanking_repeat,
            "top_flanking_nonreference_read_count": top_flanking_reads,
        }
        row.update(classified)
        if row["gangstr_review_required"] == "no" and (
            enclosing_nonreference > 0 or flanking_nonreference > 0
        ):
            row["gangstr_review_required"] = "yes"
        output_rows.append(row)

    fieldnames = [
        "sample_id",
        "group",
        "vcf_path",
        "vcf_status",
        "chrom",
        "pos",
        "filter",
        "GT",
        "DP",
        "Q",
        "REPCN",
        "REPCI",
        "RC",
        "ENCLREADS",
        "FLNKREADS",
        "RU",
        "REF",
        "allele1_middle_repeat_count",
        "allele2_middle_repeat_count",
        "allele1_middle_repeat_ci",
        "allele2_middle_repeat_ci",
        "allele1_total_repeat_count",
        "allele2_total_repeat_count",
        "allele1_delta_vs_reference",
        "allele2_delta_vs_reference",
        "middle_repeat_counts",
        "total_repeat_counts",
        "delta_vs_reference",
        "gangstr_interpretation",
        "enclosing_counts",
        "enclosing_total_reads",
        "enclosing_nonreference_reads",
        "enclosing_nonreference_fraction",
        "top_enclosing_nonreference_repeat_count",
        "top_enclosing_nonreference_read_count",
        "flanking_counts",
        "flanking_total_reads",
        "flanking_nonreference_reads",
        "flanking_nonreference_fraction",
        "top_flanking_nonreference_repeat_count",
        "top_flanking_nonreference_read_count",
        "gangstr_review_required",
    ]
    write_tsv(args.results_root / "gangstr_calls.tsv", output_rows, fieldnames)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
