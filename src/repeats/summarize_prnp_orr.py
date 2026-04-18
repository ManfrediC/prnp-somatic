#!/usr/bin/env python3
"""Summarize per-sample ExpansionHunter outputs for the PRNP ORR screen."""

from __future__ import annotations

# Standard-library only on purpose so this summarizer can run in the lean repeat
# workflow environment without extra runtime dependencies.
import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Argument parsing and simple file readers
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize per-sample ExpansionHunter PRNP ORR outputs."
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--reference-total-repeats", required=True, type=int)
    parser.add_argument("--variable-repeat-offset", required=True, type=int)
    parser.add_argument("--reviewer-enabled", choices=("0", "1"), default="1")
    return parser.parse_args()


def read_manifest(path: Path) -> list[dict[str, str]]:
    # The manifest emitted by the shell runner is the authoritative cohort list.
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# ---------------------------------------------------------------------------
# Small parsing helpers for ExpansionHunter VCF fields
# ---------------------------------------------------------------------------

def parse_info_field(field: str) -> dict[str, str]:
    # INFO fields are semicolon-delimited key=value pairs.
    info = {}
    for item in field.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
    return info


def split_allele_values(raw: str) -> list[str]:
    # ExpansionHunter emits per-allele values with '/', '|', or ',' separators
    # depending on the field; normalize all of them to a simple list.
    if raw in {".", ""}:
        return []
    return [part for part in re.split(r"[/|,]", raw) if part != ""]


def parse_repeat_counts(raw: str) -> list[int]:
    # Keep only integer allele counts and silently drop placeholders.
    counts = []
    for token in split_allele_values(raw):
        try:
            counts.append(int(token))
        except ValueError:
            pass
    return counts


def parse_confidence_intervals(raw: str) -> list[tuple[int, int] | None]:
    # Confidence intervals are stored as "start-end" per allele. Any malformed
    # token is treated as missing and will later contribute to an uncertain call.
    intervals: list[tuple[int, int] | None] = []
    for token in split_allele_values(raw):
        match = re.fullmatch(r"(-?\d+)-(-?\d+)", token)
        if not match:
            intervals.append(None)
            continue
        intervals.append((int(match.group(1)), int(match.group(2))))
    return intervals


# ---------------------------------------------------------------------------
# Per-sample file parsing
# ---------------------------------------------------------------------------

def parse_sample_vcf(vcf_path: Path) -> dict[str, str]:
    # A missing or malformed VCF is surfaced explicitly so the final summary
    # tables can distinguish failed_qc from true reference-like calls.
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
            # The PRNP repeat catalog contains a single locus, so the first data
            # line is the only record we expect to summarize for each sample.
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

    # Pull out the specific INFO and FORMAT fields that are useful for final
    # cohort review and manuscript-facing table exports.
    info = parse_info_field(record_data["info_raw"])
    record_data["info_REF"] = info.get("REF", "")
    record_data["info_RL"] = info.get("RL", "")
    record_data["info_VARID"] = info.get("VARID", "")
    record_data["info_REPID"] = info.get("REPID", "")
    record_data["info_RU"] = info.get("RU", "")
    record_data["info_END"] = info.get("END", "")

    keys = record_data["format"].split(":")
    values = record_data["sample_field"].split(":")
    fmt = dict(zip(keys, values))
    for key in ["GT", "SO", "REPCN", "REPCI", "ADSP", "ADFL", "ADIR", "LC"]:
        record_data[key] = fmt.get(key, "")
    return record_data


def parse_sample_json(json_path: Path) -> dict[str, str]:
    # The JSON is mainly used as a presence/provenance check here, not as the
    # primary source of genotype values.
    if not json_path.exists():
        return {"json_status": "missing"}
    try:
        with json_path.open(encoding="utf-8") as handle:
            data = json.load(handle)
    except json.JSONDecodeError:
        return {"json_status": "malformed"}

    result = {"json_status": "ok"}
    result["json_locus_count"] = str(len(data.get("LocusResults", {})))
    return result


# ---------------------------------------------------------------------------
# Interpretation helpers
# ---------------------------------------------------------------------------

def classify_sample(
    record: dict[str, str],
    reference_total_repeats: int,
    variable_repeat_offset: int,
) -> dict[str, str]:
    # ExpansionHunter models the variable middle repeat block. To recover the
    # biologically meaningful total PRNP ORR repeat count, add the fixed repeat
    # offset that represents the non-variable flanking repeat elements.
    repcn = parse_repeat_counts(record.get("REPCN", ""))
    repci = parse_confidence_intervals(record.get("REPCI", ""))
    total_counts = [count + variable_repeat_offset for count in repcn]

    filter_value = record.get("filter", "")
    if record.get("vcf_status") != "ok":
        interpretation = "failed_qc"
    elif len(repcn) != 2:
        interpretation = "failed_qc"
    else:
        # A narrow, exact confidence interval is the most conservative signal of
        # a stable reference-like call in this one-locus screen.
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

        # Classification is intentionally conservative: any non-PASS filter or
        # non-exact confidence interval becomes uncertain rather than reference.
        if filter_value not in {"PASS", "."}:
            interpretation = "uncertain"
        elif any(total > reference_total_repeats for total in total_counts):
            interpretation = "candidate_OPRI"
        elif any(total < reference_total_repeats for total in total_counts):
            interpretation = "candidate_OPRD"
        elif has_wide_ci:
            interpretation = "uncertain"
        else:
            interpretation = "reference"

    # These columns are designed to be directly written into the cohort TSVs.
    row = {
        "middle_repeat_counts": ",".join(str(value) for value in repcn),
        "total_repeat_counts": ",".join(str(value) for value in total_counts),
        "delta_vs_reference": ",".join(
            str(value - reference_total_repeats) for value in total_counts
        ),
        "interpretation": interpretation,
        "review_required": "yes"
        if interpretation != "reference"
        else "no",
    }

    # Emit per-allele fields as fixed columns so downstream filtering and manual
    # review do not need to repeatedly split list-like strings.
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


def reviewer_artifact(prefix: Path) -> str:
    # Prefer a human-viewable image when available, while still tolerating runs
    # where only one plot format was generated.
    candidates = sorted(prefix.parent.glob(f"{prefix.name}*"))
    for path in candidates:
        if path.suffix.lower() in {".svg", ".png", ".pdf"}:
            return str(path)
    return ""


def pick_reference_calibrators(
    rows: list[dict[str, str]],
) -> dict[str, set[str]]:
    # The plan calls for reviewing a fixed subset of reference-like samples from
    # each group to calibrate manual interpretation.
    selected: dict[str, set[str]] = defaultdict(set)
    by_group: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        if row["interpretation"] == "reference":
            by_group[row["group"]].append(row["sample_id"])
    for group, sample_ids in by_group.items():
        for sample_id in sorted(sample_ids)[:2]:
            selected[group].add(sample_id)
    return selected


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    # Always create the parent directory so the summarizer can run against a
    # freshly prepared live results tree.
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


# ---------------------------------------------------------------------------
# Main summarization flow
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()
    manifest_rows = read_manifest(args.manifest)
    raw_root = args.results_root / "raw" / "expansionhunter"
    review_root = args.results_root / "review" / "reviewer"

    # Build the one-row-per-sample primary call table by combining the manifest,
    # VCF-derived genotype fields, JSON presence, and repeat interpretation.
    sample_rows: list[dict[str, str]] = []
    for manifest_row in manifest_rows:
        sample_id = manifest_row["sample_id"]
        sample_group = manifest_row["group"]
        sample_raw_dir = raw_root / sample_id
        sample_review_dir = review_root / sample_id
        vcf_path = sample_raw_dir / f"{sample_id}.vcf"
        json_path = sample_raw_dir / f"{sample_id}.json"
        record = parse_sample_vcf(vcf_path)
        record.update(parse_sample_json(json_path))
        classified = classify_sample(
            record,
            reference_total_repeats=args.reference_total_repeats,
            variable_repeat_offset=args.variable_repeat_offset,
        )

        row = {
            "sample_id": sample_id,
            "group": sample_group,
            "bam_path": manifest_row["bam_path"],
            "vcf_path": str(vcf_path),
            "json_path": str(json_path),
            "vcf_status": record.get("vcf_status", ""),
            "json_status": record.get("json_status", ""),
            "chrom": record.get("chrom", ""),
            "pos": record.get("pos", ""),
            "filter": record.get("filter", ""),
            "GT": record.get("GT", ""),
            "SO": record.get("SO", ""),
            "REPCN": record.get("REPCN", ""),
            "REPCI": record.get("REPCI", ""),
            "LC": record.get("LC", ""),
            "ADSP": record.get("ADSP", ""),
            "ADFL": record.get("ADFL", ""),
            "ADIR": record.get("ADIR", ""),
            "RU": record.get("info_RU", ""),
            "VARID": record.get("info_VARID", ""),
        }
        row.update(classified)
        row["review_artifact"] = reviewer_artifact(
            sample_review_dir / f"{sample_id}.PRNP_ORR"
        )
        sample_rows.append(row)

    # Build the manual review queue. Candidates and uncertain calls are always
    # reviewed, and a small reference subset from each group is added for
    # calibration.
    calibrators = pick_reference_calibrators(sample_rows)
    review_rows: list[dict[str, str]] = []
    for row in sample_rows:
        review_reason = "none"
        review_required = row["review_required"]
        if row["interpretation"] in {"candidate_OPRI", "candidate_OPRD", "uncertain"}:
            review_reason = row["interpretation"]
            review_required = "yes"
        elif row["sample_id"] in calibrators.get(row["group"], set()):
            review_reason = "reference_calibration_subset"
            review_required = "yes"

        # Distinguish between "review was not requested" and "review artifact was
        # never generated". In the current workflow REViewer images are created
        # for all samples, but only a subset is queued for manual adjudication.
        if review_required == "yes":
            review_status = (
                "not_reviewed"
                if args.reviewer_enabled == "1" and row["review_artifact"]
                else "not_generated"
            )
        elif args.reviewer_enabled == "1" and row["review_artifact"]:
            review_status = "generated_not_requested"
        else:
            review_status = "not_generated"

        review_rows.append(
            {
                "sample_id": row["sample_id"],
                "group": row["group"],
                "review_required": review_required,
                "review_reason": review_reason,
                "review_artifact": row["review_artifact"],
                "review_status": review_status,
                "review_notes": "none",
            }
        )

    # Candidate table is intentionally narrower than sample_calls.tsv: it keeps
    # only samples still computationally consistent with OPRI or OPRD.
    candidate_rows = [
        row
        for row in sample_rows
        if row["interpretation"] in {"candidate_OPRI", "candidate_OPRD"}
    ]

    # Summarize interpretation counts by group for the top-level cohort table.
    summary_counter: Counter[tuple[str, str]] = Counter()
    for row in sample_rows:
        summary_counter[(row["group"], row["interpretation"])] += 1

    summary_rows: list[dict[str, str]] = []
    for group in sorted({row["group"] for row in sample_rows}):
        total = sum(
            count
            for (counter_group, _), count in summary_counter.items()
            if counter_group == group
        )
        for interpretation in [
            "reference",
            "candidate_OPRI",
            "candidate_OPRD",
            "uncertain",
            "failed_qc",
        ]:
            count = summary_counter.get((group, interpretation), 0)
            fraction = f"{(count / total):.6f}" if total else ""
            summary_rows.append(
                {
                    "group": group,
                    "interpretation": interpretation,
                    "count": str(count),
                    "total_group_samples": str(total),
                    "fraction": fraction,
                }
            )

    # Keep output schemas explicit and stable so downstream scripts and manual
    # review can depend on fixed column names.
    sample_fieldnames = [
        "sample_id",
        "group",
        "bam_path",
        "vcf_path",
        "json_path",
        "vcf_status",
        "json_status",
        "chrom",
        "pos",
        "filter",
        "GT",
        "SO",
        "REPCN",
        "REPCI",
        "LC",
        "ADSP",
        "ADFL",
        "ADIR",
        "RU",
        "VARID",
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
        "interpretation",
        "review_required",
        "review_artifact",
    ]
    review_fieldnames = [
        "sample_id",
        "group",
        "review_required",
        "review_reason",
        "review_artifact",
        "review_status",
        "review_notes",
    ]
    summary_fieldnames = [
        "group",
        "interpretation",
        "count",
        "total_group_samples",
        "fraction",
    ]

    # Write all live-run summary artifacts expected by the shell workflow.
    write_tsv(args.results_root / "sample_calls.tsv", sample_rows, sample_fieldnames)
    write_tsv(args.results_root / "sample_review.tsv", review_rows, review_fieldnames)
    write_tsv(args.results_root / "candidate_calls.tsv", candidate_rows, sample_fieldnames)
    write_tsv(args.results_root / "cohort_summary.tsv", summary_rows, summary_fieldnames)
    return 0


if __name__ == "__main__":
    # Preserve a conventional CLI entrypoint for direct invocation from the
    # shell workflow or ad hoc local reruns.
    raise SystemExit(main())
