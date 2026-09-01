#!/usr/bin/env python3
"""Calculate read-level recovery from the two mixture count files."""

import argparse
import csv
import hashlib
import subprocess
import sys
from pathlib import Path

# Use the marker set by stage 5 and the raw counts written by stage 6.
ROOT = Path(__file__).resolve().parents[2]
MARKERS = ROOT / "results2/spikein/markers/informative_markers.tsv"
MARKER_SETTINGS = ROOT / "results2/spikein/markers/run_settings.tsv"
SAMPLES = ROOT / "src2/spikein/samples.tsv"
READCOUNTS = ROOT / "results2/spikein/readcount_qc/mixtures/readcounts"

# Apply the agreed technical read-level recovery thresholds.
MIN_DEPTH = 100
MIN_ALT = 10
MIN_ALT_STRAND = 3
MIN_MEAN_BQ = 20
MIN_MEAN_MQ = 20


def read_tsv(path):
    # Read a small project TSV as named fields.
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_count_file(path):
    # Read the A/C/G/T counts and allele qualities at each requested site.
    sites = {}
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            # Separate the site fields from the per-allele measurements.
            chrom, position, ref, reported_depth, *fields = row
            alleles = {}
            for field in fields:
                # Retain canonical bases and verify their strand counts.
                base, count, mq, bq, _, forward, reverse, *_unused = field.split(":")
                if base.upper() not in "ACGT":
                    continue
                count, forward, reverse = int(count), int(forward), int(reverse)
                if forward + reverse != count:
                    raise ValueError(f"Invalid strand counts: {path}, {chrom}:{position}, {base}")
                alleles[base.upper()] = {
                    "count": count, "mean_bq": float(bq), "mean_mq": float(mq),
                    "forward": forward, "reverse": reverse,
                }

            # Reject duplicate sites or incomplete A/C/G/T measurements.
            site = (chrom, int(position))
            if site in sites or set(alleles) != set("ACGT"):
                raise ValueError(f"Invalid A/C/G/T records: {path}, {chrom}:{position}")
            sites[site] = {
                "ref": ref.upper(), "reported_depth": int(reported_depth),
                "alleles": alleles,
            }
    return sites


def evaluate(marker, sample, site):
    # Calculate counts, VAF and read-level recovery.
    ref, alt = marker["ref"].upper(), marker["alt"].upper()
    if site["ref"] != ref:
        raise ValueError(f"Reference mismatch: {marker['marker_id']}, {sample['sample_id']}")

    # Use summed A/C/G/T counts as the denominator for direct VAF.
    ref_data, alt_data = site["alleles"][ref], site["alleles"][alt]
    depth = sum(allele["count"] for allele in site["alleles"].values())
    vaf = alt_data["count"] / depth if depth else "NA"
    evaluable = depth >= MIN_DEPTH

    # Read-level recovery is unavailable below the depth threshold.
    reasons = []
    if alt_data["count"] < MIN_ALT: reasons.append("alt_count")
    if alt_data["forward"] < MIN_ALT_STRAND: reasons.append("alt_forward")
    if alt_data["reverse"] < MIN_ALT_STRAND: reasons.append("alt_reverse")
    if alt_data["count"] and alt_data["mean_bq"] < MIN_MEAN_BQ: reasons.append("alt_mean_bq")
    if alt_data["count"] and alt_data["mean_mq"] < MIN_MEAN_MQ: reasons.append("alt_mean_mq")
    recovered = not reasons if evaluable else "NA"
    recovery_reasons = ";".join(reasons) or "none"
    if not evaluable:
        recovery_reasons = "insufficient_depth"

    # Expected AF is available only when the fraction is established.
    if sample["donor_fraction"] == "NA":
        expected_af = "NA"
    else:
        fraction = float(sample["donor_fraction"])
        expected_af = fraction * float(marker["donor_af"]) + (1 - fraction) * float(marker["wt_af"])

    # Keep marker identity, direct evidence and statuses in one audit row.
    return {
        "marker_id": marker["marker_id"], "chromosome": marker["chromosome"],
        "position": marker["position"], "ref": marker["ref"], "alt": marker["alt"],
        "category": marker["category"], "sample_role": sample["role"],
        "sample_id": sample["sample_id"], "expected_af": expected_af,
        "reported_depth": site["reported_depth"], "usable_depth": depth,
        "ref_count": ref_data["count"], "alt_count": alt_data["count"], "alt_vaf": vaf,
        "alt_forward": alt_data["forward"], "alt_reverse": alt_data["reverse"],
        "alt_mean_bq": alt_data["mean_bq"] if alt_data["count"] else "NA",
        "alt_mean_mq": alt_data["mean_mq"] if alt_data["count"] else "NA",
        "evaluable_at_depth": evaluable, "alt_read_present": alt_data["count"] > 0,
        "read_level_recovered": recovered,
        "recovery_failure_reasons": recovery_reasons,
    }


def write_tsv(path, rows):
    # Write deterministic columns and refuse an existing file.
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def derive_empirical_lod(rows):
    # Use the lowest recovered marker VAF across both mixtures.
    recovered_results = []
    for result in rows:
        if result["read_level_recovered"] is True:
            recovered_results.append(result)
    if not recovered_results:
        raise ValueError("Cannot derive an empirical LoD without a recovered marker")
    source = min(recovered_results, key=lambda result: result["alt_vaf"])
    vaf = source["alt_vaf"]

    # Retain the observation that defines the empirical value.
    return {
        "estimate": "minimum_observed_recovered_vaf",
        "scope": "all_markers_across_high_and_low_mixtures",
        "empirical_lod_vaf": f"{vaf:.10f}",
        "supporting_marker": source["marker_id"],
        "supporting_mixture": source["sample_role"],
        "supporting_sample": source["sample_id"],
        "alt_count": source["alt_count"],
        "usable_depth": source["usable_depth"],
    }


def main():
    # New files go to a fresh directory below results2/spikein.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results2/spikein/read_recovery"))
    output = (ROOT / parser.parse_args().output_dir).resolve()
    if not output.is_relative_to(ROOT / "results2/spikein") or output.exists():
        raise SystemExit("Use a new output directory below results2/spikein")

    # Require the unchanged marker table frozen by stage 5.
    markers = read_tsv(MARKERS)
    settings = {row["key"]: row["value"] for row in read_tsv(MARKER_SETTINGS)}
    marker_hash = hashlib.sha256(MARKERS.read_bytes()).hexdigest()
    if marker_hash != settings["informative_markers_sha256"]:
        raise ValueError("Informative marker table differs from the set frozen by stage 5")
    marker_sites = {(row["chromosome"], int(row["position"])) for row in markers}
    if not markers or len(marker_sites) != len(markers):
        raise ValueError("Informative markers must be non-empty and have unique sites")

    # Require the high and low mixture roles from the manifest.
    samples = [row for row in read_tsv(SAMPLES) if row["role"] in ("high", "low")]
    samples.sort(key=lambda row: ("high", "low").index(row["role"]))
    if [row["role"] for row in samples] != ["high", "low"]:
        raise ValueError("Expected one high and one low mixture in samples.tsv")

    # Evaluate every marker in both mixtures.
    counts = {}
    for sample in samples:
        path = READCOUNTS / f"{sample['sample_id']}.txt"
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty read-count file: {path}")
        counts[sample["role"]] = read_count_file(path)
        if set(counts[sample["role"]]) != marker_sites:
            raise ValueError(f"Read-count sites differ from the frozen marker set: {sample['sample_id']}")

    rows = []
    for marker in markers:
        site = (marker["chromosome"], int(marker["position"]))
        for sample in samples:
            rows.append(evaluate(marker, sample, counts[sample["role"]][site]))

    # Derive the empirical LoD from all recovered marker observations.
    empirical_lod = derive_empirical_lod(rows)

    # Write mixture recovery, the empirical LoD and concise run settings.
    output.mkdir(parents=True)
    write_tsv(output / "mixture_read_recovery.tsv", rows)
    write_tsv(output / "empirical_lod.tsv", [empirical_lod])

    # Record the fixed recovery thresholds separately from the derived LoD.
    run_settings = {
        "command": " ".join(sys.argv),
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "informative_markers_sha256": marker_hash,
        "min_depth": MIN_DEPTH, "min_alt": MIN_ALT,
        "min_alt_strand": MIN_ALT_STRAND, "min_mean_bq": MIN_MEAN_BQ,
        "min_mean_mq": MIN_MEAN_MQ,
    }
    write_tsv(output / "run_settings.tsv",
              [{"key": key, "value": value} for key, value in run_settings.items()])
    print(f"Evaluated {len(markers)} markers in both mixtures: {output}")
    print(f"Empirical LoD: {empirical_lod['alt_count']}/{empirical_lod['usable_depth']} "
          f"({100 * float(empirical_lod['empirical_lod_vaf']):.3f}%)")


if __name__ == "__main__":
    main()
