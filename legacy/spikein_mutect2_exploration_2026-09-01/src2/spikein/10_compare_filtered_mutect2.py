#!/usr/bin/env python3
"""Compare filtered calls from the four completed Mutect2 runs."""

import csv
import gzip
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
REF_FASTA = REPO_ROOT / "resources/references/snv/chr2_chr4_chr20.fasta"
MARKERS = REPO_ROOT / "results2/spikein/markers/informative_markers.tsv"
JOINT_VCF = REPO_ROOT / "results2/spikein/discovery/joint.vcf.gz"
SAMPLES = REPO_ROOT / "src2/spikein/samples.tsv"
RUNS = (
    "original_defaults",
    "no_alignment_start_cap",
    "max_population_af_1",
    "both_changes",
)


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_filtered_records(path):
    total = 0
    passed = 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            total += 1
            passed += line.split("\t", 7)[6] == "PASS"
    return total, passed


def atomised_alleles(path):
    command = [
        "bcftools",
        "norm",
        "-a",
        "-m-any",
        "-f",
        str(REF_FASTA),
        "-Ov",
        str(path),
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    alleles = {}
    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        key = (fields[0], int(fields[1]), fields[3], fields[4])
        filters = fields[6]
        # Treat repeated representations as one allele; any PASS record takes precedence.
        if key not in alleles or filters == "PASS":
            alleles[key] = filters
        elif alleles[key] != "PASS":
            reasons = set(alleles[key].split(";")) | set(filters.split(";"))
            alleles[key] = ";".join(sorted(reasons))
    return alleles


def write_comparison(output_root):
    markers = read_tsv(MARKERS)
    marker_keys = {
        (row["chromosome"], int(row["position"]), row["ref"], row["alt"])
        for row in markers
    }
    pure_source_keys = set(atomised_alleles(JOINT_VCF))
    samples = [
        (row["role"], row["sample_id"])
        for row in read_tsv(SAMPLES)
        if row["role"] in {"high", "low"}
    ]

    fields = [
        "run",
        "sample_role",
        "sample_id",
        "total_records",
        "pass_records",
        "filtered_out_records",
        "nonmarker_alleles",
        "nonmarker_with_pure_source_support",
        "nonmarker_without_pure_source_support",
        "pass_nonmarker_alleles",
        "pass_nonmarker_with_pure_source_support",
        "pass_nonmarker_without_pure_source_support",
        "marker_id",
        "marker_category",
        "marker_emitted",
        "marker_filter_status",
        "marker_filter_reasons",
    ]

    output = output_root / "comparison.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for run in RUNS:
            for role, sample in samples:
                filtered = output_root / "filtered" / run / f"{sample}.filtered.vcf.gz"
                total, passed = count_filtered_records(filtered)
                alleles = atomised_alleles(filtered)
                nonmarkers = set(alleles) - marker_keys
                pass_nonmarkers = {key for key in nonmarkers if alleles[key] == "PASS"}

                counts = {
                    "total_records": total,
                    "pass_records": passed,
                    "filtered_out_records": total - passed,
                    "nonmarker_alleles": len(nonmarkers),
                    "nonmarker_with_pure_source_support": len(nonmarkers & pure_source_keys),
                    "nonmarker_without_pure_source_support": len(nonmarkers - pure_source_keys),
                    "pass_nonmarker_alleles": len(pass_nonmarkers),
                    "pass_nonmarker_with_pure_source_support": len(pass_nonmarkers & pure_source_keys),
                    "pass_nonmarker_without_pure_source_support": len(pass_nonmarkers - pure_source_keys),
                }

                for marker in markers:
                    key = (
                        marker["chromosome"],
                        int(marker["position"]),
                        marker["ref"],
                        marker["alt"],
                    )
                    filter_reasons = alleles.get(key, "NA")
                    writer.writerow(
                        {
                            "run": run,
                            "sample_role": role,
                            "sample_id": sample,
                            **counts,
                            "marker_id": marker["marker_id"],
                            "marker_category": marker["category"],
                            "marker_emitted": key in alleles,
                            "marker_filter_status": (
                                "PASS"
                                if filter_reasons == "PASS"
                                else "filtered"
                                if key in alleles
                                else "NA"
                            ),
                            "marker_filter_reasons": filter_reasons,
                        }
                    )


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 10_compare_filtered_mutect2.py OUTPUT_ROOT")
    write_comparison(Path(sys.argv[1]).resolve())


if __name__ == "__main__":
    main()
