#!/usr/bin/env python3
"""Derive the empirical LoD from markers that pass the complete pipeline."""

import csv
import subprocess
from pathlib import Path

# Use the direct counts from stage 7 and the final calls from stage 9.
ROOT = Path(__file__).resolve().parents[2]
RECOVERY = ROOT / "results2/spikein/read_recovery/mixture_read_recovery.tsv"
FILTERED = ROOT / "results2/spikein/filtermutectcalls/filtered"
SAMPLES = ROOT / "src2/spikein/samples.tsv"
REFERENCE = ROOT / "resources/references/snv/chr2_chr4_chr20.fasta"
OUTPUT = ROOT / "results2/spikein/read_recovery/empirical_lod.tsv"


def read_tsv(path):
    # Read a small project TSV as named fields.
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_calls(path, target_alleles):
    # Atomise compound records before matching exact marker alleles.
    result = subprocess.run(
        ["bcftools", "norm", "--atomize", "-m", "-any", "-f", REFERENCE, "-Ov", path],
        check=True, capture_output=True, text=True,
    )
    calls = {}
    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 10:
            raise ValueError(f"Expected one sample in {path}")
        chrom, position, _identifier, ref, alt, _qual, status, _info, form, sample = fields
        key = (chrom, position, ref, alt)
        if key not in target_alleles:
            continue
        if key in calls:
            raise ValueError(f"Duplicate atomised allele in {path}: {':'.join(key)}")
        sample_values = dict(zip(form.split(":"), sample.split(":"), strict=True))
        calls[key] = {"status": status, "mutect2_vaf": sample_values.get("AF", "NA")}
    return calls


def write_result(row):
    # Replace only the canonical LoD table after the result is complete.
    temporary = OUTPUT.with_suffix(".tsv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, list(row), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    temporary.replace(OUTPUT)


def main():
    # Require one read-level result for each marker and mixture.
    recovery = read_tsv(RECOVERY)
    targets = {}
    seen = set()
    for row in recovery:
        key = (row["marker_id"], row["sample_id"])
        if key in seen:
            raise ValueError(f"Duplicate recovery row: {', '.join(key)}")
        seen.add(key)
        allele = (row["chromosome"], row["position"], row["ref"], row["alt"])
        targets.setdefault(row["sample_id"], set()).add(allele)

    # Match the authoritative high and low sample roles to their filtered VCFs.
    samples = [row for row in read_tsv(SAMPLES) if row["role"] in ("high", "low")]
    samples.sort(key=lambda row: ("high", "low").index(row["role"]))
    if [row["role"] for row in samples] != ["high", "low"]:
        raise ValueError("Expected one high and one low mixture in samples.tsv")
    calls = {}
    for sample in samples:
        vcf = FILTERED / f"{sample['sample_id']}.filtered.vcf.gz"
        calls[sample["sample_id"]] = read_calls(vcf, targets[sample["sample_id"]])

    # Retain read-level recoveries that also pass the final exact-allele call.
    qualifying = []
    for row in recovery:
        allele = (row["chromosome"], row["position"], row["ref"], row["alt"])
        call = calls[row["sample_id"]].get(allele)
        if row["read_level_recovered"] == "True" and call and call["status"] == "PASS":
            depth = int(row["usable_depth"])
            count_vaf = int(row["alt_count"]) / depth
            if abs(count_vaf - float(row["alt_vaf"])) > 1e-12:
                raise ValueError(f"Count and VAF mismatch: {row['marker_id']}, {row['sample_id']}")
            qualifying.append((count_vaf, row, call))
    if not qualifying:
        raise ValueError("No read-level recovery passed Mutect2 and FilterMutectCalls")

    # The lowest direct-count VAF defines the empirical complete-pipeline LoD.
    vaf, source, call = min(qualifying, key=lambda result: result[0])
    result = {
        "estimate": "minimum_bam_readcount_vaf_passing_complete_pipeline",
        "scope": "all_markers_across_high_and_low_mixtures",
        "empirical_lod_vaf": f"{vaf:.10f}",
        "supporting_marker": source["marker_id"],
        "marker_category": source["category"],
        "supporting_mixture": source["sample_role"],
        "supporting_sample": source["sample_id"],
        "bam_readcount_alt_count": source["alt_count"],
        "bam_readcount_usable_depth": source["usable_depth"],
        "bam_readcount_fraction": f"{source['alt_count']}/{source['usable_depth']}",
        "bam_readcount_vaf": source["alt_vaf"],
        "mutect2_vaf": call["mutect2_vaf"],
        "filtermutectcalls_status": call["status"],
    }
    write_result(result)
    print(f"Empirical complete-pipeline LoD: {result['bam_readcount_fraction']} "
          f"({100 * vaf:.3f}%)")


if __name__ == "__main__":
    main()
