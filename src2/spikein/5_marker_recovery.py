#!/usr/bin/env python3
"""Finalise informative markers using the two pure-source count tables."""

import argparse
import csv
import hashlib
from pathlib import Path
import subprocess


# Use the latest provisional markers and validated pure-source metrics.
ROOT = Path(__file__).resolve().parents[2]
CANDIDATES = ROOT / "results2/spikein/discovery/candidates_wt_vaf_only/candidate_markers.tsv"
DONOR_METRICS = ROOT / "results2/spikein/readcount_qc/pure/metrics/A100_1to2_metrics.tsv"
WT_METRICS = ROOT / "results2/spikein/readcount_qc/pure/metrics/NA100_undil_metrics.tsv"

# Apply the thresholds.
MIN_DEPTH = 100
MIN_DONOR_ALT = 10
MIN_DONOR_ALT_STRAND = 3
MIN_MEAN_BQ = 20
MIN_MEAN_MQ = 20
HET_MIN_AF = 0.30
HET_MAX_AF = 0.70
HOM_ALT_MIN_AF = 0.90
MAX_WT_AF = 0.001


def read_dict_tsv(path):
    # Read a small project TSV as named fields.
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_metrics(path):
    # Index each allele by genomic site and reject duplicate records.
    metrics = {}
    for row in read_dict_tsv(path):
        key = (row["CHROM"], int(row["POS"]), row["BASE"].upper())
        if key in metrics:
            raise ValueError(f"Duplicate metric row: {key[0]}:{key[1]}:{key[2]}")
        metrics[key] = row
    return metrics


def get_site_metrics(marker, metrics, role):
    # Require complete A/C/G/T rows and matching reference bases.
    site = (marker["chromosome"], int(marker["position"]))
    ref, alt = marker["ref"].upper(), marker["alt"].upper()
    alleles = {base: metrics.get((*site, base)) for base in "ACGT"}
    if any(row is None for row in alleles.values()):
        raise ValueError(f"Missing {role} A/C/G/T metrics: {marker['marker_id']}")
    if {row["REF"].upper() for row in alleles.values()} != {ref}:
        raise ValueError(f"Mismatched {role} reference: {marker['marker_id']}")

    # Use stage 4 depth after checking agreement across the four rows.
    ref_row, alt_row = alleles[ref], alleles[alt]
    depth = int(ref_row["ACGT_DEPTH"])
    if {int(row["ACGT_DEPTH"]) for row in alleles.values()} != {depth}:
        raise ValueError(f"Inconsistent {role} A/C/G/T depth: {marker['marker_id']}")
    quality_row = alt_row if role == "donor" else ref_row
    if quality_row["MEAN_BQ"] == "NA" or quality_row["MEAN_MQ"] == "NA":
        raise ValueError(f"Missing required {role} quality: {marker['marker_id']}")

    # Return the direct counts and qualities needed for marker selection.
    alt_count = int(alt_row["COUNT"])
    return {
        "reported_depth": int(ref_row["DEPTH"]), "usable_depth": depth,
        "ref_count": int(ref_row["COUNT"]), "alt_count": alt_count,
        "af": alt_count / depth, "alt_forward": int(alt_row["FWD"]),
        "alt_reverse": int(alt_row["REV"]),
        "mean_bq": float(quality_row["MEAN_BQ"]),
        "mean_mq": float(quality_row["MEAN_MQ"]),
    }


def evaluate_marker(marker, donor, wt):
    # Apply direct depth, support, strand, quality and WT background checks.
    reasons = []
    if marker["candidate_qc_pass"] != "True": reasons.append("candidate_qc")
    if donor["usable_depth"] < MIN_DEPTH: reasons.append("donor_depth")
    if wt["usable_depth"] < MIN_DEPTH: reasons.append("wt_depth")
    if donor["alt_count"] < MIN_DONOR_ALT: reasons.append("donor_alt_count")
    if donor["alt_forward"] < MIN_DONOR_ALT_STRAND: reasons.append("donor_alt_forward")
    if donor["alt_reverse"] < MIN_DONOR_ALT_STRAND: reasons.append("donor_alt_reverse")
    if donor["mean_bq"] < MIN_MEAN_BQ: reasons.append("donor_alt_mean_bq")
    if donor["mean_mq"] < MIN_MEAN_MQ: reasons.append("donor_alt_mean_mq")
    if wt["mean_bq"] < MIN_MEAN_BQ: reasons.append("wt_ref_mean_bq")
    if wt["mean_mq"] < MIN_MEAN_MQ: reasons.append("wt_ref_mean_mq")
    if wt["af"] > MAX_WT_AF: reasons.append("wt_af")

    # Apply the donor-balance threshold for the stage-2 donor genotype.
    donor_gt = sorted(marker["donor_gt"].replace("|", "/").split("/"))
    if donor_gt == ["0", "1"] and not HET_MIN_AF <= donor["af"] <= HET_MAX_AF:
        reasons.append("donor_heterozygous_balance")
    elif donor_gt == ["1", "1"] and donor["af"] < HOM_ALT_MIN_AF:
        reasons.append("donor_homozygous_alt_balance")
    elif donor_gt not in (["0", "1"], ["1", "1"]):
        raise ValueError(f"Unsupported donor GT: {marker['marker_id']}")

    # Keep source genotypes and direct evidence together in one row.
    return {
        "marker_id": marker["marker_id"], "chromosome": marker["chromosome"],
        "position": marker["position"], "ref": marker["ref"], "alt": marker["alt"],
        "category": marker["category"], "dbsnp_id": marker["dbsnp_id"],
        "donor_gt": marker["donor_gt"], "donor_gq": marker["donor_gq"],
        "donor_dp": marker["donor_dp"], "donor_ad": marker["donor_ad"],
        "wt_gt": marker["wt_gt"], "wt_gq": marker["wt_gq"],
        "wt_dp": marker["wt_dp"], "wt_ad": marker["wt_ad"],
        "candidate_qc_pass": marker["candidate_qc_pass"],
        "candidate_qc_reasons": marker["qc_reasons"],
        "donor_reported_depth": donor["reported_depth"],
        "donor_usable_depth": donor["usable_depth"], "donor_ref_count": donor["ref_count"],
        "donor_alt_count": donor["alt_count"], "donor_af": donor["af"],
        "donor_alt_forward": donor["alt_forward"], "donor_alt_reverse": donor["alt_reverse"],
        "donor_alt_mean_bq": donor["mean_bq"], "donor_alt_mean_mq": donor["mean_mq"],
        "wt_reported_depth": wt["reported_depth"], "wt_usable_depth": wt["usable_depth"],
        "wt_ref_count": wt["ref_count"], "wt_alt_count": wt["alt_count"], "wt_af": wt["af"],
        "wt_ref_mean_bq": wt["mean_bq"], "wt_ref_mean_mq": wt["mean_mq"],
        "marker_qc_pass": not reasons, "exclusion_reasons": ";".join(reasons) or "none",
    }


def write_tsv(path, rows):
    # Preserve result-field order and refuse existing files.
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    # Confine new files to a fresh directory below results2/spikein.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results2/spikein/markers"))
    output = (ROOT / parser.parse_args().output_dir).resolve()
    if not output.is_relative_to(ROOT / "results2/spikein") or output.exists():
        raise SystemExit("Use a new output directory below results2/spikein")

    # Require one unique marker per ID and genomic site.
    candidates = read_dict_tsv(CANDIDATES)
    candidates.sort(key=lambda row: (row["chromosome"], int(row["position"]), row["ref"], row["alt"]))
    marker_ids = [row["marker_id"] for row in candidates]
    marker_sites = {(row["chromosome"], int(row["position"])) for row in candidates}
    if len(candidates) != len(set(marker_ids)) or len(candidates) != len(marker_sites):
        raise ValueError("Candidate markers must have unique IDs and sites")

    # Require both pure-source tables to cover exactly the provisional sites.
    donor_metrics, wt_metrics = read_metrics(DONOR_METRICS), read_metrics(WT_METRICS)
    if ({key[:2] for key in donor_metrics} != marker_sites or
            {key[:2] for key in wt_metrics} != marker_sites):
        raise ValueError("Metric tables must contain exactly the candidate sites")

    # Evaluate every marker without reading mixture evidence.
    results = []
    for marker in candidates:
        donor = get_site_metrics(marker, donor_metrics, "donor")
        wt = get_site_metrics(marker, wt_metrics, "WT")
        results.append(evaluate_marker(marker, donor, wt))

    # Require one heterozygous A117V control that passes direct-count QC.
    controls = [row for row in results if row["category"] == "a117v_positive_control"]
    if len(controls) != 1:
        raise ValueError("Expected exactly one A117V positive control")
    if controls[0]["donor_gt"] not in ("0/1", "1/0", "0|1", "1|0"):
        raise ValueError("A117V donor genotype is not heterozygous")
    if not controls[0]["marker_qc_pass"]:
        raise ValueError("A117V pure-source QC failed")
    informative = [row for row in results if row["marker_qc_pass"]]

    # Write the complete audit and its passing informative-marker subset.
    output.mkdir(parents=True)
    write_tsv(output / "marker_qc.tsv", results)
    write_tsv(output / "informative_markers.tsv", informative)
    informative_hash = hashlib.sha256((output / "informative_markers.tsv").read_bytes()).hexdigest()

    # Record thresholds and hashes needed to reproduce the fixed marker set.
    settings = {
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "min_depth": MIN_DEPTH, "min_donor_alt": MIN_DONOR_ALT,
        "min_donor_alt_strand": MIN_DONOR_ALT_STRAND, "min_mean_bq": MIN_MEAN_BQ,
        "min_mean_mq": MIN_MEAN_MQ, "het_min_af": HET_MIN_AF,
        "het_max_af": HET_MAX_AF, "hom_alt_min_af": HOM_ALT_MIN_AF,
        "max_wt_af": MAX_WT_AF,
        "candidate_markers_sha256": hashlib.sha256(CANDIDATES.read_bytes()).hexdigest(),
        "donor_metrics_sha256": hashlib.sha256(DONOR_METRICS.read_bytes()).hexdigest(),
        "wt_metrics_sha256": hashlib.sha256(WT_METRICS.read_bytes()).hexdigest(),
        "informative_markers_sha256": informative_hash,
    }
    write_tsv(output / "run_settings.tsv",
              [{"key": key, "value": value} for key, value in settings.items()])

    # Report other SNPs separately from the A117V positive control.
    other_count = sum(row["category"] == "other_snp" for row in informative)
    print(f"Finalised {other_count} other SNPs plus A117V: {output}")
    print("No mixture data were read and no recovery statuses were calculated.")


if __name__ == "__main__":
    main()
