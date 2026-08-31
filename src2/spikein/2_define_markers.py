#!/usr/bin/env python3
"""Select provisional donor SNPs from the joint pure-source VCF."""

import argparse
import csv
import hashlib
from pathlib import Path
import shlex
import subprocess
import sys

# Agreed thresholds; AF is ALT AD divided by total AD.
MIN_DP = 100
MIN_GQ = 30
MIN_DONOR_AD = 10
HET_MIN_AF = 0.30
HET_MAX_AF = 0.70
HOM_ALT_MIN_AF = 0.90
MAX_WT_ALT = 2
MAX_WT_AF = 0.001

ROOT = Path(__file__).resolve().parents[2]
VCF = ROOT / "results2/spikein/discovery/joint.vcf.gz"
MANIFEST = ROOT / "src2/spikein/samples.tsv"
REFERENCE = ROOT / "resources/references/snv/chr2_chr4_chr20.fasta"
A117V_SOURCE = ROOT / "results/sequencing_qc/a117v_spike_in_lod_table.csv"
DISCOVERY_SETTINGS = ROOT / "results2/spikein/discovery/run_settings.tsv"

# Fixed output columns, including unavailable genotypes for an absent control record.
HEADER = [
    "marker_id", "chromosome", "position", "ref", "alt", "category", "dbsnp_id",
    "vcf_filter", "reference_match",
    "donor_gt", "donor_gq", "donor_dp", "donor_ad", "donor_ref_count", "donor_alt_count", "donor_af",
    "wt_gt", "wt_gq", "wt_dp", "wt_ad", "wt_ref_count", "wt_alt_count", "wt_af",
    "candidate_qc_pass", "qc_reasons", "validation_status",
]


def sample_fields(text):
    # bcftools supplies a fixed field order; preserve missing and multiallelic AD.
    gt, gq, dp, ad = text.split(":")
    ref_count = alt_count = af = None
    depths = ad.split(",")
    if len(depths) == 2 and "." not in depths:
        ref_count, alt_count = map(int, depths)
        if min(ref_count, alt_count) < 0:
            raise ValueError("Negative allele depth")
        if ref_count + alt_count > 0:
            af = alt_count / (ref_count + alt_count)
    return dict(gt=gt, gq=None if gq == "." else int(gq), dp=None if dp == "." else int(dp),
                ad=ad, ref_count=ref_count, alt_count=alt_count, af=af)


def genotype_qc(donor, wt):
    # Check both source genotypes, allowing phased and reversed heterozygotes.
    donor_gt = sorted(donor["gt"].replace("|", "/").split("/"))
    reasons = []
    if donor_gt not in (["0", "1"], ["1", "1"]):
        reasons.append("donor_gt")
    if wt["gt"].replace("|", "/") != "0/0":
        reasons.append("wt_gt")
    for role, values in (("donor", donor), ("wt", wt)):
        if values["dp"] is None or values["dp"] < MIN_DP:
            reasons.append(f"{role}_dp")
        if values["gq"] is None or values["gq"] < MIN_GQ:
            reasons.append(f"{role}_gq")
        if values["af"] is None:
            reasons.append(f"{role}_ad_unavailable")

    # Require donor balance consistent with GT and minimal WT ALT support.
    if donor["af"] is not None:
        if donor["alt_count"] < MIN_DONOR_AD:
            reasons.append("donor_alt_count")
        if donor_gt == ["0", "1"]:
            if donor["ref_count"] < MIN_DONOR_AD or not HET_MIN_AF <= donor["af"] <= HET_MAX_AF:
                reasons.append("donor_heterozygous_balance")
        elif donor_gt == ["1", "1"] and donor["af"] < HOM_ALT_MIN_AF:
            reasons.append("donor_homozygous_alt_balance")
    if wt["af"] is not None and (wt["alt_count"] > MAX_WT_ALT or wt["af"] > MAX_WT_AF):
        reasons.append("wt_alt_support")
    return reasons


def write_tsv(path, header, rows):
    # Follow the existing TSV convention: fixed columns, NA values and Unix newlines.
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(["NA" if row[field] is None else row[field] for field in header])


def main():
    # Keep inputs fixed and write only to a new results2 directory.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results2/spikein/discovery/candidates"))
    output = (ROOT / parser.parse_args().output_dir).resolve()
    if not output.is_relative_to(ROOT / "results2/spikein") or output.exists():
        raise ValueError("Use a new output directory below results2/spikein")

    # Check the full sample list before requesting donor-then-WT columns.
    with MANIFEST.open() as handle:
        sources = [row for row in csv.DictReader(handle, delimiter="\t")
                   if row["role"] in ("pure_donor", "pure_wt")]
    samples = {row["role"]: row["sample_id"] for row in sources}
    if len(sources) != 2 or set(samples) != {"pure_donor", "pure_wt"} or len(set(samples.values())) != 2:
        raise ValueError("Manifest must identify two distinct pure samples")
    vcf_samples = subprocess.check_output(["bcftools", "query", "-l", str(VCF)], text=True).splitlines()
    if len(vcf_samples) != 2 or set(vcf_samples) != set(samples.values()):
        raise ValueError("Joint VCF must contain exactly the pure donor and pure WT")

    # Read the control allele, not historical mixture counts.
    with A117V_SOURCE.open() as handle:
        controls = {(row["chromosome"], int(row["position"]), row["ref"], row["alt"], row["dbsnp_id"])
                    for row in csv.DictReader(handle) if row["variant"] == "A117V"}
    if len(controls) != 1:
        raise ValueError("Historical A117V table must define one allele")
    control_chrom, control_position, control_ref, control_alt, control_rsid = controls.pop()
    control_key = (control_chrom, control_position, control_ref, control_alt)

    # Check every original REF. Discard norm output; selection uses the unchanged VCF.
    if not Path(str(REFERENCE) + ".fai").is_file():
        raise ValueError("The existing reference FASTA index is required")
    ref_check = ["bcftools", "norm", "-c", "e", "-f", str(REFERENCE), "-Ou", str(VCF)]
    subprocess.run(ref_check, check=True, stdout=subprocess.DEVNULL)
    control_query = ["samtools", "faidx", str(REFERENCE), f"{control_chrom}:{control_position}-{control_position}"]
    control_base = subprocess.check_output(control_query, text=True).splitlines()[1]
    if control_base != control_ref:
        raise ValueError("A117V disagrees with the reference")

    # Stage 1 restricted this VCF to capture targets. Query it without splitting records.
    query = ["bcftools", "query", "-s", f"{samples['pure_donor']},{samples['pure_wt']}",
             "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER[\t%GT:%GQ:%DP:%AD]\n", str(VCF)]
    audit = []
    for line in subprocess.check_output(query, text=True).splitlines():
        chrom, position, identifier, ref, alt, vcf_filter, donor_text, wt_text = line.split("\t")
        position = int(position)
        donor, wt = sample_fields(donor_text), sample_fields(wt_text)
        is_snv = ref in "ACGT" and alt in "ACGT" and len(ref) == len(alt) == 1 and ref != alt
        is_control = (chrom, position, ref, alt) == control_key
        reasons = genotype_qc(donor, wt)
        if not is_snv:
            reasons.insert(0, "not_biallelic_snv")
        category = "other_snp" if is_snv else "excluded_variant"
        rsid = identifier if identifier != "." else None
        if is_control:
            category, rsid = "a117v_positive_control", control_rsid
        row = {
            "marker_id": f"{chrom}:{position}:{ref}>{alt}",
            "chromosome": chrom, "position": position, "ref": ref, "alt": alt,
            "category": category, "dbsnp_id": rsid, "vcf_filter": vcf_filter, "reference_match": True,
            "donor_gt": donor["gt"], "donor_gq": donor["gq"], "donor_dp": donor["dp"],
            "donor_ad": donor["ad"], "donor_ref_count": donor["ref_count"],
            "donor_alt_count": donor["alt_count"], "donor_af": donor["af"],
            "wt_gt": wt["gt"], "wt_gq": wt["gq"], "wt_dp": wt["dp"],
            "wt_ad": wt["ad"], "wt_ref_count": wt["ref_count"],
            "wt_alt_count": wt["alt_count"], "wt_af": wt["af"],
            "candidate_qc_pass": not reasons, "qc_reasons": ";".join(reasons) or "none",
            "validation_status": "pending_pure_read_counts" if not reasons or is_control else "excluded",
        }
        audit.append(row)

    # A117V must still be counted if its exact VCF record is absent or fails QC.
    if not any(row["category"] == "a117v_positive_control" for row in audit):
        row = dict.fromkeys(HEADER)
        row.update(marker_id=f"{control_chrom}:{control_position}:{control_ref}>{control_alt}",
                   chromosome=control_chrom, position=control_position, ref=control_ref, alt=control_alt,
                   category="a117v_positive_control", dbsnp_id=control_rsid, reference_match=True,
                   candidate_qc_pass=False, qc_reasons="exact_a117v_record_absent",
                   validation_status="pending_pure_read_counts")
        audit.append(row)

    # Sort the tables; the counting stage will read coordinates from the candidate table.
    audit.sort(key=lambda row: (row["chromosome"], row["position"], row["ref"], row["alt"]))
    candidates = [row for row in audit if row["candidate_qc_pass"] or row["category"] == "a117v_positive_control"]
    # Record selection settings; stage 1 already holds detailed input and tool provenance.
    settings = dict(donor_sample=samples["pure_donor"], wt_sample=samples["pure_wt"],
                    min_dp=MIN_DP, min_gq=MIN_GQ, min_donor_ad=MIN_DONOR_AD,
                    het_min_af=HET_MIN_AF, het_max_af=HET_MAX_AF, hom_alt_min_af=HOM_ALT_MIN_AF,
                    max_wt_alt=MAX_WT_ALT, max_wt_af=MAX_WT_AF, command=shlex.join([sys.executable, *sys.argv]),
                    git_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                    stage1_settings=str(DISCOVERY_SETTINGS.relative_to(ROOT)),
                    stage1_settings_sha256=hashlib.sha256(DISCOVERY_SETTINGS.read_bytes()).hexdigest())

    # Write only after the input checks and provenance collection have succeeded.
    output.mkdir(parents=True)
    write_tsv(output / "candidate_markers.tsv", HEADER, candidates)
    write_tsv(output / "candidate_qc.tsv", HEADER, audit)
    write_tsv(output / "run_settings.tsv", ["key", "value"],
              [{"key": key, "value": value} for key, value in settings.items()])
    other_count = sum(row["category"] == "other_snp" for row in candidates)
    print(f"Selected {other_count} other candidate SNPs plus A117V: {output}")
    print("Pure-sample read validation is pending; no final marker set was written.")


if __name__ == "__main__":
    main()
