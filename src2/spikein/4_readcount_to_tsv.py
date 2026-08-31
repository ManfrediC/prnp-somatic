#!/usr/bin/env python3
"""Convert pure-source bam-readcount output into per-allele TSV tables.
"""

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
READCOUNTS = ROOT / "results2/spikein/readcount_qc/pure/readcounts"
HEADER = [
    "CHROM", "POS", "REF", "DEPTH", "ACGT_DEPTH", "BASE", "COUNT",
    "MEAN_BQ", "MEAN_MQ", "FWD", "REV",
]


def parse_bam_readcount(input_file: Path, output_file: Path) -> None:
    # Retain one row per allele, including zero counts, N and indel events.
    with input_file.open(encoding="utf-8") as infile, output_file.open(
        "x", newline="", encoding="utf-8"
    ) as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)

        for row in csv.reader(infile, delimiter="\t"):
            chrom, position, ref, depth, *allele_fields = row
            position, depth = int(position), int(depth)
            alleles = [field.split(":") for field in allele_fields]

            # SNV VAF uses A/C/G/T counts; retain the reported depth separately.
            acgt_depth = sum(int(parts[1]) for parts in alleles if parts[0].upper() in ("A", "C", "G", "T"))
            for parts in alleles:
                base, count, mq, bq, _, fwd, rev, *unused = parts
                count, fwd, rev = int(count), int(fwd), int(rev)
                bq, mq = float(bq), float(mq)
                if min(count, fwd, rev) < 0 or fwd + rev != count:
                    raise ValueError(f"Invalid counts: {input_file}, {chrom}:{position}, {base}")

                # The raw order is MQ then BQ. Missing SM tags make the next
                # single-end quality field unusable, so it is not exported.
                # With no supporting reads, mean qualities are unavailable.
                writer.writerow([
                    chrom, position, ref.upper(), depth, acgt_depth, base.upper(), count,
                    bq if count else "NA", mq if count else "NA", fwd, rev,
                ])


def main() -> None:
    # Match stage 2's output option; never overwrite earlier results.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results2/spikein/readcount_qc/pure/metrics"))
    output = (ROOT / parser.parse_args().output_dir).resolve()
    if not output.is_relative_to(ROOT / "results2/spikein") or output.exists():
        raise SystemExit("Use a new output directory below results2/spikein")

    # Require both pure samples from the manifest. Mixtures are not read.
    with (ROOT / "src2/spikein/samples.tsv").open(encoding="utf-8") as manifest:
        samples = [row for row in csv.DictReader(manifest, delimiter="\t")
                   if row["role"] in ("pure_donor", "pure_wt")]
    if sorted(row["role"] for row in samples) != ["pure_donor", "pure_wt"]:
        raise SystemExit("Expected one pure_donor and one pure_wt in samples.tsv")
    inputs = sorted(READCOUNTS / f"{row['sample_id']}.txt" for row in samples)
    for input_file in inputs:
        if not input_file.is_file() or input_file.stat().st_size == 0:
            raise SystemExit(f"Missing or empty read-count file: {input_file}")

    # Preserve sample names and input site order in the converted tables.
    output.mkdir(parents=True)
    for input_file in inputs:
        output_file = output / f"{input_file.stem}_metrics.tsv"
        parse_bam_readcount(input_file, output_file)
        print(f"Processed {input_file.relative_to(ROOT)} -> {output_file.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
