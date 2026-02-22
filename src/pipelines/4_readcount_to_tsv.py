#!/usr/bin/env python3
import argparse
import csv
import glob
import os


HEADER = [
    "CHROM",
    "POS",
    "REF",
    "BASE",
    "COUNT",
    "MEAN_BQ",
    "MEAN_MQ",
    "MEAN_POS",
    "FWD",
    "REV",
    "FRAC",
    "MMLQ",
    "MEAN_READ_POS",
    "FWD_BQ",
    "REV_BQ",
    "FRD",
    "SB_RATIO",
]


def parse_bam_readcount(input_file: str, output_file: str) -> None:
    with open(input_file, "r", encoding="utf-8") as infile, open(
        output_file, "w", newline="", encoding="utf-8"
    ) as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)

        for row in reader:
            if not row or row[0].startswith("#"):
                continue

            chrom = row[0]
            pos = row[1]
            ref = row[2]

            for allele_field in row[4:]:
                if not allele_field:
                    continue
                parts = allele_field.split(":")
                base = parts[0]
                metrics = parts[1:]

                if len(metrics) < 13:
                    metrics += [""] * (13 - len(metrics))

                out = [
                    chrom,
                    pos,
                    ref,
                    base,
                    metrics[0],   # COUNT
                    metrics[1],   # MEAN_BQ
                    metrics[2],   # MEAN_MQ
                    metrics[3],   # MEAN_POS
                    metrics[4],   # FWD
                    metrics[5],   # REV
                    metrics[6],   # FRAC
                    metrics[7],   # MMLQ
                    metrics[8],   # MEAN_READ_POS
                    metrics[9],   # FWD_BQ
                    metrics[10],  # REV_BQ
                    metrics[11],  # FRD
                    metrics[12],  # SB_RATIO
                ]
                writer.writerow(out)


def process_directory(input_dir: str, output_dir: str) -> int:
    os.makedirs(output_dir, exist_ok=True)
    input_files = sorted(glob.glob(os.path.join(input_dir, "*.txt")))

    for input_file in input_files:
        sample_name = os.path.splitext(os.path.basename(input_file))[0]
        output_path = os.path.join(output_dir, f"{sample_name}_metrics.tsv")
        parse_bam_readcount(input_file, output_path)
        print(f"Processed {input_file} -> {output_path}")

    return len(input_files)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert bam-readcount output files into flat TSV metrics tables."
    )
    parser.add_argument("--input-dir", required=True, help="Directory containing *.txt readcount files")
    parser.add_argument("--output-dir", required=True, help="Directory to write *_metrics.tsv files")
    args = parser.parse_args()

    n_files = process_directory(args.input_dir, args.output_dir)
    print(f"Done. Processed {n_files} readcount file(s).")


if __name__ == "__main__":
    main()
