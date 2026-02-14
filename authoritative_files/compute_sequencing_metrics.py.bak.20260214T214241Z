#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def run(cmd: List[str]) -> str:
    """Run a command and return stdout (strip trailing newline). Raise on non-zero exit."""
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {p.returncode}): {' '.join(cmd)}\nSTDERR:\n{p.stderr}"
        )
    return p.stdout.rstrip("\n")


def parse_hex_flag(s: str) -> int:
    s = s.strip()
    if s.lower().startswith("0x"):
        return int(s, 16)
    return int(s)


def bed_total_bases(bed_path: Path) -> int:
    total = 0
    with bed_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            chrom, start, end, *rest = line.split("\t")
            total += int(end) - int(start)
    return total


def load_output_columns(schema_path: Path) -> List[str]:
    """
    Supports either:
      A) a one-line header file containing tab-separated output column names, OR
      B) a schema table whose first column is 'column_name' (with a header row).
    """
    lines = [ln.rstrip("\n") for ln in schema_path.open()]
    lines = [ln for ln in lines if ln.strip()]

    if not lines:
        raise RuntimeError(f"Schema/header file is empty: {schema_path}")

    # Case A: single line and it looks like the output header
    if len(lines) == 1 and "\t" in lines[0] and "sample_id" in lines[0]:
        return lines[0].split("\t")

    # Case B: tabular schema; first row is schema header, first column is column_name
    reader = csv.reader(lines, delimiter="\t")
    header = next(reader)
    if not header or header[0] not in ("column_name", "name"):
        # Fallback: treat first non-empty line as header
        return lines[0].split("\t")

    cols = []
    for row in reader:
        if not row:
            continue
        cols.append(row[0])
    if not cols:
        raise RuntimeError(f"No columns parsed from schema file: {schema_path}")
    return cols


def discover_bam(sample_id: str, input_dir: Path) -> Tuple[Path, str]:
    """
    Tries short name first (<sample_id>.bam), then long GATK name
    (<sample_id>.bwa.picard.markedDup.recal.bam).
    Returns (bam_path, bam_style).
    """
    short = input_dir / f"{sample_id}.bam"
    if short.exists():
        return short, "short"

    long_ = input_dir / f"{sample_id}.bwa.picard.markedDup.recal.bam"
    if long_.exists():
        return long_, "long"

    raise FileNotFoundError(f"Could not find BAM for sample_id={sample_id} in {input_dir}")


def picard_metrics_path_for(sample_id: str, input_dir: Path, bam_style: str) -> Optional[Path]:
    if bam_style == "long":
        p = input_dir / f"{sample_id}.bwa.picard.markedDup.metrics"
        return p if p.exists() else None

    candidates = [
        input_dir / f"{sample_id}.picard.markedDup.metrics",
        input_dir / f"{sample_id}.markedDup.metrics",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def parse_picard_markdup_metrics(metrics_path: Path) -> Tuple[Optional[float], Optional[int]]:
    """
    Parses Picard MarkDuplicates metrics file.
    Returns (PERCENT_DUPLICATION, ESTIMATED_LIBRARY_SIZE) or (None, None) if not found.
    """
    pct_dup = None
    est_lib = None

    lines = metrics_path.read_text().splitlines()

    for i, ln in enumerate(lines):
        if ln.startswith("LIBRARY") and "PERCENT_DUPLICATION" in ln and "ESTIMATED_LIBRARY_SIZE" in ln:
            headers = ln.split("\t")
            if i + 1 < len(lines):
                vals = lines[i + 1].split("\t")
                if len(vals) >= len(headers):
                    d = dict(zip(headers, vals))
                    try:
                        raw = d.get("PERCENT_DUPLICATION")
                        pct_dup = float(raw) if raw not in (None, "") else None
                    except ValueError:
                        pct_dup = None
                    try:
                        raw = d.get("ESTIMATED_LIBRARY_SIZE")
                        if raw not in (None, "", "NA"):
                            est_lib = int(float(raw))
                    except ValueError:
                        est_lib = None
            break

    return pct_dup, est_lib


def samtools_view_count(bam: Path, extra_args: List[str]) -> int:
    out = run(["samtools", "view", "-c"] + extra_args + [str(bam)])
    return int(out)


def depth_stats_from_bed(
    bam: Path,
    bed: Path,
    bases_total: int,
    min_mapq: int,
    exclude_flags_hex: str,
    cov_thr_x: Tuple[int, int, int],
) -> Dict[str, object]:
    """
    Computes depth distribution across all bases in a BED.
    """
    excl_int = parse_hex_flag(exclude_flags_hex)

    cmd = [
        "samtools", "depth",
        "-a",
        "-Q", str(min_mapq),
        "-G", str(excl_int),
        "-b", str(bed),
        str(bam),
    ]
    depths: List[int] = []
    for ln in run(cmd).splitlines():
        if not ln:
            continue
        parts = ln.split("\t")
        if len(parts) < 3:
            continue
        depths.append(int(parts[2]))

    if bases_total > len(depths):
        depths.extend([0] * (bases_total - len(depths)))

    if bases_total == 0:
        return {
            "mean": "NA", "median": "NA", "p20": "NA", "fold80": "NA",
            "pct_ge": ["NA", "NA", "NA"], "bases_ge": [0, 0, 0],
        }

    depths.sort()
    total_depth = sum(depths)
    mean = total_depth / bases_total

    mid = bases_total // 2
    if bases_total % 2 == 1:
        median = float(depths[mid])
    else:
        median = (depths[mid - 1] + depths[mid]) / 2.0

    p20_index = int(0.20 * (bases_total - 1))
    p20 = float(depths[p20_index])

    fold80 = "NA"
    if mean > 0 and p20 > 0:
        fold80 = mean / p20

    bases_ge = []
    pct_ge = []
    for thr in cov_thr_x:
        n = sum(1 for d in depths if d >= thr)
        bases_ge.append(n)
        pct_ge.append((n / bases_total) * 100.0)

    return {
        "mean": mean,
        "median": median,
        "p20": p20,
        "fold80": fold80,
        "pct_ge": pct_ge,
        "bases_ge": bases_ge,
    }


def fmt_float(x: object, decimals: int = 3) -> str:
    if x is None or x == "NA":
        return "NA"
    return f"{float(x):.{decimals}f}"


def fmt_int(x: object) -> str:
    if x is None or x == "NA":
        return "NA"
    return str(int(x))


def parse_cov_thr(s: str) -> Tuple[int, int, int]:
    parts = [p.strip() for p in s.split(",") if p.strip()]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("cov-thr must be three comma-separated integers, e.g. 100,500,1250")
    try:
        a, b, c = (int(parts[0]), int(parts[1]), int(parts[2]))
    except ValueError as e:
        raise argparse.ArgumentTypeError("cov-thr values must be integers") from e
    return (a, b, c)


def build_argparser() -> argparse.ArgumentParser:
    script_dir = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(
        description="Compute per-sample sequencing QC metrics from BAMs listed in a manifest."
    )

    # Portable defaults: look next to the script
    p.add_argument("--manifest", type=Path, default=script_dir / "manifest.tsv",
                   help="Path to manifest.tsv (default: alongside script).")
    p.add_argument("--schema", type=Path, default=script_dir / "sequencing_metrics_per_sample.schema.tsv",
                   help="Schema/header file defining output columns (default: alongside script).")

    # Keep your current absolute paths as defaults, but allow overrides
    p.add_argument("--targets-bed", type=Path, default=Path("/home/mcarta/databases/targets.bed"),
                   help="Targets BED file.")
    p.add_argument("--coding-bed", type=Path, default=Path("/home/mcarta/databases/prnp_coding.bed"),
                   help="PRNP coding BED file.")
    p.add_argument("--reference-fasta", type=Path, default=Path("/home/mcarta/databases/chr2_chr4_chr20.fasta"),
                   help="Reference FASTA (reported in output; not otherwise used).")

    p.add_argument("--min-mapq", type=int, default=20, help="Minimum MAPQ for depth (-Q) and counts.")
    p.add_argument("--cov-thr", type=parse_cov_thr, default=(100, 500, 1250),
                   help="Coverage thresholds as three comma-separated integers (default: 100,500,1250).")

    p.add_argument("--exclude-flags-common-hex", type=str, default="0xB04",
                   help="Samtools flag mask for common exclusions, written to output as-is.")
    p.add_argument("--exclude-flags-unique-hex", type=str, default="0xF04",
                   help="Samtools flag mask for unique exclusions (typically includes duplicates).")

    p.add_argument("--limit", type=int, default=0,
                   help="Process only the first N manifest rows (0 = all). Useful for fast smoke tests.")

    return p


def main() -> None:
    args = build_argparser().parse_args()

    manifest: Path = args.manifest
    schema: Path = args.schema
    targets_bed: Path = args.targets_bed
    coding_bed: Path = args.coding_bed
    reference_fasta: Path = args.reference_fasta

    cov_thr_x: Tuple[int, int, int] = args.cov_thr
    min_mapq: int = args.min_mapq
    excl_common_hex: str = args.exclude_flags_common_hex
    excl_unique_hex: str = args.exclude_flags_unique_hex
    limit: int = args.limit

    out_cols = load_output_columns(schema)

    target_bases_total = bed_total_bases(targets_bed)
    coding_bases_total = bed_total_bases(coding_bed)

    sys.stdout.write("\t".join(out_cols) + "\n")

    with manifest.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"sample_id", "group", "batch", "input_dir"}
        if not required.issubset(reader.fieldnames or []):
            raise RuntimeError(f"Manifest missing required columns {sorted(required)}; found {reader.fieldnames}")

        for i, row in enumerate(reader, start=1):
            if limit and i > limit:
                break

            sample_id = row["sample_id"].strip()
            group = row["group"].strip()
            batch = row["batch"].strip()
            input_dir = Path(row["input_dir"].strip())

            bam_path, bam_style = discover_bam(sample_id, input_dir)

            reads_total = samtools_view_count(bam_path, [])
            reads_mapped = samtools_view_count(bam_path, ["-F", "0x4"])
            reads_primary_mapped = samtools_view_count(bam_path, ["-F", excl_common_hex])
            reads_qc_fail = samtools_view_count(bam_path, ["-f", "0x200"])
            reads_duplicate = samtools_view_count(bam_path, ["-f", "0x400", "-F", excl_common_hex])

            on_target_reads = samtools_view_count(
                bam_path,
                ["-F", excl_common_hex, "-L", str(targets_bed)]
            )

            on_target_fraction = "NA"
            on_target_percent = "NA"
            if reads_primary_mapped > 0:
                on_target_fraction = on_target_reads / reads_primary_mapped
                on_target_percent = on_target_fraction * 100.0

            targ = depth_stats_from_bed(
                bam=bam_path, bed=targets_bed, bases_total=target_bases_total,
                min_mapq=min_mapq, exclude_flags_hex=excl_unique_hex, cov_thr_x=cov_thr_x
            )
            cod = depth_stats_from_bed(
                bam=bam_path, bed=coding_bed, bases_total=coding_bases_total,
                min_mapq=min_mapq, exclude_flags_hex=excl_unique_hex, cov_thr_x=cov_thr_x
            )

            pic_path = picard_metrics_path_for(sample_id, input_dir, bam_style)
            pct_dup, est_lib = (None, None)
            if pic_path is not None:
                pct_dup, est_lib = parse_picard_markdup_metrics(pic_path)

            st = bam_path.stat()
            bam_mtime_iso = dt.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")

            out: Dict[str, str] = {
                "sample_id": sample_id,
                "group": group,
                "batch": batch,
                "input_dir": str(input_dir),
                "bam_path": str(bam_path),
                "bam_style": bam_style,
                "targets_bed": str(targets_bed),
                "coding_bed": str(coding_bed),
                "reference_fasta": str(reference_fasta),
                "target_bases_total": str(target_bases_total),
                "coding_bases_total": str(coding_bases_total),
                "cov_thr_x1": str(cov_thr_x[0]),
                "cov_thr_x2": str(cov_thr_x[1]),
                "cov_thr_x3": str(cov_thr_x[2]),
                "min_mapq": str(min_mapq),
                "exclude_flags_common_hex": excl_common_hex,
                "exclude_flags_unique_hex": excl_unique_hex,
                "reads_total": str(reads_total),
                "reads_mapped": str(reads_mapped),
                "reads_primary_mapped": str(reads_primary_mapped),
                "reads_duplicate": str(reads_duplicate),
                "reads_qc_fail": str(reads_qc_fail),
                "on_target_reads": str(on_target_reads),
                "on_target_fraction": fmt_float(on_target_fraction, 6) if on_target_fraction != "NA" else "NA",
                "on_target_percent": fmt_float(on_target_percent, 3) if on_target_percent != "NA" else "NA",
                "target_mean_depth": fmt_float(targ["mean"], 3),
                "target_median_depth": fmt_float(targ["median"], 3),
                "target_p20_depth": fmt_float(targ["p20"], 3),
                "target_fold80": fmt_float(targ["fold80"], 3) if targ["fold80"] != "NA" else "NA",
                "pct_bases_ge_x1": fmt_float(targ["pct_ge"][0], 3),
                "pct_bases_ge_x2": fmt_float(targ["pct_ge"][1], 3),
                "pct_bases_ge_x3": fmt_float(targ["pct_ge"][2], 3),
                "bases_covered_ge_x1": str(targ["bases_ge"][0]),
                "bases_covered_ge_x2": str(targ["bases_ge"][1]),
                "bases_covered_ge_x3": str(targ["bases_ge"][2]),
                "coding_mean_depth": fmt_float(cod["mean"], 3),
                "coding_median_depth": fmt_float(cod["median"], 3),
                "coding_p20_depth": fmt_float(cod["p20"], 3),
                "coding_fold80": fmt_float(cod["fold80"], 3) if cod["fold80"] != "NA" else "NA",
                "coding_pct_bases_ge_x1": fmt_float(cod["pct_ge"][0], 3),
                "coding_pct_bases_ge_x2": fmt_float(cod["pct_ge"][1], 3),
                "coding_pct_bases_ge_x3": fmt_float(cod["pct_ge"][2], 3),
                "coding_bases_covered_ge_x1": str(cod["bases_ge"][0]),
                "coding_bases_covered_ge_x2": str(cod["bases_ge"][1]),
                "coding_bases_covered_ge_x3": str(cod["bases_ge"][2]),
                "picard_metrics_path": str(pic_path) if pic_path is not None else "NA",
                "pct_duplication": fmt_float(pct_dup, 6) if pct_dup is not None else "NA",
                "estimated_library_size": fmt_int(est_lib) if est_lib is not None else "NA",
                "bam_mtime_iso": bam_mtime_iso,
                "bam_size_bytes": str(st.st_size),
                "pipeline_run_id": os.environ.get("PIPELINE_RUN_ID", "NA"),
            }

            sys.stdout.write("\t".join(out.get(c, "NA") for c in out_cols) + "\n")


if __name__ == "__main__":
    main()
