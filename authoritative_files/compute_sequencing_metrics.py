#!/usr/bin/env python3
from __future__ import annotations

import csv
import datetime as dt
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---- authoritative inputs ----
MANIFEST = Path("/add/authoritative_files/manifest.tsv")
SCHEMA_HEADER_FILE = Path("/add/authoritative_files/sequencing_metrics_per_sample.schema.tsv")

TARGETS_BED = Path("/home/mcarta/databases/targets.bed")
CODING_BED = Path("/home/mcarta/databases/prnp_coding.bed")  # hg38 chr20:4699221-4699982 (0-based BED)
REFERENCE_FASTA = Path("/home/mcarta/databases/chr2_chr4_chr20.fasta")

# Thresholds (must match schema)
COV_THR_X = (100, 500, 1250)
MIN_MAPQ = 20

# Flag sets (hex strings are written to output as-is)
# common: unmapped + secondary + QC-fail + supplementary
EXCLUDE_FLAGS_COMMON_HEX = "0xB04"
# unique: common + duplicates
EXCLUDE_FLAGS_UNIQUE_HEX = "0xF04"


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
    Tries short name first (<sample_id>.bam), then long GATK name (<sample_id>.bwa.picard.markedDup.recal.bam).
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

    # for short BAMs, only use a metrics file if it happens to exist
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

    # Find the table that includes PERCENT_DUPLICATION
    # Typical Picard structure: header lines, then a row header line with fields.
    for i, ln in enumerate(lines):
        if ln.startswith("LIBRARY") and "PERCENT_DUPLICATION" in ln and "ESTIMATED_LIBRARY_SIZE" in ln:
            headers = ln.split("\t")
            if i + 1 < len(lines):
                vals = lines[i + 1].split("\t")
                if len(vals) >= len(headers):
                    d = dict(zip(headers, vals))
                    try:
                        pct_dup = float(d.get("PERCENT_DUPLICATION")) if d.get("PERCENT_DUPLICATION") not in (None, "") else None
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


def depth_stats_from_bed(bam: Path, bed: Path, bases_total: int, min_mapq: int, exclude_flags_hex: str) -> Dict[str, object]:
    """
    Computes depth distribution across all bases in a BED.
    """
    excl_int = parse_hex_flag(exclude_flags_hex)

    # NOTE: -Q = minimum mapping quality

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

    # If samtools doesn't emit all positions for some reason, pad with zeros (defensive).
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

    # median
    mid = bases_total // 2
    if bases_total % 2 == 1:
        median = float(depths[mid])
    else:
        median = (depths[mid - 1] + depths[mid]) / 2.0

    # 20th percentile depth (i.e. depth at which 80% of bases are >= that depth)
    p20_index = int(0.20 * (bases_total - 1))
    p20 = float(depths[p20_index])

    fold80 = "NA"
    if mean > 0 and p20 > 0:
        fold80 = mean / p20

    bases_ge = []
    pct_ge = []
    for thr in COV_THR_X:
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


def main() -> None:
    out_cols = load_output_columns(SCHEMA_HEADER_FILE)

    target_bases_total = bed_total_bases(TARGETS_BED)
    coding_bases_total = bed_total_bases(CODING_BED)

    # Emit header
    sys.stdout.write("\t".join(out_cols) + "\n")

    # Read manifest
    with MANIFEST.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"sample_id", "group", "batch", "input_dir"}
        if not required.issubset(reader.fieldnames or []):
            raise RuntimeError(f"Manifest missing required columns {sorted(required)}; found {reader.fieldnames}")

        for row in reader:
            sample_id = row["sample_id"].strip()
            group = row["group"].strip()
            batch = row["batch"].strip()
            input_dir = Path(row["input_dir"].strip())

            bam_path, bam_style = discover_bam(sample_id, input_dir)

            # read counts (samtools view uses -q for MAPQ; depth uses -Q)
            reads_total = samtools_view_count(bam_path, [])
            reads_mapped = samtools_view_count(bam_path, ["-F", "0x4"])
            reads_primary_mapped = samtools_view_count(bam_path, ["-F", EXCLUDE_FLAGS_COMMON_HEX])
            reads_qc_fail = samtools_view_count(bam_path, ["-f", "0x200"])
            reads_duplicate = samtools_view_count(bam_path, ["-f", "0x400", "-F", EXCLUDE_FLAGS_COMMON_HEX])

            # On-target reads (primary mapped; no duplicate removal here)
            on_target_reads = samtools_view_count(
                bam_path,
                ["-F", EXCLUDE_FLAGS_COMMON_HEX, "-L", str(TARGETS_BED)]
            )

            on_target_fraction = "NA"
            on_target_percent = "NA"
            if reads_primary_mapped > 0:
                on_target_fraction = on_target_reads / reads_primary_mapped
                on_target_percent = on_target_fraction * 100.0

            # Depth stats (exclude duplicates via EXCLUDE_FLAGS_UNIQUE_HEX)
            targ = depth_stats_from_bed(
                bam=bam_path, bed=TARGETS_BED, bases_total=target_bases_total,
                min_mapq=MIN_MAPQ, exclude_flags_hex=EXCLUDE_FLAGS_UNIQUE_HEX
            )
            cod = depth_stats_from_bed(
                bam=bam_path, bed=CODING_BED, bases_total=coding_bases_total,
                min_mapq=MIN_MAPQ, exclude_flags_hex=EXCLUDE_FLAGS_UNIQUE_HEX
            )

            # Picard metrics (if present)
            pic_path = picard_metrics_path_for(sample_id, input_dir, bam_style)
            pct_dup, est_lib = (None, None)
            if pic_path is not None:
                pct_dup, est_lib = parse_picard_markdup_metrics(pic_path)

            st = bam_path.stat()
            bam_mtime_iso = dt.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")

            # Build output row dict (strings only)
            out: Dict[str, str] = {
                "sample_id": sample_id,
                "group": group,
                "batch": batch,
                "input_dir": str(input_dir),
                "bam_path": str(bam_path),
                "bam_style": bam_style,
                "targets_bed": str(TARGETS_BED),
                "coding_bed": str(CODING_BED),
                "reference_fasta": str(REFERENCE_FASTA),
                "target_bases_total": str(target_bases_total),
                "coding_bases_total": str(coding_bases_total),
                "cov_thr_x1": str(COV_THR_X[0]),
                "cov_thr_x2": str(COV_THR_X[1]),
                "cov_thr_x3": str(COV_THR_X[2]),
                "min_mapq": str(MIN_MAPQ),
                "exclude_flags_common_hex": EXCLUDE_FLAGS_COMMON_HEX,
                "exclude_flags_unique_hex": EXCLUDE_FLAGS_UNIQUE_HEX,
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

            # Emit row in locked column order
            sys.stdout.write("\t".join(out.get(c, "NA") for c in out_cols) + "\n")


if __name__ == "__main__":
    main()
