#!/usr/bin/env python3
"""Render the sequencing2 PRNP call table from the final CJD TSV."""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
from pathlib import Path

TABLE_DIR = Path(__file__).resolve().parent
REPO_ROOT = TABLE_DIR.parents[1]
DEFAULT_INPUT = (
    REPO_ROOT
    / "results2"
    / "sequencing2"
    / "results"
    / "mutect2_cjd_dilutions_with_pon"
    / "variant_qc"
    / "cjd"
    / "filtered_prnp_variants.tsv"
)
DEFAULT_SETTINGS = DEFAULT_INPUT.parent / "run_settings.tsv"
DEFAULT_LOD = REPO_ROOT / "results2" / "spikein" / "read_recovery" / "empirical_lod.tsv"
DEFAULT_OUTPUT = TABLE_DIR / "table_prnp_somatic_snv_summary.tex"
DEFAULT_PREVIEW = TABLE_DIR / "table_prnp_somatic_snv_summary_preview.tex"
EXPECTED_KEYS = {
    ("CJD2", "chr20", "4691920", "G", "A"),
    ("CJD23", "chr20", "4691920", "G", "A"),
    ("CJD23", "chr20", "4694249", "T", "C"),
}
EXPECTED_LOD = Decimal(26) / Decimal(3891)
EXPECTED_LOD_TEXT = "0.006682086867129272"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--settings", type=Path, default=DEFAULT_SETTINGS)
    parser.add_argument("--lod", type=Path, default=DEFAULT_LOD)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--preview", type=Path, default=DEFAULT_PREVIEW)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Input table has no data rows: {path}")
    return rows


def validate_inputs(
    rows: list[dict[str, str]],
    settings_rows: list[dict[str, str]],
    lod_rows: list[dict[str, str]],
) -> Decimal:
    required = {
        "sample",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "FILTER",
        "DP",
        "REF_count",
        "ALT_count",
        "SB_refF",
        "SB_refR",
        "SB_altF",
        "SB_altR",
        "AAF",
        "gene",
        "mutation_type",
        "population_frequency",
        "MEAN_BQ",
        "MEAN_MQ",
    }
    missing = sorted(required - set(rows[0]))
    if missing:
        raise ValueError(f"CJD TSV is missing columns: {', '.join(missing)}")

    keys = {
        (row["sample"], row["CHROM"], row["POS"], row["REF"], row["ALT"])
        for row in rows
    }
    if keys != EXPECTED_KEYS:
        raise ValueError(f"Unexpected sequencing2 CJD membership: {sorted(keys)!r}")

    settings = {row["key"]: row["value"] for row in settings_rows}
    if (
        settings.get("enable_aaf_filter") != "TRUE"
        or settings.get("aaf_filter_applied") != "TRUE"
    ):
        raise ValueError("The CJD AAF filter is not recorded as enabled and applied")
    if settings.get("aaf_threshold") != EXPECTED_LOD_TEXT:
        raise ValueError(
            f"Unexpected CJD AAF threshold: {settings.get('aaf_threshold')!r}"
        )

    if len(lod_rows) != 1:
        raise ValueError("The empirical LoD table must contain exactly one row")
    lod = lod_rows[0]
    if (
        lod.get("bam_readcount_fraction") != "26/3891"
        or lod.get("bam_readcount_vaf") != EXPECTED_LOD_TEXT
        or lod.get("filtermutectcalls_status") != "PASS"
    ):
        raise ValueError("The empirical LoD evidence does not match 26/3891 PASS")

    threshold = Decimal(settings["aaf_threshold"])
    if abs(threshold - EXPECTED_LOD) > Decimal("1e-18"):
        raise ValueError("The recorded threshold does not equal 26/3891")
    for row in rows:
        if row["FILTER"] != "PASS" or row["gene"] != "PRNP":
            raise ValueError(
                f"Unexpected final call classification: {row['sample']} {row['POS']}"
            )
        if Decimal(row["AAF"]) <= threshold:
            raise ValueError(
                f"Final call does not exceed the LoD: {row['sample']} {row['POS']}"
            )
    return threshold


def format_row(row: dict[str, str]) -> str:
    population = row["population_frequency"].strip() or "not described"
    cells = [
        row["sample"],
        f"Chr{row['CHROM'].removeprefix('chr')}",
        row["POS"],
        f"{float(row['AAF']) * 100:.2f}",
        row["REF"],
        row["ALT"],
        row["mutation_type"].lower(),
        population,
        str(int(float(row["DP"]))),
        str(int(float(row["REF_count"]))),
        str(int(float(row["ALT_count"]))),
        f"{float(row['MEAN_BQ']):.2f}",
        f"{float(row['MEAN_MQ']):.2f}",
        row["SB_refF"],
        row["SB_refR"],
        row["SB_altF"],
        row["SB_altR"],
    ]
    return "      " + " & ".join(cells) + r" \\"


def render_table(rows: list[dict[str, str]]) -> str:
    ordered = sorted(
        rows, key=lambda row: (-float(row["AAF"]), row["sample"], row["POS"])
    )
    rendered_rows = "\n".join(format_row(row) for row in ordered)
    return f"""% Generated by manuscript2/tables2/make_prnp_somatic_snv_summary_tex.py
\\begin{{table}}[t]
  \\centering
  \\resizebox{{\\textwidth}}{{!}}{{
    \\begin{{tabular}}{{lllcccccrrrrrrrrr}}
      \\toprule
      Sample & Chromosome & Position & VAF (\\%) & REF & ALT & Mutation type
      & Population frequency & Read depth & REF count & ALT count & Mean BQ & Mean MQ
      & REF fwd & REF rev & ALT fwd & ALT rev \\\\
      \\midrule
{rendered_rows}
      \\bottomrule
    \\end{{tabular}}
  }}
  \\vspace{{0.5em}}
  \\caption{{\\textbf{{PRNP intronic calls retained by sequencing2.}} Three sample-level calls,
  representing two genomic sites, passed all configured filters, including the empirical
  complete-pipeline LoD of 26/3891 (0.668\\%). VAFs ranged from 0.82\\% to 0.95\\%.
  Abbreviations: VAF, variant allele fraction; REF, reference allele; ALT, alternative allele;
  BQ, base quality; MQ, mapping quality; fwd, forward; rev, reverse.}}
  \\label{{tab:sequencing2_prnp_somatic_snv_summary}}
\\end{{table}}

\\clearpage
"""


def render_preview(table_name: str) -> str:
    return f"""\\documentclass{{article}}
\\usepackage[margin=1cm]{{geometry}}
\\usepackage{{booktabs}}
\\usepackage{{graphicx}}

\\begin{{document}}
\\input{{{table_name}}}
\\end{{document}}
"""


def main() -> int:
    args = parse_args()
    rows = read_tsv(args.input)
    validate_inputs(rows, read_tsv(args.settings), read_tsv(args.lod))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_table(rows), encoding="utf-8")
    args.preview.write_text(render_preview(args.output.name), encoding="utf-8")
    print(f"Wrote {args.output}")
    print(f"Wrote {args.preview}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
