"""Convert the lean ddPCR positive-well TeX rows into a styled booktabs table.

Input : manuscript/tables/ddpcr_e200k_positive_well_results/ddpcr_e200k_positive_well_results_rows.tex
        One ``a & b & ... \\\\`` line per row; the line ``%%SECTION:pooled%%``
        marks the start of the pooled block. Written by
        src/ddpcr/figures/ddpcr_positive_well_fractional_abundance.R.
Output: manuscript/tables/supplement/table_ddpcr_e200k_positive_well_results.tex
        Complete table environment in the style of
        manuscript/tables/main/table_prnp_somatic_snv_summary.tex.

All presentation decisions (headers, column widths, caption, label) live here;
the R step emits data only.
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
ROWS_TEX = REPO_ROOT / "manuscript" / "tables" / "ddpcr_e200k_positive_well_results" / "ddpcr_e200k_positive_well_results_rows.tex"
OUT_TEX = REPO_ROOT / "manuscript" / "tables" / "supplement" / "table_ddpcr_e200k_positive_well_results.tex"

SECTION_MARKER = "%%SECTION:pooled%%"

HEADERS = [
    "Sample", "Brain region", "Replicate",
    "Accepted droplets", "Mutant-positive droplets", "FA (\\%)",
    "95\\% CI (\\%)", "LoB FA (\\%)", "Above LoB", "Above LoD",
]
COL_WIDTHS_CM = [1.0, 1.4, 1.2, 1.6, 1.9, 0.9, 1.9, 1.2, 1.1, 1.1]

CAPTION = (
    "\\textbf{Detailed ddPCR results for the E200K-positive samples CJD4 pons "
    "and CJD21 thalamus.} Well rows report individual ddPCR wells from "
    "independent DNA preparations; pooled rows combine all wells of one "
    "sample. FA: fractional abundance of the E200K allele; 95\\% CI: "
    "Poisson-based confidence interval; LoB FA: measurement-specific limit of "
    "blank on the FA scale; Above LoB: mutant-positive droplet count above "
    "the 95th percentile of expected blank noise; Above LoD: lower 95\\% CI "
    "bound above the E200K limit of detection (0.067\\%). Abbreviations: FA, "
    "fractional abundance; CI, confidence interval; LoB, limit of blank; "
    "LoD, limit of detection."
)
LABEL = "tab:ddpcr_e200k_positive_well_results"


def read_rows(path):
    """Split the lean rows file into (well_rows, pooled_rows) cell lists."""
    blocks = {"wells": [], "pooled": []}
    block = "wells"
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line == SECTION_MARKER:
            block = "pooled"
            continue
        if not line.endswith("\\\\"):
            raise ValueError(f"Row does not end with '\\\\': {line!r}")
        cells = [cell.strip() for cell in line[:-2].split(" & ")]
        if len(cells) != len(HEADERS):
            raise ValueError(f"Expected {len(HEADERS)} cells: {line!r}")
        blocks[block].append(cells)
    if not blocks["wells"] or not blocks["pooled"]:
        raise ValueError("Both a well block and a pooled block are required.")
    return blocks["wells"], blocks["pooled"]


def col_spec():
    """Booktabs column spec: centred p-columns, one per header."""
    return "".join(f">{{\\centering\\arraybackslash}}p{{{w}cm}}" for w in COL_WIDTHS_CM)


def render(well_rows, pooled_rows):
    """Assemble the complete table environment."""

    def data_row(cells):
        # \rule{0pt}{1.2em} keeps the row height of the reference tables.
        return "\\rule{0pt}{1.2em}" + " & ".join(cells) + " \\\\"

    header_row = " & ".join(f"\\textbf{{{h}}}" for h in HEADERS) + " \\\\"
    lines = [
        "\\begin{table}[htbp]",
        "\t\\centering",
        "\t{\\normalsize",
        "\t\t\\begin{adjustbox}{max width=\\textwidth}",
        f"\t\t\t\\begin{{tabular}}{{{col_spec()}}}",
        "\t\t\t\t\\toprule",
        f"\t\t\t\t{header_row}",
        "\t\t\t\t\\midrule",
        *[f"\t\t\t\t{data_row(r)}" for r in well_rows],
        "\t\t\t\t\\midrule",
        *[f"\t\t\t\t{data_row(r)}" for r in pooled_rows],
        "\t\t\t\t\\bottomrule",
        "\t\t\t\\end{tabular}",
        "\t\t\\end{adjustbox}",
        "\t}",
        "",
        "\t\\vspace{0.5em}",
        "",
        "\t\\begin{adjustbox}{max width=\\textwidth}",
        "\t\t\\begin{minipage}{\\linewidth}",
        f"\t\t\t\\caption{{{CAPTION}}}",
        f"\t\t\t\\label{{{LABEL}}}",
        "\t\t\\end{minipage}",
        "\t\\end{adjustbox}",
        "",
        "\\end{table}",
        "",
        "\\clearpage",
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    well_rows, pooled_rows = read_rows(ROWS_TEX)
    OUT_TEX.write_text(render(well_rows, pooled_rows), encoding="utf-8")
    print(f"Wrote {OUT_TEX}")
