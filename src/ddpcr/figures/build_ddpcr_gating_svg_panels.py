#!/usr/bin/env python3
"""Assemble ddPCR gating figure panels from individual SVG plots."""

from __future__ import annotations

import copy
import csv
import re
import shutil
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parents[3]
POSITIVE_DIR = PROJECT_ROOT / "manuscript" / "figures" / "ddpcr_gating_lob_lod_positive"
STRATEGY_DIR = PROJECT_ROOT / "manuscript" / "figures" / "ddpcr_gating_strategy"

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)

URL_REF_RE = re.compile(r"url\(#([^)]+)\)")

# Keep legend colours in sync with the R scripts that generate the individual
# gating SVGs.
CLASS_COLOURS = {
    "Reference-only": "#0072B2",
    "Mutant-only": "#D55E00",
    "Double-positive": "#CC79A7",
    "Double-negative": "#9CA3AF",
    "Gated/unassigned": "#E69F00",
    "Rejected/unassigned": "#E5E7EB",
}


def read_manifest(path: Path) -> list[dict[str, str]]:
    """Read the R-generated plot manifest that lists individual panel SVGs."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def text(
    parent: ET.Element,
    x: float,
    y: float,
    value: str,
    size: int = 18,
    weight: str = "normal",
    fill: str = "#111827",
) -> None:
    """Add plain SVG text using the common panel typography."""
    node = ET.SubElement(
        parent,
        f"{{{SVG_NS}}}text",
        {
            "x": f"{x:g}",
            "y": f"{y:g}",
            "font-family": "Arial, Helvetica, sans-serif",
            "font-size": str(size),
            "font-weight": weight,
            "fill": fill,
        },
    )
    node.text = value


def circle(parent: ET.Element, cx: float, cy: float, r: float, fill: str, stroke: str = "none") -> None:
    """Add a legend marker circle."""
    ET.SubElement(
        parent,
        f"{{{SVG_NS}}}circle",
        {
            "cx": f"{cx:g}",
            "cy": f"{cy:g}",
            "r": f"{r:g}",
            "fill": fill,
            "stroke": stroke,
        },
    )


def line(
    parent: ET.Element,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    stroke: str = "#111827",
    width: float = 1.4,
    dasharray: str | None = None,
) -> None:
    """Add a straight SVG line, optionally dashed for auto thresholds."""
    attrs = {
        "x1": f"{x1:g}",
        "y1": f"{y1:g}",
        "x2": f"{x2:g}",
        "y2": f"{y2:g}",
        "stroke": stroke,
        "stroke-width": f"{width:g}",
        "stroke-linecap": "round",
    }
    if dasharray:
        attrs["stroke-dasharray"] = dasharray
    ET.SubElement(parent, f"{{{SVG_NS}}}line", attrs)


def draw_strategy_legend(root: ET.Element, centre_x: float, y: float) -> None:
    """Draw the shared droplet-class and gate legend below strategy panels."""
    x = centre_x - 410
    text(root, x, y + 12, "Droplet class", size=12, weight="bold")
    class_items = list(CLASS_COLOURS.items())
    class_x = x + 110
    class_gap_x = 136
    class_gap_y = 22
    for index, (label, colour) in enumerate(class_items):
        col = index % 3
        row = index // 3
        item_x = class_x + col * class_gap_x
        item_y = y + 10 + row * class_gap_y
        circle(root, item_x, item_y - 4, 3.5, colour, stroke="#9CA3AF" if label == "Rejected/unassigned" else "none")
        text(root, item_x + 12, item_y, label, size=10)

    gate_x = class_x + 3 * class_gap_x + 38
    text(root, gate_x, y + 12, "Gate", size=12, weight="bold")
    line(root, gate_x + 52, y + 6, gate_x + 92, y + 6, width=1.2)
    text(root, gate_x + 102, y + 10, "Final QuantaSoft gate", size=10)
    line(root, gate_x + 52, y + 28, gate_x + 92, y + 28, width=1.2, dasharray="2 3")
    text(root, gate_x + 102, y + 32, "Auto threshold", size=10)


def svg_viewbox(root: ET.Element) -> tuple[float, float, float, float]:
    """Return an SVG viewBox, falling back to width and height attributes."""
    view_box = root.attrib.get("viewBox")
    if view_box:
        parts = [float(part) for part in re.split(r"[,\s]+", view_box.strip()) if part]
        if len(parts) == 4:
            return parts[0], parts[1], parts[2], parts[3]
    width = float(re.sub(r"[^0-9.]+$", "", root.attrib["width"]))
    height = float(re.sub(r"[^0-9.]+$", "", root.attrib["height"]))
    return 0.0, 0.0, width, height


def prefix_svg_ids(node: ET.Element, prefix: str) -> None:
    """Namespace embedded SVG IDs so repeated glyphs and clips do not collide."""
    for child in node.iter():
        if "id" in child.attrib:
            child.attrib["id"] = f"{prefix}{child.attrib['id']}"
        for key, value in list(child.attrib.items()):
            if key == "id":
                continue
            if isinstance(value, str):
                value = URL_REF_RE.sub(lambda match: f"url(#{prefix}{match.group(1)})", value)
                if value.startswith("#"):
                    value = f"#{prefix}{value[1:]}"
                child.attrib[key] = value


def inline_svg(
    parent: ET.Element,
    source_svg: Path,
    x: float,
    y: float,
    width: float,
    height: float,
    prefix: str,
) -> None:
    """Copy one source SVG into the output panel while preserving its aspect ratio."""
    source_root = ET.parse(source_svg).getroot()
    min_x, min_y, source_width, source_height = svg_viewbox(source_root)
    scale = min(width / source_width, height / source_height)
    tx = x + (width - source_width * scale) / 2 - min_x * scale
    ty = y + (height - source_height * scale) / 2 - min_y * scale

    wrapper = ET.SubElement(
        parent,
        f"{{{SVG_NS}}}g",
        {"transform": f"matrix({scale:g} 0 0 {scale:g} {tx:g} {ty:g})"},
    )
    for child in list(source_root):
        cloned = copy.deepcopy(child)
        prefix_svg_ids(cloned, prefix)
        wrapper.append(cloned)


def write_panel(
    rows: list[dict[str, str]],
    output_svg: Path,
    ncols: int,
    cell_width: int,
    cell_height: int,
    title: str,
    common_legend: str | None = None,
) -> Path:
    """Assemble a grid of individual SVG plots into one labelled panel."""
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    ncols = min(ncols, max(1, len(rows)))
    margin_x = 30
    margin_y = 56
    gap_x = 18
    gap_y = 28
    legend_height = 66 if common_legend == "strategy" else 0
    nrows = (len(rows) + ncols - 1) // ncols
    width = margin_x * 2 + ncols * cell_width + (ncols - 1) * gap_x
    height = margin_y + nrows * cell_height + (nrows - 1) * gap_y + legend_height + 24

    root = ET.Element(
        f"{{{SVG_NS}}}svg",
        {
            "width": str(width),
            "height": str(height),
            "viewBox": f"0 0 {width} {height}",
            "version": "1.1",
        },
    )
    text(root, margin_x, 30, title, size=22, weight="bold")

    # Lay out source SVGs in reading order and add panel letters in the margin.
    for index, row in enumerate(rows):
        grid_row, grid_col = divmod(index, ncols)
        x = margin_x + grid_col * (cell_width + gap_x)
        y = margin_y + grid_row * (cell_height + gap_y)
        letter = chr(ord("A") + index)
        inline_svg(root, Path(row["svg_path"]), x + 24, y, cell_width - 24, cell_height, f"p{index}_")
        text(root, x + 2, y + 24, letter, size=18, weight="bold")

    # Strategy panels use one shared legend to avoid repeating legends in every cell.
    if common_legend == "strategy":
        draw_strategy_legend(root, width / 2, height - legend_height + 10)

    ET.ElementTree(root).write(output_svg, encoding="utf-8", xml_declaration=True)
    return output_svg


def export_pdf(svg_path: Path) -> Path:
    """Export the assembled SVG panel to PDF using Inkscape."""
    pdf_path = svg_path.with_suffix(".pdf")
    inkscape = shutil.which("inkscape") or "/usr/bin/inkscape"
    if not Path(inkscape).exists():
        raise RuntimeError("Inkscape is required to export PDF panels")
    subprocess.run(
        [
            inkscape,
            svg_path.as_posix(),
            "--export-type=pdf",
            f"--export-filename={pdf_path.as_posix()}",
        ],
        check=True,
    )
    return pdf_path


def sort_positive(rows: list[dict[str, str]], plot_kind: str) -> list[dict[str, str]]:
    """Keep LoB+LoD+ panels in the R manifest order."""
    selected = [row for row in rows if row["plot_kind"] == plot_kind]
    return sorted(selected, key=lambda row: int(row["plot_order"]))


def sort_strategy(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Order control panels by mutation and gating stage."""
    assay_order = {"D178N": 0, "E200K": 1, "P102L": 2}
    stage_order = {"NTC": 0, "WT": 1, "Positive control": 2, "Final adjusted gate": 3}
    return sorted(rows, key=lambda row: (assay_order[row["assay"]], stage_order[row["stage"]]))


def main() -> None:
    """Build all combined gating panels from individual plot manifests."""
    # Load manifests written by the R figure-asset script.
    positive_rows = read_manifest(POSITIVE_DIR / "plot_manifest.csv")
    strategy_rows = read_manifest(STRATEGY_DIR / "plot_manifest.csv")
    positive_merged = sort_positive(positive_rows, "merged")
    positive_faceted = sort_positive(positive_rows, "faceted")

    outputs = []
    # Positive panels are only meaningful when both merged and faceted views exist.
    if positive_merged and positive_faceted:
        outputs.extend(
            [
                write_panel(
                    positive_merged,
                    POSITIVE_DIR / "ddpcr_lob_lod_positive_merged_panel.svg",
                    ncols=3,
                    cell_width=430,
                    cell_height=395,
                    title="LoB+LoD+ sample-region ddPCR gating: merged wells",
                ),
                write_panel(
                    positive_faceted,
                    POSITIVE_DIR / "ddpcr_lob_lod_positive_faceted_panel.svg",
                    ncols=3,
                    cell_width=430,
                    cell_height=395,
                    title="LoB+LoD+ sample-region ddPCR gating: contributing wells",
                ),
            ]
        )
    else:
        print("No LoB+LoD+ sample-region rows; skipping positive gating panels")

    # The gating-strategy panel is always written because it comes from control wells.
    outputs.append(
        write_panel(
            sort_strategy(strategy_rows),
            STRATEGY_DIR / "ddpcr_gating_strategy_panel.svg",
            ncols=4,
            cell_width=345,
            cell_height=315,
            title="ddPCR gating strategy from QuantaSoft JSON thresholds",
            common_legend="strategy",
        ),
    )

    # Export every assembled SVG to PDF and report both artefacts.
    for svg_path in outputs:
        pdf_path = export_pdf(svg_path)
        print(f"Wrote {svg_path}")
        print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    # Keep imports side-effect free; only build panels when executed as a script.
    main()
