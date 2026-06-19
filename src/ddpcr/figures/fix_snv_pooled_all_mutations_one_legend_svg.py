#!/usr/bin/env python3
"""Normalise the final participant-pooled ddPCR SVG figure.

The R plotting script emits this SVG directly. This helper preserves the
manually adjusted bottom legend boxes, bold legend headings, and positioning,
then exports the matching manuscript PDF.
"""

from __future__ import annotations

import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path
from xml.etree import ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = PROJECT_ROOT / "manuscript" / "figures" / "ddpcr_fractional_abundance_pooled"
INPUT_SVG = OUT_DIR / "SNV_pooled_all_mutations_one_legend.svg"
OUTPUT_PDF = OUT_DIR / "SNV_pooled_all_mutations_one_legend.pdf"

SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)

LOD_GREEN = "#009E73"


def _ns_tag(tag: str) -> str:
    return f"{{{SVG_NS}}}{tag}"


def _iter_nodes_with_parent(root: ET.Element) -> Iterator[tuple[ET.Element, int, ET.Element]]:
    def walk(parent: ET.Element) -> Iterator[tuple[ET.Element, int, ET.Element]]:
        for index, child in enumerate(list(parent)):
            yield parent, index, child
            yield from walk(child)

    yield from walk(root)


def _find_by_id(root: ET.Element, node_id: str) -> ET.Element | None:
    for node in root.iter():
        if node.attrib.get("id") == node_id:
            return node
    return None


def _text_content(node: ET.Element) -> str:
    return "".join(node.itertext()).strip()


def _numeric_attrs(node: ET.Element, keys: tuple[str, ...]) -> list[float]:
    values = []
    for key in keys:
        value = node.attrib.get(key)
        if value is None:
            continue
        try:
            values.append(float(value))
        except ValueError:
            continue
    return values


def _node_in_box(node: ET.Element, x_min: float, x_max: float, y_min: float, y_max: float) -> bool:
    x_values = _numeric_attrs(node, ("x", "cx", "x1", "x2"))
    y_values = _numeric_attrs(node, ("y", "cy", "y1", "y2"))
    return (
        bool(x_values)
        and bool(y_values)
        and any(x_min <= value <= x_max for value in x_values)
        and any(y_min <= value <= y_max for value in y_values)
    )


def _set_attributes(node: ET.Element, attrs: dict[str, str]) -> None:
    for key, value in attrs.items():
        node.set(key, value)


def _group_nodes_by_box(
    root: ET.Element,
    group_id: str,
    attrs: dict[str, str],
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
) -> None:
    existing_group = _find_by_id(root, group_id)
    if existing_group is not None:
        _set_attributes(existing_group, attrs)
        return

    matches = [
        (parent, index, child)
        for parent, index, child in _iter_nodes_with_parent(root)
        if _node_in_box(child, x_min, x_max, y_min, y_max)
    ]
    if not matches:
        return

    parents = {id(parent) for parent, _index, _child in matches}
    if len(parents) != 1:
        raise RuntimeError(f"Cannot group legend nodes for {group_id}: nodes have different parents")

    parent = matches[0][0]
    ordered = sorted(matches, key=lambda item: item[1])
    group = ET.Element(_ns_tag("g"), {"id": group_id, **attrs})
    for _parent, _index, child in ordered:
        group.append(child)
    for _parent, index, _child in sorted(ordered, key=lambda item: item[1], reverse=True):
        del parent[index]
    parent.insert(ordered[0][1], group)


def _group_existing_nodes(root: ET.Element, group_id: str, child_ids: list[str], attrs: dict[str, str]) -> None:
    existing_group = _find_by_id(root, group_id)
    if existing_group is not None:
        _set_attributes(existing_group, attrs)
        return

    matches = [
        (parent, index, child)
        for parent, index, child in _iter_nodes_with_parent(root)
        if child.attrib.get("id") in child_ids
    ]
    if not matches:
        return

    parents = {id(parent) for parent, _index, _child in matches}
    if len(parents) != 1:
        raise RuntimeError(f"Cannot group legend nodes for {group_id}: nodes have different parents")

    parent = matches[0][0]
    ordered = sorted(matches, key=lambda item: item[1])
    group = ET.Element(_ns_tag("g"), {"id": group_id, **attrs})
    for _parent, _index, child in ordered:
        group.append(child)
    for _parent, index, _child in sorted(ordered, key=lambda item: item[1], reverse=True):
        del parent[index]
    parent.insert(ordered[0][1], group)


def _bold_legend_title(node: ET.Element | None) -> None:
    if node is None:
        return

    tspan_tag = _ns_tag("tspan")
    tspans = [child for child in list(node) if child.tag == tspan_tag]
    if tspans:
        tspan = tspans[0]
    else:
        text = "".join(node.itertext())
        node.text = None
        for child in list(node):
            node.remove(child)
        tspan = ET.SubElement(node, tspan_tag)
        tspan.text = text

    tspan.set(
        "style",
        "font-style:normal;font-variant:normal;font-weight:bold;"
        "font-stretch:normal;font-family:Arial;"
        "-inkscape-font-specification:'Arial Bold'",
    )


def _bold_legend_title_by_text(root: ET.Element, label: str) -> None:
    for node in root.iter(_ns_tag("text")):
        if _text_content(node) == label:
            _bold_legend_title(node)


def _set_rect_style_by_position(
    root: ET.Element,
    x: float,
    y: float,
    attrs: dict[str, str],
) -> None:
    for node in root.iter(_ns_tag("rect")):
        try:
            node_x = float(node.attrib.get("x", "nan"))
            node_y = float(node.attrib.get("y", "nan"))
        except ValueError:
            continue
        if abs(node_x - x) < 0.02 and abs(node_y - y) < 0.02:
            _set_attributes(node, attrs)


def _replace_lod_blue_with_green(root: ET.Element) -> None:
    for node in root.iter():
        for key, value in list(node.attrib.items()):
            node.set(key, value.replace("#0072B2", LOD_GREEN).replace("#0072b2", LOD_GREEN.lower()))


def _normalise_bottom_legends(root: ET.Element) -> None:
    above_lob_ids = [
        "rect1516", "text1518", "rect1520", "circle1522", "rect1524",
        "circle1526", "text1528", "text1530",
    ]
    above_lod_ids = [
        "rect1532", "text1534", "rect1536", "circle1538", "rect1540",
        "circle1542", "text1544", "text1546",
    ]

    _group_existing_nodes(root, "g2530", above_lob_ids, {"transform": "translate(17.953534)"})
    _group_existing_nodes(root, "g2519", above_lod_ids, {"transform": "translate(25.746463)"})
    _group_nodes_by_box(root, "g2530", {"transform": "translate(17.953534)"}, 160.0, 307.0, 818.0, 846.0)
    _group_nodes_by_box(root, "g2519", {"transform": "translate(25.746463)"}, 340.0, 488.0, 818.0, 846.0)

    for node_id, attrs in {
        "rect1516": {
            "x": "160.33",
            "y": "818.38",
            "width": "146.42999",
            "height": "27.24",
            "style": "fill:#ffffff;stroke:#000000;stroke-width:0.97;stroke-opacity:1",
        },
        "rect1532": {
            "x": "340.63",
            "y": "818.38",
            "width": "147.03",
            "height": "27.24",
            "style": "fill:#ffffff;stroke:#000000;stroke-width:0.97;stroke-opacity:1",
        },
    }.items():
        node = _find_by_id(root, node_id)
        if node is not None:
            _set_attributes(node, attrs)

    _set_rect_style_by_position(
        root,
        160.33,
        818.38,
        {
            "x": "160.33",
            "y": "818.38",
            "width": "146.42999",
            "height": "27.24",
            "style": "fill:#ffffff;stroke:#000000;stroke-width:0.97;stroke-opacity:1",
        },
    )
    _set_rect_style_by_position(
        root,
        340.63,
        818.38,
        {
            "x": "340.63",
            "y": "818.38",
            "width": "147.03",
            "height": "27.24",
            "style": "fill:#ffffff;stroke:#000000;stroke-width:0.97;stroke-opacity:1",
        },
    )
    _bold_legend_title(_find_by_id(root, "text1518"))
    _bold_legend_title(_find_by_id(root, "text1534"))
    _bold_legend_title_by_text(root, "Above LoB")
    _bold_legend_title_by_text(root, "Above LoD")


def _r_string(value: Path) -> str:
    return str(value).replace("\\", "/").replace('"', '\\"')


def _export_pdf_from_svg() -> None:
    expression = (
        "if (!requireNamespace('rsvg', quietly = TRUE)) "
        "stop('The rsvg package is required to export the fixed SVG to PDF'); "
        f'rsvg::rsvg_pdf("{_r_string(INPUT_SVG)}", "{_r_string(OUTPUT_PDF)}")'
    )
    subprocess.run(
        ["Rscript", "-e", expression],
        cwd=PROJECT_ROOT,
        check=True,
    )


def main() -> int:
    if not INPUT_SVG.exists():
        raise FileNotFoundError(f"Expected R output SVG is missing: {INPUT_SVG}")

    root = ET.parse(INPUT_SVG).getroot()
    _replace_lod_blue_with_green(root)
    _normalise_bottom_legends(root)

    ET.ElementTree(root).write(
        INPUT_SVG,
        encoding="utf-8",
        xml_declaration=True,
    )
    _export_pdf_from_svg()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
