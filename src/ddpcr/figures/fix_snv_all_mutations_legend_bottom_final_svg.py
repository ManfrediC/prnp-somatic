#!/usr/bin/env python3
"""Normalise mutation labels in the final ddPCR SVG figure.

The R plotting script emits this SVG directly. This helper enforces the
manually adjusted label placement/typography observed in
SNV_all_mutations_legend_bottom_final_tmp.svg for the three mutation labels.
"""

from __future__ import annotations

import sys
import subprocess
from collections.abc import Iterator
from pathlib import Path
from xml.etree import ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = PROJECT_ROOT / "manuscript" / "figures" / "ddpcr_fractional_abundance"
INPUT_SVG = OUT_DIR / "SNV_all_mutations_legend_bottom_final.svg"
OUTPUT_PDF = OUT_DIR / "SNV_all_mutations_legend_bottom_final.pdf"

SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)

LABEL_SPEC = {
    "D178N": {
        "group": {},
        "text": {
            "x": "770.00",
            "y": "90.25",
            "text-anchor": "end",
            "style": (
                "font-variant:normal;"
                "font-weight:bold;"
                "font-size:30px;"
                "font-family:Arial;"
                "-inkscape-font-specification:Arial-Bold;"
                "writing-mode:lr-tb;"
                "fill:#000000;"
                "fill-opacity:1;"
                "fill-rule:nonzero;"
                "stroke:none"
            ),
            "tspan": {
                "x": "770.00",
                "y": "90.25",
                "style": "font-size:30px",
            },
        },
    },
    "E200K": {
        "group": {},
        "text": {
            "x": "770.00",
            "y": "333.30",
            "text-anchor": "end",
            "style": (
                "font-variant:normal;"
                "font-weight:bold;"
                "font-size:30px;"
                "font-family:Arial;"
                "-inkscape-font-specification:Arial-Bold;"
                "writing-mode:lr-tb;"
                "fill:#000000;"
                "fill-opacity:1;"
                "fill-rule:nonzero;"
                "stroke:none"
            ),
            "tspan": {
                "x": "770.00",
                "y": "333.30",
                "style": "font-size:30px",
            },
        },
    },
    "P102L": {
        "group": {},
        "text": {
            "x": "770.00",
            "y": "576.59",
            "text-anchor": "end",
            "style": (
                "font-variant:normal;"
                "font-weight:bold;"
                "font-size:30px;"
                "font-family:Arial;"
                "-inkscape-font-specification:Arial-Bold;"
                "writing-mode:lr-tb;"
                "fill:#000000;"
                "fill-opacity:1;"
                "fill-rule:nonzero;"
                "stroke:none"
            ),
            "tspan": {
                "x": "770.00",
                "y": "576.59",
                "style": "font-size:30px",
            },
        },
    },
}


def _ns_tag(tag: str) -> str:
    return f"{{{SVG_NS}}}{tag}"


def _get_text_content(node: ET.Element) -> str:
    text_parts = "".join(node.itertext()).strip()
    return text_parts


def _set_attributes(node: ET.Element, attrs: dict[str, str], skip_none: bool = True) -> None:
    for key, value in attrs.items():
        if value is None and skip_none:
            continue
        node.set(key, value)


def _iter_label_nodes_with_parent(root: ET.Element) -> Iterator[tuple[ET.Element, int, ET.Element, str]]:
    text_tag = _ns_tag("text")

    def walk(parent: ET.Element) -> Iterator[tuple[ET.Element, int, ET.Element, str]]:
        for index, child in enumerate(list(parent)):
            label = child.attrib.get("data-ddpcr-label")
            if label in LABEL_SPEC:
                yield parent, index, child, label
                continue
            if child.tag == text_tag:
                label = _get_text_content(child)
                if label in LABEL_SPEC:
                    yield parent, index, child, label
                    continue
            yield from walk(child)

    yield from walk(root)


def _replacement_label_group(label: str) -> ET.Element:
    spec = LABEL_SPEC[label]
    group_attrs = {"data-ddpcr-label": label, **spec["group"]}
    group = ET.Element(_ns_tag("g"), group_attrs)
    text_attrs = {key: value for key, value in spec["text"].items() if key != "tspan"}
    text_node = ET.SubElement(group, _ns_tag("text"), text_attrs)
    tspan = ET.SubElement(text_node, _ns_tag("tspan"), spec["text"]["tspan"])
    tspan.text = label
    return group


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
        raise FileNotFoundError(
            f"Expected R output SVG is missing: {INPUT_SVG}"
        )

    root = ET.parse(INPUT_SVG).getroot()
    matched_nodes = []
    replacements = []

    for parent, index, node, label in _iter_label_nodes_with_parent(root):
        replacements.append((parent, index, node, label))
        matched_nodes.append(label)

    if len(matched_nodes) != len(LABEL_SPEC):
        missing = set(LABEL_SPEC) - set(matched_nodes)
        raise RuntimeError(
            f"Expected 3 mutation labels, but found {len(matched_nodes)} "
            f"({', '.join(sorted(matched_nodes))}). Missing: {', '.join(sorted(missing))}"
        )

    for parent, _index, text_node, label in replacements:
        parent.remove(text_node)
        parent.append(_replacement_label_group(label))

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
