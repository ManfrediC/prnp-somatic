#!/usr/bin/env python3
"""Normalise mutation labels in the final ddPCR SVG figure.

The R plotting script emits this SVG directly. This helper preserves the
manually adjusted label placement and Helvetica-bold typography for the three
mutation labels.
"""

from __future__ import annotations

import sys
import re
import subprocess
import tempfile
from collections.abc import Iterator
from pathlib import Path
from xml.etree import ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = PROJECT_ROOT / "manuscript" / "figures" / "ddpcr_fractional_abundance"
INPUT_SVG = OUT_DIR / "SNV_all_mutations_legend_bottom_final.svg"
OUTPUT_PDF = OUT_DIR / "SNV_all_mutations_legend_bottom_final.pdf"
PAGE_HEIGHT = 864.0

SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)

HELVETICA_LABEL_STYLE = (
    "font-variant:normal;"
    "font-weight:bold;"
    "font-size:30px;"
    "font-family:Helvetica;"
    "-inkscape-font-specification:Helvetica-Bold;"
    "writing-mode:lr-tb;"
    "fill:#000000;"
    "fill-opacity:1;"
    "fill-rule:nonzero;"
    "stroke:none"
)
TRANSLATE_RE = re.compile(
    r"^translate\(\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+))"
    r"(?:[\s,]+([+-]?(?:\d+(?:\.\d*)?|\.\d+)))?\s*\)$"
)

LABEL_SPEC = {
    "D178N": {
        "group": {
            "transform": "translate(0,-11.002067)",
            "id": "g1911",
        },
        "text": {
            "x": "663.41998",
            "y": "90.25",
            "style": HELVETICA_LABEL_STYLE,
            "id": "text1909",
            "tspan": {
                "x": "663.41998",
                "y": "90.25",
                "style": "font-size:30px",
                "id": "tspan1907",
            },
        },
    },
    "E200K": {
        "group": {
            "transform": "translate(5.7458586,-12.757078)",
            "id": "g3764",
        },
        "text": {
            "x": "656.56",
            "y": "330.35001",
            "style": HELVETICA_LABEL_STYLE,
            "id": "text3762",
            "tspan": {
                "x": "656.56",
                "y": "330.35001",
                "style": "font-size:30px",
                "id": "tspan3760",
            },
        },
    },
    "P102L": {
        "group": {
            "transform": "translate(5.1552824,-18.157041)",
            "id": "g5649",
        },
        "text": {
            "x": "661.34003",
            "y": "573.60999",
            "style": HELVETICA_LABEL_STYLE,
            "id": "text5647",
            "tspan": {
                "x": "661.34003",
                "y": "573.60999",
                "style": "font-size:30px",
                "id": "tspan5645",
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


def _normalise_bottom_legends(root: ET.Element) -> None:
    brain_region_ids = [
        "rect5904", "text5906", "rect5908", "circle5910", "line5912",
        "rect5914", "circle5916", "line5918", "rect5920", "circle5922",
        "line5924", "rect5926", "circle5928", "line5930", "rect5932",
        "circle5934", "line5936", "rect5938", "circle5940", "line5942",
        "rect5944", "circle5946", "line5948", "text5950", "text5952",
        "text5954", "text5956", "text5958", "text5960", "text5962",
    ]
    above_lob_ids = [
        "rect5964", "text5966", "rect5968", "circle5970", "rect5972",
        "circle5974", "text5976", "text5978",
    ]

    _group_existing_nodes(root, "g4014", brain_region_ids, {})
    _group_existing_nodes(root, "g3981", above_lob_ids, {"transform": "translate(0,11.13)"})
    _group_nodes_by_box(root, "g4014", {}, 121.0, 497.0, 807.0, 857.0)
    _group_nodes_by_box(root, "g3981", {"transform": "translate(0,11.13)"}, 524.0, 671.0, 807.0, 835.0)

    for node_id, attrs in {
        "rect5904": {
            "x": "121.23",
            "y": "807.25",
            "width": "375.22",
            "height": "49.5",
            "style": "fill:#ffffff;stroke:#000000;stroke-width:0.97;stroke-opacity:1",
        },
        "rect5964": {
            "x": "524.34003",
            "y": "807.25",
            "width": "146.42999",
            "height": "27.24",
            "style": "fill:#ffffff;stroke:#000000;stroke-width:0.97;stroke-opacity:1",
        },
    }.items():
        node = _find_by_id(root, node_id)
        if node is not None:
            _set_attributes(node, attrs)

    _set_rect_style_by_position(
        root,
        121.23,
        807.25,
        {
            "x": "121.23",
            "y": "807.25",
            "width": "375.22",
            "height": "49.5",
            "style": "fill:#ffffff;stroke:#000000;stroke-width:0.97;stroke-opacity:1",
        },
    )
    _set_rect_style_by_position(
        root,
        524.34,
        807.25,
        {
            "x": "524.34003",
            "y": "807.25",
            "width": "146.42999",
            "height": "27.24",
            "style": "fill:#ffffff;stroke:#000000;stroke-width:0.97;stroke-opacity:1",
        },
    )
    _bold_legend_title(_find_by_id(root, "text5906"))
    _bold_legend_title(_find_by_id(root, "text5966"))
    _bold_legend_title_by_text(root, "Brain region")
    _bold_legend_title_by_text(root, "Above LoB")


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
    group_attrs = {}
    if "transform" in spec["group"]:
        group_attrs["transform"] = spec["group"]["transform"]
    group_attrs["data-ddpcr-label"] = label
    group_attrs.update(
        {
            key: value
            for key, value in spec["group"].items()
            if key != "transform"
        }
    )
    group = ET.Element(_ns_tag("g"), group_attrs)
    text_attrs = {key: value for key, value in spec["text"].items() if key != "tspan"}
    text_node = ET.SubElement(group, _ns_tag("text"), text_attrs)
    tspan = ET.SubElement(text_node, _ns_tag("tspan"), spec["text"]["tspan"])
    tspan.text = label
    return group


def _r_string(value: Path) -> str:
    return str(value).replace("\\", "/").replace('"', '\\"')


def _write_label_free_svg(path: Path) -> None:
    """Write a temporary copy of the canonical SVG without mutation labels."""
    root = ET.parse(INPUT_SVG).getroot()
    replacements = list(_iter_label_nodes_with_parent(root))
    if len(replacements) != len(LABEL_SPEC):
        found = [label for _parent, _index, _node, label in replacements]
        raise RuntimeError(
            f"Expected 3 mutation labels in fixed SVG, found {len(found)} "
            f"({', '.join(sorted(found))})"
        )

    for parent, _index, node, _label in replacements:
        parent.remove(node)

    ET.ElementTree(root).write(
        path,
        encoding="utf-8",
        xml_declaration=True,
    )


def _export_base_pdf_from_svg(input_svg: Path, output_pdf: Path) -> None:
    expression = (
        "if (!requireNamespace('rsvg', quietly = TRUE)) "
        "stop('The rsvg package is required to export the fixed SVG to PDF'); "
        f'rsvg::rsvg_pdf("{_r_string(input_svg)}", "{_r_string(output_pdf)}")'
    )
    subprocess.run(
        ["Rscript", "-e", expression],
        cwd=PROJECT_ROOT,
        check=True,
    )


def _pdf_string(value: str) -> str:
    """Escape a short ASCII label for a PDF text string literal."""
    return value.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def _translated_svg_position(spec: dict[str, dict[str, str]]) -> tuple[float, float]:
    x = float(spec["text"]["x"])
    y = float(spec["text"]["y"])
    transform = spec["group"].get("transform")
    if transform is None:
        return x, y

    match = TRANSLATE_RE.match(transform)
    if match is None:
        raise RuntimeError(f"Unsupported label transform: {transform}")

    x += float(match.group(1))
    if match.group(2) is not None:
        y += float(match.group(2))
    return x, y


def _overlay_helvetica_bold_labels(base_pdf: Path) -> None:
    """Overlay mutation labels using the standard PDF Helvetica-Bold font."""
    try:
        import pikepdf
        from pikepdf import Dictionary, Name
    except ImportError as exc:
        raise RuntimeError("pikepdf is required to overlay Helvetica-Bold PDF labels") from exc

    with pikepdf.Pdf.open(base_pdf) as pdf:
        page = pdf.pages[0]
        contents = page.obj[Name("/Contents")]
        if isinstance(contents, pikepdf.Array):
            base_content = b"\n".join(stream.read_bytes() for stream in contents)
        else:
            base_content = contents.read_bytes()
        page.obj[Name("/Contents")] = pdf.make_stream(b"q\n" + base_content + b"\nQ\n")

        font_name = page.add_resource(
            Dictionary(
                Type=Name("/Font"),
                Subtype=Name("/Type1"),
                BaseFont=Name("/Helvetica-Bold"),
                Encoding=Name("/WinAnsiEncoding"),
            ),
            Name("/Font"),
            Name("/FHelveticaBold"),
            replace_existing=True,
        )

        commands = ["q", "BT", f"{font_name} 30 Tf", "0 0 0 rg"]
        for label, spec in LABEL_SPEC.items():
            svg_x, svg_y = _translated_svg_position(spec)
            x = svg_x
            y = PAGE_HEIGHT - svg_y
            commands.append(f"1 0 0 1 {x:.2f} {y:.2f} Tm")
            commands.append(f"({_pdf_string(label)}) Tj")
        commands.extend(["ET", "Q"])
        page.contents_add(("\n".join(commands) + "\n").encode("ascii"))
        pdf.save(OUTPUT_PDF)


def _export_pdf_from_svg() -> None:
    with tempfile.TemporaryDirectory(prefix="ddpcr_label_fix_") as tmp_dir:
        tmp_svg = Path(tmp_dir) / "label_free.svg"
        tmp_pdf = Path(tmp_dir) / "label_free.pdf"
        _write_label_free_svg(tmp_svg)
        _export_base_pdf_from_svg(tmp_svg, tmp_pdf)
        _overlay_helvetica_bold_labels(tmp_pdf)


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

    for parent, index, text_node, label in replacements:
        parent.remove(text_node)
        parent.insert(index, _replacement_label_group(label))

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
