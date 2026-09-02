#!/usr/bin/env python3
"""Apply the reviewed Figure 8 edits to the sequencing2 lollipop SVG.

The original reviewed SVG remains the layout specification. Its axes, labels,
colours, exon track and typography are retained. The original R head paths are
reused, the obsolete CJD6 elements are removed, and the retained heads, stems,
legend and LoD line are moved to the sequencing2 values.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
XML_NS = "http://www.w3.org/XML/1998/namespace"
INKSCAPE_NS = "http://www.inkscape.org/namespaces/inkscape"
SODIPODI_NS = "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)
ET.register_namespace("xml", XML_NS)
ET.register_namespace("inkscape", INKSCAPE_NS)
ET.register_namespace("sodipodi", SODIPODI_NS)

EXPECTED = [
    ("CJD2", 4691920, 0.8864499788940482, "G", "A"),
    ("CJD23", 4691920, 0.9463113171108536, "G", "A"),
    ("CJD23", 4694249, 0.8172362555720654, "T", "C"),
]
EXPECTED_LOD_PERCENT = 100 * 26 / 3891
OLD_LOD_PERCENT = 0.81
OLD_LOD_Y = 178.32896
Y_PER_VAF_PERCENT = -157.85986363636362
ANNOTATION_STYLE = (
    "font-variant:normal;font-weight:normal;font-family:Arial;"
    "-inkscape-font-specification:Arial;writing-mode:lr-tb;"
    "fill:#000000;fill-opacity:1;fill-rule:nonzero;stroke:none"
)

HEAD_LAYOUT = {
    "blue": {
        "head_id": "path3080",
        "solid_id": "path3076",
        "dotted_id": "path3078",
        "old_vaf": 0.55,
        "new_vaf": EXPECTED[0][2],
        "old_y": 209.90185,
    },
    "green_uncertain": {
        "head_id": "path3092",
        "solid_id": "path3088",
        "dotted_id": "path3090",
        "old_vaf": 0.60,
        "new_vaf": EXPECTED[1][2],
        "old_y": 202.00602,
    },
    "green_high": {
        "head_id": "path3098",
        "solid_id": "path3094",
        "dotted_id": "path3096",
        "old_vaf": 0.82,
        "new_vaf": EXPECTED[2][2],
        "old_y": 167.27685,
    },
}

OBSOLETE_IDS = {
    "path3082",
    "path3084",
    "path3086",
    "path1013-2",
    "path1019-0",
    "text871-7-8-0",
    "use3124",
    "use3126",
    "use3128",
    "use3130",
    "use3132",
    "use3134",
    "use3136",
    "use3140",
    "use3142",
    "use3144",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    required = {"sample", "position", "vaf", "ref", "alt"}
    if not rows or not required.issubset(rows[0]):
        missing = sorted(required.difference(rows[0] if rows else set()))
        raise ValueError(
            f"Lollipop data is missing required columns: {', '.join(missing)}"
        )
    return rows


def validate_rows(rows: list[dict[str, str]]) -> None:
    observed = [
        (row["sample"], int(row["position"]), float(row["vaf"]), row["ref"], row["alt"])
        for row in rows
    ]
    if len(observed) != len(EXPECTED):
        raise ValueError(
            f"Expected {len(EXPECTED)} lollipop rows, found {len(observed)}"
        )
    for actual, expected in zip(observed, EXPECTED, strict=True):
        same_identity = actual[:2] == expected[:2] and actual[3:] == expected[3:]
        if not same_identity or abs(actual[2] - expected[2]) > 1e-12:
            raise ValueError(
                f"R lollipop data differs from sequencing2: {actual!r} != {expected!r}"
            )


def local_name(tag: str) -> str:
    if not isinstance(tag, str):
        return ""
    return tag.rsplit("}", 1)[-1]


def style_map(style: str) -> dict[str, str]:
    return {
        key.strip().lower(): value.strip().lower()
        for item in style.split(";")
        if item.strip() and ":" in item
        for key, value in [item.split(":", 1)]
    }


def with_full_opacity(style: str) -> str:
    values = style_map(style)
    values["fill-opacity"] = "1"
    values["stroke-opacity"] = "1"
    values.pop("opacity", None)
    return ";".join(f"{key}:{value}" for key, value in values.items())


def colour_family(style: str) -> str:
    fill = style_map(style).get("fill", "")
    if fill in {"#0070b0", "#0072b2", "rgb(0%,44.705882%,69.803922%)"}:
        return "blue"
    if fill in {"#ff8733", "#ff8833", "rgb(100%,53.333333%,20%)"}:
        return "orange"
    if fill in {"#009c73", "#009e73", "rgb(0%,61.960784%,45.098039%)"}:
        return "green"
    return "other"


def opacity(style: str) -> float:
    values = style_map(style)
    for key in ("fill-opacity", "stroke-opacity", "opacity"):
        if key in values:
            try:
                return float(values[key])
            except ValueError:
                pass
    return 1.0


def first_point(d: str) -> tuple[float, float]:
    pattern = r"[Mm]\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)[, ]+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
    match = re.search(pattern, d)
    if not match:
        raise ValueError(f"Could not find the first point in SVG path: {d[:80]!r}")
    return float(match.group(1)), float(match.group(2))


def parse_viewbox(root: ET.Element) -> tuple[float, float]:
    values = [
        float(value)
        for value in re.split(r"[ ,]+", root.attrib.get("viewBox", "").strip())
        if value
    ]
    if len(values) != 4 or values[2] <= 0 or values[3] <= 0:
        raise ValueError("SVG viewBox must contain a positive width and height")
    return values[2], values[3]


def path_candidates(root: ET.Element) -> list[ET.Element]:
    width, _ = parse_viewbox(root)
    return [
        element
        for element in root.iter()
        if local_name(element.tag) == "path"
        and len(element.attrib.get("d", "")) > 1000
        and "fill-opacity" in element.attrib.get("style", "")
        and colour_family(element.attrib.get("style", "")) != "other"
        and first_point(element.attrib["d"])[0] < width * 0.60
    ]


def head_key(element: ET.Element) -> str:
    family = colour_family(element.attrib.get("style", ""))
    if family == "green" and opacity(element.attrib.get("style", "")) < 0.9:
        return "green_uncertain"
    if family == "green":
        return "green_high"
    return family


def keyed_heads(root: ET.Element) -> dict[str, ET.Element]:
    result: dict[str, ET.Element] = {}
    for element in path_candidates(root):
        key = head_key(element)
        if key in result:
            raise ValueError(f"Duplicate lollipop head key {key!r}")
        result[key] = element
    required = {"blue", "orange", "green_uncertain", "green_high"}
    missing = required.difference(result)
    if missing:
        raise ValueError(
            f"Lollipop SVG is missing head paths: {', '.join(sorted(missing))}"
        )
    return result


def by_id(root: ET.Element, identifier: str) -> ET.Element:
    matches = [
        element for element in root.iter() if element.attrib.get("id") == identifier
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one SVG element with id {identifier!r}, found {len(matches)}"
        )
    return matches[0]


def remove_ids(root: ET.Element, identifiers: set[str]) -> None:
    parents = {child: parent for parent in root.iter() for child in parent}
    for identifier in sorted(identifiers):
        element = by_id(root, identifier)
        parents[element].remove(element)


def set_stems(root: ET.Element, key: str, target_y: float) -> None:
    layout = HEAD_LAYOUT[key]
    bottom = target_y + 15.78645
    if key == "blue":
        solid = (
            f"M 355.84387,339.15184 V 317.25081 L 338.00534,304.62059 V {bottom:.5f}"
        )
        dotted = f"M 338.00534,{target_y - 0.92330:.5f} V 138.86018"
    elif key == "green_uncertain":
        solid = (
            f"M 355.84387,339.15184 V 317.25081 L 373.67721,304.62059 V {bottom:.5f}"
        )
        dotted = f"M 373.67721,{target_y - 0.99346:.5f} V 138.86018"
    else:
        solid = f"M 462.04179,339.15184 V {bottom:.5f}"
        dotted = f"M 462.04179,{target_y:.5f} V 138.86018"
    by_id(root, str(layout["solid_id"])).attrib["d"] = solid
    by_id(root, str(layout["dotted_id"])).attrib["d"] = dotted


def update_legend(root: ET.Element) -> None:
    blue = by_id(root, "path1007-7")
    green = by_id(root, "path1025-3")
    blue.attrib["transform"] = "matrix(1.3333333,0,0,-1.3333333,0,407.51815)"
    green.attrib["transform"] = "matrix(1.3333333,0,0,-1.3333333,0,407.51815)"
    blue.attrib["style"] = with_full_opacity(blue.attrib.get("style", ""))
    green.attrib["style"] = with_full_opacity(green.attrib.get("style", ""))

    cjd23 = by_id(root, "text871-7")
    cjd23.attrib["y"] = "-247.3419"
    for child in cjd23:
        if local_name(child.tag) == "tspan":
            child.attrib["y"] = "-247.3419"


def add_required_annotations(root: ET.Element, lod_percent: float) -> float:
    annotation_layer = by_id(root, "layer1")
    existing_text = {
        "".join(element.itertext()).strip()
        for element in root.iter()
        if local_name(element.tag) == "text"
    }
    required = {"VAF (%)", f"LoD ({lod_percent:.2f}%)"}
    duplicates = required.intersection(existing_text)
    if duplicates:
        raise ValueError(
            f"Reviewed SVG already contains: {', '.join(sorted(duplicates))}"
        )

    root.attrib["viewBox"] = "-18 0 882 432"
    root.attrib["width"] = "882"
    line_y = OLD_LOD_Y + Y_PER_VAF_PERCENT * (lod_percent - OLD_LOD_PERCENT)
    by_id(root, "path925-1").attrib["d"] = f"M 86.211186,{line_y:.5f} H 815.9897"

    y_axis_label = ET.SubElement(
        annotation_layer,
        f"{{{SVG_NS}}}text",
        {
            "id": "figure8-y-axis-label",
            "x": "-2",
            "y": "225",
            "transform": "rotate(-90 -2 225)",
            "text-anchor": "middle",
            "style": f"{ANNOTATION_STYLE};font-size:16px",
        },
    )
    y_axis_label.text = "VAF (%)"

    lod_label = ET.SubElement(
        annotation_layer,
        f"{{{SVG_NS}}}text",
        {
            "id": "figure8-lod-label",
            "x": "817.344",
            "y": f"{line_y - 16.32896:.5f}",
            "text-anchor": "end",
            "style": f"{ANNOTATION_STYLE};font-size:16px;font-weight:bold",
        },
    )
    lod_label.text = f"LoD ({lod_percent:.2f}%)"
    return line_y


def find_inkscape(explicit: str | None) -> str:
    candidates = [explicit, os.environ.get("INKSCAPE_PATH")]
    candidates.extend(
        [
            "/mnt/c/Program Files/Inkscape/bin/inkscape.com",
            "/mnt/c/Program Files/Inkscape/bin/inkscape.exe",
            r"C:\Program Files\Inkscape\bin\inkscape.com",
            r"C:\Program Files\Inkscape\bin\inkscape.exe",
        ]
    )
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return candidate
    raise FileNotFoundError(
        "Inkscape was not found; cannot regenerate the lollipop PDF"
    )


def render_pdf(inkscape: str, svg_path: Path, pdf_path: Path) -> None:
    def inkscape_path(path: Path) -> str:
        resolved = str(path.resolve())
        match = re.fullmatch(r"/mnt/([A-Za-z])/(.*)", resolved)
        if match and inkscape.lower().endswith((".com", ".exe")):
            return f"{match.group(1).upper()}:/{match.group(2)}"
        return resolved

    result = subprocess.run(
        [
            inkscape,
            inkscape_path(svg_path),
            "--export-type=pdf",
            f"--export-filename={inkscape_path(pdf_path)}",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Inkscape PDF export failed ({result.returncode}): {result.stderr.strip()}"
        )


def build(args: argparse.Namespace) -> None:
    rows = read_rows(args.input)
    validate_rows(rows)
    if abs(args.lod_percent - EXPECTED_LOD_PERCENT) > 1e-12:
        raise ValueError(f"Unexpected LoD percentage: {args.lod_percent!r}")
    for source in (args.baseline, args.reviewed_template):
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(
                f"Lollipop SVG source is missing or empty: {source}"
            )

    baseline_root = ET.parse(args.baseline).getroot()
    output_root = ET.parse(args.reviewed_template).getroot()
    baseline_width, baseline_height = parse_viewbox(baseline_root)
    reviewed_width, reviewed_height = parse_viewbox(output_root)
    sx = reviewed_width / baseline_width
    sy = reviewed_height / baseline_height
    if abs(sx - sy) > 1e-9:
        raise ValueError(f"Baseline/reviewed SVGs have unequal scales: {sx} and {sy}")

    baseline_heads = keyed_heads(baseline_root)
    remove_ids(output_root, OBSOLETE_IDS)
    update_legend(output_root)
    changes: list[dict[str, object]] = []
    for key, layout in HEAD_LAYOUT.items():
        baseline_head = baseline_heads[key]
        reviewed_head = by_id(output_root, str(layout["head_id"]))
        baseline_origin = first_point(baseline_head.attrib["d"])
        reviewed_origin = first_point(reviewed_head.attrib["d"])
        target_y = float(layout["old_y"]) + Y_PER_VAF_PERCENT * (
            float(layout["new_vaf"]) - float(layout["old_vaf"])
        )
        tx = reviewed_origin[0] - baseline_origin[0] * sx
        ty = target_y - baseline_origin[1] * sy
        reviewed_head.attrib["d"] = baseline_head.attrib["d"]
        reviewed_head.attrib["transform"] = (
            f"matrix({sx:.7g},0,0,{sy:.7g},{tx:.7g},{ty:.7g})"
        )
        reviewed_head.attrib["style"] = with_full_opacity(
            reviewed_head.attrib.get("style", "")
        )
        set_stems(output_root, key, target_y)
        changes.append(
            {
                "key": key,
                "sample": "CJD2" if key == "blue" else "CJD23",
                "position": 4694249 if key == "green_high" else 4691920,
                "vaf_percent": layout["new_vaf"],
                "target_head_y": target_y,
                "baseline_path_sha256": hashlib.sha256(
                    baseline_head.attrib["d"].encode()
                ).hexdigest(),
                "scale": [sx, sy],
                "translation": [tx, ty],
            }
        )

    lod_line_y = add_required_annotations(output_root, args.lod_percent)
    output_root.attrib[f"{{{SODIPODI_NS}}}docname"] = "SNV_lollipop.svg"
    output_root.insert(
        0,
        ET.Comment(
            "Original R head paths and reviewed layout retained; sequencing2 values applied"
        ),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    svg_path = args.output_dir / "SNV_lollipop.svg"
    pdf_path = args.output_dir / "SNV_lollipop.pdf"
    diff_path = args.output_dir / "SNV_lollipop_edit_diff.json"
    ET.ElementTree(output_root).write(svg_path, encoding="utf-8", xml_declaration=True)
    report = {
        "baseline_svg": str(args.baseline),
        "reviewed_template_svg": str(args.reviewed_template),
        "canonical_svg": str(svg_path),
        "canonical_pdf": str(pdf_path),
        "data_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "lod_percent": args.lod_percent,
        "lod_line_y": lod_line_y,
        "removed_ids": sorted(OBSOLETE_IDS),
        "head_edits": changes,
    }
    diff_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    render_pdf(find_inkscape(args.inkscape), svg_path, pdf_path)

    for path in (svg_path, pdf_path, diff_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Lollipop output is missing or empty: {path}")
    print(f"Recovered baseline: {args.baseline}")
    print(f"Applied {len(changes)} sequencing2 head edits and removed CJD6")
    print(f"Wrote {svg_path} ({svg_path.stat().st_size} bytes)")
    print(f"Wrote {pdf_path} ({pdf_path.stat().st_size} bytes)")
    print(f"Wrote {diff_path} ({diff_path.stat().st_size} bytes)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--reviewed-template", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--lod-percent", type=float, required=True)
    parser.add_argument("--inkscape")
    build(parser.parse_args())


if __name__ == "__main__":
    main()
