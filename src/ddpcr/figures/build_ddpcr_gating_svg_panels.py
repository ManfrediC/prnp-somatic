#!/usr/bin/env python3
"""Assemble ddPCR gating figure panels from individual SVG plots."""

from __future__ import annotations

import copy
import csv
import math
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET


# ---- relative paths ----
PROJECT_ROOT = Path(__file__).resolve().parents[3]
PANEL_DIR = PROJECT_ROOT / "results" / "ddPCR" / "panels"
POSITIVE_SCATTERPLOT_DIR = PROJECT_ROOT / "results" / "ddPCR" / "scatterplots" / "lob_lod_positive"
POSITIVE_MANIFEST = POSITIVE_SCATTERPLOT_DIR / "plot_manifest.csv"
STRATEGY_DIR = PROJECT_ROOT / "manuscript" / "figures" / "ddpcr_gating_strategy"
STRATEGY_MANIFEST = STRATEGY_DIR / "plot_manifest.csv"
TOOL_ENV = PROJECT_ROOT / "env" / "ddpcr.figure-tools.env"

# ---- SVG namespaces and path patterns ----
SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)

# Match internal SVG references and the two absolute-path styles emitted by
# Windows and WSL runs of the ddPCR figure workflow.
URL_REF_RE = re.compile(r"url\(#([^)]+)\)")
WSL_MOUNT_RE = re.compile(r"^/mnt/([A-Za-z])/(.*)$")
WINDOWS_DRIVE_RE = re.compile(r"^([A-Za-z]):[\\/](.*)$")
NUMBER_RE = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")

# ---- shared visual settings ----

POSITIVE_PANEL_WIDTH = 1243.3007
POSITIVE_PANEL_HEIGHT = 814.26233
POSITIVE_PANEL_WIDTH_MARGIN = 8
STRATEGY_PANEL_LABEL_OFFSET_X = 4.0
STRATEGY_PANEL_LABEL_OFFSET_Y = 20.8
POSITIVE_PANEL_TRANSFORMS = [
    "matrix(1.03333,0,0,1.03333,0.45103926,41.8125)",
    "matrix(1.03333,0,0,1.03333,398.45104,41.8125)",
    "matrix(1.03333,0,0,1.03333,796.45104,41.8125)",
    "matrix(1.03333,0,0,1.03333,0.45103926,441.8125)",
    "matrix(1.03333,0,0,1.03333,398.45104,441.8125)",
]
POSITIVE_PANEL_LEGEND_TRANSFORM = "matrix(1.4835665,0,0,1.4835665,420.93552,367.84585)"
POSITIVE_PANEL_LABELS = [
    ("A", 6.1182089, 62.601131),
    ("B", 402.55573, 62.601131),
    ("C", 801.10779, 62.861546),
    ("D", 4.5765443, 462.60114),
    ("E", 402.56613, 462.60114),
]
POSITIVE_WELL_STATUS_LABEL_X = 287.58289
POSITIVE_WELL_STATUS_LABEL_Y = 63.478531
POSITIVE_WELL_STATUS_LABEL_SIZE = "14px"
POSITIVE_WELL_STATUS_LABEL_LINE_HEIGHT = 22
POSITIVE_WELL_STATUS_LABEL_STYLE = "font-size:18.0646px"
# Keep legend colours in sync with the R scripts that generate the individual
# gating SVGs.
CLASS_COLOURS = {
    "Double-positive": "#E69F00",
    "Reference-only": "#009E73",
    "Mutant-only": "#CC79A7",
    "Double-negative": "#5684E9",
}
STRATEGY_PANEL_TITLE = "ddPCR gating strategy based on controls"


# ---- manifest input ----

def read_manifest(path: Path) -> list[dict[str, str]]:
    """Read the R-generated plot manifest that lists individual panel SVGs."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


# ---- local tool configuration ----

def read_env_file(path: Path) -> dict[str, str]:
    """Read a minimal KEY=VALUE env file without changing process state."""
    values: dict[str, str] = {}
    if not path.exists():
        return values

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip().strip("\"'")
    return values


def configured_inkscape_paths() -> list[str]:
    """Return optional local Inkscape candidates from env/ddpcr.figure-tools.env."""
    config = read_env_file(TOOL_ENV)
    raw_paths = config.get("INKSCAPE_PATHS", "")
    return [path.strip() for path in raw_paths.split(";") if path.strip()]


# ---- low-level SVG drawing helpers ----

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


def styled_text(
    parent: ET.Element,
    x: float,
    y: float,
    value: str,
    size: str,
    weight: str = "bold",
    style: str | None = None,
) -> None:
    """Add text with the accepted panel typography and Inkscape-facing style."""
    attrs = {
        "x": str(x),
        "y": str(y),
        "font-family": "Arial, Helvetica, sans-serif",
        "font-size": size,
        "font-weight": weight,
        "fill": "#111827",
    }
    if style:
        attrs["style"] = style
    node = ET.SubElement(parent, f"{{{SVG_NS}}}text", attrs)
    node.text = value


def right_aligned_text(
    parent: ET.Element,
    x: float,
    y: float,
    value: str,
    size: str = "12px",
    weight: str = "bold",
    fill: str = "#111827",
    style: str | None = None,
) -> None:
    """Add compact right-aligned SVG text for in-plot status labels."""
    attrs = {
        "x": f"{x:g}",
        "y": f"{y:g}",
        "font-family": "Arial, Helvetica, sans-serif",
        "font-size": size,
        "font-weight": weight,
        "fill": fill,
        "text-anchor": "end",
    }
    if style:
        attrs["style"] = style
    node = ET.SubElement(parent, f"{{{SVG_NS}}}text", attrs)
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


# ---- shared legend ----

def draw_strategy_legend(root: ET.Element, centre_x: float, y: float) -> None:
    """Draw the shared droplet-class and gate legend below strategy panels."""
    # Keep the legend compact enough to sit below the 12-panel strategy grid.
    x = centre_x - 345
    text(root, x, y + 12, "Droplet class", size=12, weight="bold")
    class_items = list(CLASS_COLOURS.items())
    class_x = x + 110
    class_gap_x = 150
    for index, (label, colour) in enumerate(class_items):
        item_x = class_x + index * class_gap_x
        item_y = y + 10
        circle(root, item_x, item_y - 4, 3.5, colour)
        text(root, item_x + 12, item_y, label, size=10)

    gate_x = class_x + len(class_items) * class_gap_x + 26
    text(root, gate_x, y + 12, "Gate", size=12, weight="bold")
    line(root, gate_x + 52, y + 6, gate_x + 92, y + 6, width=1.2, dasharray="0.75 2.25")
    text(root, gate_x + 102, y + 10, "Cluster-boundary gate", size=10)


# ---- SVG embedding helpers ----

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


def parse_matrix(transform: str) -> tuple[float, float, float]:
    """Return scale and x/y translations for a pure matrix(a,0,0,a,e,f) transform."""
    values = [float(value.group(0)) for value in NUMBER_RE.finditer(transform)]
    if len(values) != 6:
        raise ValueError(f"Unexpected transform format: {transform}")
    a, b, c, d, e, f = values
    if abs(b) > 1e-9 or abs(c) > 1e-9:
        raise ValueError(f"Expected axis-aligned panel transform, got: {transform}")
    if abs(a - d) > 1e-9:
        raise ValueError(f"Expected uniform scale in panel transform, got: {transform}")
    return a, e, f


def parse_finite_number(value: str) -> float | None:
    """Parse numeric attributes with optional unit suffixes."""
    cleaned = value.strip()
    match = re.fullmatch(r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)", cleaned)
    if match:
        return float(match.group(1))
    # Handle values like "12px" and "10%".
    try:
        return float(re.sub(r"[^0-9.+-eE].*$", "", cleaned))
    except ValueError:
        return None


def update_bounds(current: float | None, candidate: float | None) -> float | None:
    """Track the max bound from optional numeric candidates."""
    if candidate is None or not math.isfinite(candidate):
        return current
    return candidate if current is None else max(current, candidate)


def svg_content_max_x(root: ET.Element) -> float | None:
    """Estimate rightmost x coordinate from source SVG drawing coordinates."""
    max_x: float | None = None
    for node in root.iter():
        tag = node.tag.rsplit("}", maxsplit=1)[-1]

        if tag == "path":
            coordinates = [float(match.group(0)) for match in NUMBER_RE.finditer(node.attrib.get("d", ""))]
            if coordinates:
                max_x = update_bounds(max_x, max(coordinates[0::2]))

        elif tag == "line":
            max_x = update_bounds(max_x, parse_finite_number(node.attrib.get("x1", "")))
            max_x = update_bounds(max_x, parse_finite_number(node.attrib.get("x2", "")))

        elif tag in {"polyline", "polygon"}:
            points = [float(match.group(0)) for match in NUMBER_RE.finditer(node.attrib.get("points", ""))]
            if points:
                max_x = update_bounds(max_x, max(points[0::2]))

        elif tag == "rect":
            x = parse_finite_number(node.attrib.get("x", ""))
            width = parse_finite_number(node.attrib.get("width", ""))
            if x is not None and width is not None:
                max_x = update_bounds(max_x, x + width)
            else:
                max_x = update_bounds(max_x, x)

        elif tag == "circle":
            cx = parse_finite_number(node.attrib.get("cx", ""))
            r = parse_finite_number(node.attrib.get("r", ""))
            max_x = update_bounds(max_x, cx + r if cx is not None and r is not None else cx)

        elif tag == "ellipse":
            cx = parse_finite_number(node.attrib.get("cx", ""))
            rx = parse_finite_number(node.attrib.get("rx", ""))
            max_x = update_bounds(max_x, cx + rx if cx is not None and rx is not None else cx)

        elif tag == "text":
            text_x = parse_finite_number(node.attrib.get("x", ""))
            if text_x is not None:
                text = (node.text or "")
                approximate_chars = len(text.strip()) if text else 0
                # Conservative glyph-space fallback for labels and tick values.
                max_x = update_bounds(max_x, text_x + approximate_chars * 6.5)

        elif tag == "use":
            use_x = parse_finite_number(node.attrib.get("x", ""))
            # The glyph glyph widths vary; add a conservative fallback.
            max_x = update_bounds(max_x, use_x + 24.0 if use_x is not None else None)

    if max_x is None:
        raise ValueError("No extractable geometry in source SVG")
    return max_x


def transformed_max_x(svg_path: Path, transform: str) -> float | None:
    """Estimate max x after applying a fixed-scale panel transform."""
    source_root = ET.parse(svg_path).getroot()
    local_max_x = svg_content_max_x(source_root)
    scale, tx, _ = parse_matrix(transform)
    return scale * local_max_x + tx


def prefix_svg_ids(node: ET.Element, prefix: str) -> None:
    """Namespace embedded SVG IDs so repeated glyphs and clips do not collide."""
    for child in node.iter():
        # ggplot-generated SVGs reuse IDs such as clip paths; prefixing avoids
        # cross-panel references after multiple plots are combined.
        if "id" in child.attrib:
            child.attrib["id"] = f"{prefix}{child.attrib['id']}"
        for key, value in list(child.attrib.items()):
            if key == "id":
                continue
            if isinstance(value, str):
                value = URL_REF_RE.sub(lambda match: f"url(#{prefix}{match.group(1)})", value)
                xlink_href = f"{{{XLINK_NS}}}href"
                if key in {"href", xlink_href} and value.startswith("#"):
                    value = f"#{prefix}{value[1:]}"
                child.attrib[key] = value


def is_droplet_path(node: ET.Element) -> bool:
    """Recognise individual ggplot droplet paths from style and path shape."""
    fill = node.attrib.get("fill")
    stroke = node.attrib.get("stroke")
    d = node.attrib.get("d", "")
    return (
        node.tag == f"{{{SVG_NS}}}path"
        and node.attrib.get("stroke-width") == "0.0075"
        and fill is not None
        and fill == stroke
        and " C " in d
    )


def droplet_style_key(node: ET.Element) -> tuple[tuple[str, str], ...]:
    """Group droplet paths only when their visual style is identical."""
    return tuple(sorted((key, value) for key, value in node.attrib.items() if key != "d"))


def merge_droplet_paths(group: ET.Element) -> int:
    """Merge same-style droplet paths into compound paths for faster editing."""
    children = list(group)
    buckets: dict[tuple[tuple[str, str], ...], list[tuple[int, ET.Element]]] = defaultdict(list)
    for index, child in enumerate(children):
        if is_droplet_path(child):
            buckets[droplet_style_key(child)].append((index, child))

    replacements: dict[int, ET.Element] = {}
    skip_indexes: set[int] = set()
    merged_count = 0
    for items in buckets.values():
        if len(items) < 2:
            continue
        first_index, first_path = items[0]
        merged = copy.deepcopy(first_path)
        merged.attrib["d"] = " ".join(path.attrib["d"].strip() for _, path in items)
        merged.attrib["id"] = f"{group.attrib.get('id', 'panel')}_merged_droplets_{merged_count + 1}"
        replacements[first_index] = merged
        skip_indexes.update(index for index, _ in items[1:])
        merged_count += 1

    if replacements:
        group[:] = [replacements.get(index, child) for index, child in enumerate(children) if index not in skip_indexes]
    return merged_count


def is_plot_border(node: ET.Element) -> bool:
    """Recognise the ggplot plot-area border path by style and geometry."""
    d = node.attrib.get("d", "")
    values = [float(match.group(0)) for match in NUMBER_RE.finditer(d)]
    x_values = values[0::2]
    y_values = values[1::2]
    width = max(x_values) - min(x_values) if x_values else 0.0
    height = max(y_values) - min(y_values) if y_values else 0.0
    return (
        node.tag == f"{{{SVG_NS}}}path"
        and node.attrib.get("fill") == "none"
        and node.attrib.get("stroke") in {"rgb(20%, 20%, 20%)", "#333333"}
        and " Z " in f" {d} "
        and width > 200
        and height > 200
    )


def ensure_plot_border_drawn_last(group: ET.Element) -> None:
    """Keep the source plot border present and above merged droplets."""
    parent = {child: parent_node for parent_node in group.iter() for child in parent_node}
    borders = [child for child in group.iter(f"{{{SVG_NS}}}path") if is_plot_border(child)]
    if len(borders) != 1:
        raise RuntimeError(f"Expected exactly one plot border in {group.attrib.get('id')}; found {len(borders)}")
    border = borders[0]
    parent[border].remove(border)
    group.append(border)


def referenced_ids(root: ET.Element) -> set[str]:
    """Collect SVG IDs still referenced after filtering."""
    ids: set[str] = set()
    xlink_href = f"{{{XLINK_NS}}}href"
    for node in root.iter():
        for key, value in node.attrib.items():
            if key == "id" or not isinstance(value, str):
                continue
            ids.update(URL_REF_RE.findall(value))
            if key in {"href", xlink_href} and value.startswith("#"):
                ids.add(value[1:])
    return ids


def prune_unused_defs(root: ET.Element) -> int:
    """Remove unused definitions after plot/legend filtering."""
    used = referenced_ids(root)
    removed = 0
    for defs in list(root.iter(f"{{{SVG_NS}}}defs")):
        kept = []
        for child in list(defs):
            child_id = child.attrib.get("id")
            if child_id and child_id not in used:
                removed += 1
                continue
            kept.append(child)
        defs[:] = kept
    return removed


def append_source_svg(group: ET.Element, source_svg: Path, prefix: str) -> None:
    """Copy one panel-ready ggplot SVG component into an editable group."""
    source_root = ET.parse(source_svg).getroot()
    copied = 0
    for child in list(source_root):
        cloned = copy.deepcopy(child)
        prefix_svg_ids(cloned, prefix)
        group.append(cloned)
        copied += 1
    if copied == 0:
        raise RuntimeError(f"No SVG content copied from {source_svg}")


def inline_svg(
    parent: ET.Element,
    source_svg: Path,
    x: float,
    y: float,
    width: float,
    height: float,
    prefix: str,
) -> ET.Element:
    """Copy one source SVG into the output panel while preserving its aspect ratio."""
    source_root = ET.parse(source_svg).getroot()
    min_x, min_y, source_width, source_height = svg_viewbox(source_root)

    # Fit the source plot into its cell without distorting axis geometry.
    scale = min(width / source_width, height / source_height)
    tx = x + (width - source_width * scale) / 2 - min_x * scale
    ty = y + (height - source_height * scale) / 2 - min_y * scale

    # Copy the source SVG children rather than linking external files so the
    # resulting manuscript panel is self-contained.
    wrapper = ET.SubElement(
        parent,
        f"{{{SVG_NS}}}g",
        {"transform": f"matrix({scale:g} 0 0 {scale:g} {tx:g} {ty:g})"},
    )
    for child in list(source_root):
        cloned = copy.deepcopy(child)
        prefix_svg_ids(cloned, prefix)
        wrapper.append(cloned)
    return wrapper


# ---- LoB/LoD status labels ----

def manifest_bool(row: dict[str, str], field: str) -> bool:
    """Parse a required logical field from the R-written manifest."""
    value = (row.get(field) or "").strip().lower()
    if value in {"true", "t", "1", "yes"}:
        return True
    if value in {"false", "f", "0", "no"}:
        return False
    raise ValueError(f"Manifest row has no valid {field} value: {row}")


def lob_lod_status_lines(row: dict[str, str]) -> tuple[str, str]:
    """Return compact LoD/LoB status labels for one individual well."""
    lod = "LoD+" if manifest_bool(row, "detected_above_LoD") else "LoD-"
    lob = "LoB+" if manifest_bool(row, "detected_above_LoB") else "LoB-"
    return lob, lod


def add_lob_lod_status_label(group: ET.Element, row: dict[str, str]) -> None:
    """Draw the individual well LoD/LoB call in source plot coordinates."""
    for line_index, label in enumerate(lob_lod_status_lines(row)):
        right_aligned_text(
            group,
            POSITIVE_WELL_STATUS_LABEL_X,
            POSITIVE_WELL_STATUS_LABEL_Y + POSITIVE_WELL_STATUS_LABEL_LINE_HEIGHT * line_index,
            label,
            size=POSITIVE_WELL_STATUS_LABEL_SIZE,
            style=POSITIVE_WELL_STATUS_LABEL_STYLE,
        )


# ---- manifest path normalisation ----

def manifest_svg_path(row: dict[str, str]) -> Path:
    """Return the SVG path from either supported manifest schema."""
    svg_path = row.get("output_svg_path") or row.get("svg_path")
    if not svg_path:
        raise ValueError(f"Manifest row has no SVG path: {row}")
    path = Path(svg_path)
    if path.exists():
        return path

    # Convert WSL manifests for Windows-side execution.
    wsl_match = WSL_MOUNT_RE.match(svg_path)
    if wsl_match:
        drive, rest = wsl_match.groups()
        windows_path = Path(f"{drive.upper()}:/{rest}")
        if windows_path.exists():
            return windows_path

    # Convert Windows manifests for WSL-side execution.
    windows_match = WINDOWS_DRIVE_RE.match(svg_path)
    if windows_match:
        drive, rest = windows_match.groups()
        rest = rest.replace("\\", "/")
        wsl_path = Path(f"/mnt/{drive.lower()}/{rest}")
        if wsl_path.exists():
            return wsl_path

    return path


def write_lob_lod_positive_panel(
    rows: list[dict[str, str]],
    legend_row: dict[str, str],
    output_svg: Path,
) -> Path:
    """Build the accepted LoB/LoD panel from current ggplot source SVGs."""
    if len(rows) != 5:
        raise RuntimeError(f"Expected 5 LoB/LoD individual-well rows; found {len(rows)}")

    panel_row_paths = [manifest_svg_path(row) for row in rows]
    legend_path = manifest_svg_path(legend_row)

    panel_width = POSITIVE_PANEL_WIDTH
    right_bound: float | None = None
    width_is_estimated = True
    for source_path, transform in zip(panel_row_paths, POSITIVE_PANEL_TRANSFORMS, strict=True):
        try:
            right_bound = update_bounds(right_bound, transformed_max_x(source_path, transform))
        except Exception:
            width_is_estimated = False
            break
    if width_is_estimated:
        try:
            right_bound = update_bounds(right_bound, transformed_max_x(legend_path, POSITIVE_PANEL_LEGEND_TRANSFORM))
        except Exception:
            width_is_estimated = False
    if width_is_estimated and right_bound is not None and right_bound > 0:
        panel_width = max(POSITIVE_PANEL_WIDTH_MARGIN, math.ceil(right_bound + POSITIVE_PANEL_WIDTH_MARGIN))

    output_svg.parent.mkdir(parents=True, exist_ok=True)
    root = ET.Element(
        f"{{{SVG_NS}}}svg",
        {
            "width": str(panel_width),
            "height": str(POSITIVE_PANEL_HEIGHT),
            "viewBox": f"0 0 {panel_width} {POSITIVE_PANEL_HEIGHT}",
            "version": "1.1",
        },
    )

    for index, (row_path, transform) in enumerate(zip(panel_row_paths, POSITIVE_PANEL_TRANSFORMS, strict=True), start=1):
        group = ET.SubElement(
            root,
            f"{{{SVG_NS}}}g",
            {"id": f"editable_panel_{index}", "transform": transform},
        )
        append_source_svg(group, row_path, f"panel{index}_")
        merge_droplet_paths(group)
        add_lob_lod_status_label(group, rows[index - 1])
        ensure_plot_border_drawn_last(group)

    legend_group = ET.SubElement(
        root,
        f"{{{SVG_NS}}}g",
        {"id": "editable_panel_6", "transform": POSITIVE_PANEL_LEGEND_TRANSFORM},
    )
    append_source_svg(legend_group, legend_path, "legend_")

    for label, x, y in POSITIVE_PANEL_LABELS:
        styled_text(root, x, y, label, "18px", style="font-size:21.3333px")

    prune_unused_defs(root)

    tmp_svg = output_svg.with_name(f".{output_svg.name}.tmp")
    ET.ElementTree(root).write(str(tmp_svg), encoding="utf-8", xml_declaration=True)
    tmp_svg.replace(output_svg)
    return output_svg


# ---- panel assembly ----

def write_panel(
    rows: list[dict[str, str]],
    output_svg: Path,
    ncols: int,
    cell_width: int,
    cell_height: int,
    title: str,
    common_legend: str | None = None,
    centre_incomplete_last_row: bool = False,
) -> Path:
    """Assemble a grid of individual SVG plots into one labelled panel."""
    # Keep layout constants local so different panel types can share the same
    # assembly machinery without hiding manuscript dimensions in globals.
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
    styled_text(root, margin_x, 30, title, "22px")

    # Lay out source SVGs in reading order and add panel letters in the margin.
    for index, row in enumerate(rows):
        grid_row, grid_col = divmod(index, ncols)
        row_count = min(ncols, len(rows) - grid_row * ncols)
        row_offset = 0.0
        # Five positive wells read better as a centred 3-over-2 panel.
        if centre_incomplete_last_row and row_count < ncols:
            row_offset = (ncols - row_count) * (cell_width + gap_x) / 2
        x = margin_x + row_offset + grid_col * (cell_width + gap_x)
        y = margin_y + grid_row * (cell_height + gap_y)
        letter = chr(ord("A") + index)
        wrapper = inline_svg(root, manifest_svg_path(row), x + 24, y, cell_width - 24, cell_height, f"p{index}_")
        merge_droplet_paths(wrapper)
        ensure_plot_border_drawn_last(wrapper)
        transform = wrapper.attrib.get("transform", "matrix(1 0 0 1 0 0)")
        _, panel_x, panel_y = parse_matrix(transform)
        label_x = panel_x + STRATEGY_PANEL_LABEL_OFFSET_X
        label_y = panel_y + STRATEGY_PANEL_LABEL_OFFSET_Y
        if common_legend != "strategy":
            label_x = x + 2
            label_y = y + 24
        styled_text(root, label_x, label_y, letter, "18px", style="font-size:21.3333px")

    # Strategy panels use one shared legend to avoid repeating legends in every cell.
    if common_legend == "strategy":
        draw_strategy_legend(root, width / 2, height - legend_height + 10)

    prune_unused_defs(root)

    # Write via a same-directory temporary file so large SVG panels are not
    # left half-written if export is interrupted.
    tmp_svg = output_svg.with_name(f".{output_svg.name}.tmp")
    ET.ElementTree(root).write(str(tmp_svg), encoding="utf-8", xml_declaration=True)
    tmp_svg.replace(output_svg)
    return output_svg


# ---- PDF export ----

def export_pdf(svg_path: Path) -> Path:
    """Export the assembled SVG panel to PDF using Inkscape."""
    pdf_path = svg_path.with_suffix(".pdf")
    # Prefer PATH for portability, then fall back to local tool paths declared
    # in the ignored env-side configuration file.
    candidates = [shutil.which("inkscape"), *configured_inkscape_paths()]
    inkscape = next((candidate for candidate in candidates if candidate and Path(candidate).exists()), None)
    if not inkscape:
        raise RuntimeError(f"Inkscape is required; add INKSCAPE_PATHS to {TOOL_ENV}")
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


# ---- row ordering ----

def sort_positive_individual_scatterplots(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Order LoB+LoD+ individual-well scatterplots by sample and replicate."""
    # The scatterplot manifest also contains merged-sample rows; keep only the
    # individual-well SVGs requested for the manuscript-facing panel.
    selected = [
        row
        for row in rows
        if row["plot_kind"] == "individual_well"
        and row.get("status") == "written"
        and manifest_svg_path(row).parent.name == "individual_wells"
    ]
    participant_order = {"CJD4": 0, "CJD21": 1}
    return sorted(
        selected,
        key=lambda row: (
            participant_order.get(row["participant"], 99),
            row["positive_id"],
            int(row["replicate_index"] or 0),
            row["run_date"],
            row["well"],
        ),
    )


def positive_legend_row(rows: list[dict[str, str]]) -> dict[str, str]:
    """Return the explicit legend component written by the R scatterplot step."""
    selected = [
        row
        for row in rows
        if row["plot_kind"] == "legend"
        and row.get("status") == "written"
        and manifest_svg_path(row).parent.name == "legend"
    ]
    if len(selected) != 1:
        raise RuntimeError(f"Expected exactly one LoB/LoD legend row; found {len(selected)}")
    return selected[0]


def sort_strategy(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Order control panels by mutation and gating stage."""
    # Preserve the biological assay order and the control-to-final-gate story.
    assay_order = {"D178N": 0, "E200K": 1, "P102L": 2}
    stage_order = {"NTC": 0, "WT": 1, "Positive control": 2, "Final adjusted gate": 3}
    return sorted(rows, key=lambda row: (assay_order[row["assay"]], stage_order[row["stage"]]))


# ---- command entry point ----

def main() -> None:
    """Build all combined gating panels from individual plot manifests."""
    # Load manifests written by the R figure scripts.
    positive_rows = read_manifest(POSITIVE_MANIFEST)
    positive_individual = sort_positive_individual_scatterplots(positive_rows)
    positive_legend = positive_legend_row(positive_rows) if positive_individual else None

    outputs = []
    # Build the requested positive-well panel from the scatterplot outputs
    # written under results/ddPCR.
    if positive_individual:
        outputs.append(
            write_lob_lod_positive_panel(
                positive_individual,
                positive_legend,
                PANEL_DIR / "ddpcr_lob_lod_positive_individual_wells_panel.svg",
            )
        )
    else:
        print("No written LoB+LoD+ individual-well scatterplot rows; skipping positive panel")

    # Keep the older strategy panel available when its manifest has been
    # generated, but do not block the requested LoB/LoD scatterplot panel on it.
    if STRATEGY_MANIFEST.exists():
        strategy_rows = read_manifest(STRATEGY_MANIFEST)
        outputs.append(
            write_panel(
                sort_strategy(strategy_rows),
                PANEL_DIR / "ddpcr_gating_strategy_panel.svg",
                ncols=4,
                cell_width=345,
                cell_height=315,
                title=STRATEGY_PANEL_TITLE,
                common_legend="strategy",
            ),
        )
    else:
        print(f"No strategy manifest at {STRATEGY_MANIFEST}; skipping strategy panel")

    # Export every assembled SVG to PDF and report both artefacts.
    for svg_path in outputs:
        pdf_path = export_pdf(svg_path)
        print(f"Wrote {svg_path}")
        print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    # Keep imports side-effect free; only build panels when executed as a script.
    main()
