from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path("scratch/ddpcr_gating_method_comparison_v2")
PLOTS = ROOT / "plots"
TABLES = ROOT / "tables"

INPUTS = [
    (
        "A",
        "LoB/LoD pass matrix",
        PLOTS / "individual" / "method_summary" / "method_lob_lod_heatmap.svg",
    ),
    (
        "B",
        "Classifier failure rates",
        PLOTS / "individual" / "method_summary" / "method_failure_rates.svg",
    ),
]

OUTPUT = PLOTS / "panels" / "method_summary_panel.svg"
CHECKS = TABLES / "method_svg_panel_e2e_checks.csv"


def read_svg(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8")
    open_match = re.search(r"<svg\b([^>]*)>", text)
    close_match = re.search(r"</svg>\s*$", text)
    view_box_match = re.search(r"viewBox=['\"]([^'\"]+)['\"]", text)
    if open_match is None or close_match is None or view_box_match is None:
        raise ValueError(f"Could not parse SVG root/viewBox: {path}")

    view_box = [float(value) for value in view_box_match.group(1).split()]
    if len(view_box) != 4:
        raise ValueError(f"Unexpected SVG viewBox: {path}")

    return {
        "path": path,
        "text": text,
        "view_box": view_box,
        "body": text[open_match.end() : close_match.start()],
    }


def prefix_ids(svg_body: str, prefix: str) -> str:
    id_map: dict[str, str] = {}

    def replace_id(match: re.Match[str]) -> str:
        quote = match.group(1)
        old_id = match.group(2)
        new_id = f"{prefix}_{old_id}"
        id_map[old_id] = new_id
        return f"id={quote}{new_id}{quote}"

    svg_body = re.sub(r"id=(['\"])([^'\"]+)\1", replace_id, svg_body)

    for old_id, new_id in id_map.items():
        svg_body = svg_body.replace(f"url(#{old_id})", f"url(#{new_id})")
        svg_body = svg_body.replace(f"href='#{old_id}'", f"href='#{new_id}'")
        svg_body = svg_body.replace(f'href="#{old_id}"', f'href="#{new_id}"')
        svg_body = svg_body.replace(f"xlink:href='#{old_id}'", f"xlink:href='#{new_id}'")
        svg_body = svg_body.replace(f'xlink:href="#{old_id}"', f'xlink:href="#{new_id}"')

    return svg_body


def panel_svg(inputs: list[tuple[str, str, dict[str, object]]]) -> str:
    margin = 28.0
    title_height = 24.0
    gap = 34.0
    panel_width = max(item[2]["view_box"][2] for item in inputs)  # type: ignore[index]
    available_width = panel_width - (2 * margin)

    sections = []
    y_cursor = margin
    for panel_label, title, svg in inputs:
        view_box = svg["view_box"]  # type: ignore[assignment]
        source_width = float(view_box[2])
        source_height = float(view_box[3])
        scale = available_width / source_width
        rendered_height = source_height * scale
        body = prefix_ids(str(svg["body"]), f"panel_{panel_label}")
        sections.append(
            "\n".join(
                [
                    f"<text x='{margin:.2f}' y='{y_cursor + 12:.2f}' "
                    "font-family='Helvetica, Arial, sans-serif' font-size='14' "
                    "font-weight='bold'>"
                    f"{panel_label}. {title}</text>",
                    (
                        f"<g transform='translate({margin:.2f},{y_cursor + title_height:.2f}) "
                        f"scale({scale:.6f})'>"
                    ),
                    body,
                    "</g>",
                ]
            )
        )
        y_cursor += title_height + rendered_height + gap

    panel_height = y_cursor + margin - gap
    return "\n".join(
        [
            "<?xml version='1.0' encoding='UTF-8' ?>",
            (
                "<svg xmlns='http://www.w3.org/2000/svg' "
                "xmlns:xlink='http://www.w3.org/1999/xlink' "
                f"width='{panel_width:.2f}pt' height='{panel_height:.2f}pt' "
                f"viewBox='0 0 {panel_width:.2f} {panel_height:.2f}'>"
            ),
            "<rect width='100%' height='100%' fill='#FFFFFF'/>",
            *sections,
            "</svg>",
            "",
        ]
    )


def write_checks(check_rows: list[dict[str, object]]) -> None:
    CHECKS.parent.mkdir(parents=True, exist_ok=True)
    with CHECKS.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "passed", "details"])
        writer.writeheader()
        writer.writerows(check_rows)


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    parsed_inputs = [(label, title, read_svg(path)) for label, title, path in INPUTS]
    OUTPUT.write_text(panel_svg(parsed_inputs), encoding="utf-8")

    input_checks = [
        {
            "check": f"input_{index}_exists",
            "passed": path.exists() and path.stat().st_size > 0,
            "details": str(path),
        }
        for index, (_, _, path) in enumerate(INPUTS, start=1)
    ]
    output_text = OUTPUT.read_text(encoding="utf-8")
    checks = [
        *input_checks,
        {
            "check": "output_svg_exists",
            "passed": OUTPUT.exists() and OUTPUT.stat().st_size > 0,
            "details": str(OUTPUT),
        },
        {
            "check": "output_svg_contains_two_panels",
            "passed": output_text.count("<g transform=") == len(INPUTS),
            "details": str(len(INPUTS)),
        },
        {
            "check": "inlined_clip_ids_are_prefixed",
            "passed": "url(#cp" not in output_text and "id='cp" not in output_text,
            "details": "avoids duplicate svglite clip-path ids",
        },
    ]
    write_checks(checks)

    failed = [row for row in checks if not row["passed"]]
    if failed:
        failed_names = ", ".join(str(row["check"]) for row in failed)
        raise SystemExit(f"SVG panel E2E checks failed: {failed_names}")

    print(f"panel_svg={OUTPUT}")
    print(f"checks={CHECKS}")


if __name__ == "__main__":
    main()
