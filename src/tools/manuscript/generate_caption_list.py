#!/usr/bin/env python3
"""Generate the plain-text manuscript caption catalogue from canonical sources."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


CROSSWALK_PATH = Path("manuscript/config/figure_table_crosswalk.tsv")
OUTPUT_PATH = Path("manuscript/table_and_figure_captions.txt")

CAPTION_COMMAND = re.compile(
    r"\\caption(?![A-Za-z@])(?:\s*\[[^]]*\])?\s*\{"
    r"|\\captionof\s*\{\s*(?:figure|table)\s*\}\s*\{"
)
ARGUMENT_COMMAND = re.compile(r"\\(?:textbf|textit|emph|ensuremath)\s*\{")

SIMPLE_COMMANDS = {
    "dn": "D178N",
    "ek": "E200K",
    "pl": "P102L",
    "prnp": "PRNP",
    "prpsc": "PrPSc",
    "textgreater": ">",
    "textmu": "µ",
    "xfold": "×",
}


def _is_escaped(text: str, index: int) -> bool:
    backslashes = 0
    index -= 1
    while index >= 0 and text[index] == "\\":
        backslashes += 1
        index -= 1
    return backslashes % 2 == 1


def _strip_comments(text: str) -> str:
    lines: list[str] = []
    for line in text.splitlines():
        for index, character in enumerate(line):
            if character == "%" and not _is_escaped(line, index):
                line = line[:index]
                break
        lines.append(line)
    return "\n".join(lines)


def _braced_argument(text: str, opening_brace: int) -> tuple[str, int]:
    depth = 0
    for index in range(opening_brace, len(text)):
        character = text[index]
        if character == "{" and not _is_escaped(text, index):
            depth += 1
        elif character == "}" and not _is_escaped(text, index):
            depth -= 1
            if depth == 0:
                return text[opening_brace + 1 : index], index + 1
    raise ValueError("Unbalanced braces in caption source")


def extract_outer_caption(path: Path) -> str:
    """Return the final caption command, after any preceding subcaptions."""
    source = _strip_comments(path.read_text(encoding="utf-8"))
    matches = list(CAPTION_COMMAND.finditer(source))
    if not matches:
        raise ValueError(f"No caption found in {path}")
    match = matches[-1]
    caption, _ = _braced_argument(source, match.end() - 1)
    return caption


def _unwrap_argument_commands(text: str) -> str:
    while match := ARGUMENT_COMMAND.search(text):
        argument, end = _braced_argument(text, match.end() - 1)
        text = text[: match.start()] + argument + text[end:]
    return text


def tex_to_plain_text(text: str) -> str:
    """Convert the small TeX vocabulary used by the canonical captions."""
    text = _unwrap_argument_commands(text)

    for command, replacement in SIMPLE_COMMANDS.items():
        text = re.sub(
            rf"\\{command}(?![A-Za-z@])(?:\s*\{{\s*\}})?",
            lambda _match, value=replacement: value,
            text,
        )

    replacements = {
        r"\%": "%",
        r"\&": "&",
        r"\_": "_",
        r"\#": "#",
        r"\geq": "≥",
        r"\leq": "≤",
        r"\,": " ",
        r"\;": " ",
        r"\ ": " ",
    }
    for source, replacement in replacements.items():
        text = text.replace(source, replacement)

    text = text.replace("R^2", "R²")
    text = text.replace("``", "“").replace("''", "”")
    text = text.replace("--", "–")
    text = text.replace("~", " ")
    text = text.replace("$", "")
    text = text.replace("{", "").replace("}", "")
    text = re.sub(r"\s+", " ", text).strip()

    residual = re.search(r"\\[A-Za-z@]+|\\.|[$^{}]", text)
    if residual:
        raise ValueError(f"Unsupported TeX remains in caption near {residual.group(0)!r}: {text}")
    return text


def _external_caption_path(root: Path, row: dict[str, str]) -> Path | None:
    assets = [item.strip() for item in row["asset_output"].split(";") if item.strip()]
    workbooks = [Path(item) for item in assets if Path(item).suffix.lower() == ".xlsx"]
    if not workbooks:
        return None
    if len(workbooks) != 1:
        raise ValueError(f"Expected one workbook for {row['docx_reference']}")
    return (root / workbooks[0]).with_suffix(".caption.txt")


def _caption_for_row(root: Path, row: dict[str, str]) -> str:
    external_caption = _external_caption_path(root, row)
    if external_caption is not None:
        if not external_caption.is_file():
            raise ValueError(f"Missing workbook caption sidecar: {external_caption}")
        caption = re.sub(r"\s+", " ", external_caption.read_text(encoding="utf-8")).strip()
        if not caption:
            raise ValueError(f"Empty workbook caption sidecar: {external_caption}")
        return caption

    tex_path = root / row["tex_file"]
    return tex_to_plain_text(extract_outer_caption(tex_path))


def build_caption_catalogue(root: Path) -> str:
    crosswalk = root / CROSSWALK_PATH
    with crosswalk.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    sections = (
        ("table", "TABLE CAPTIONS", "=============="),
        ("figure", "FIGURE CAPTIONS", "==============="),
    )
    lines: list[str] = []
    for kind, heading, underline in sections:
        if lines:
            lines.append("")
        lines.extend((heading, underline, ""))
        section_rows = [row for row in rows if row["type"] == kind]
        for index, row in enumerate(section_rows):
            reference = row["docx_reference"].removeprefix("Supplementary ")
            lines.append(f"{reference}. {_caption_for_row(root, row)}")
            if index != len(section_rows) - 1:
                lines.append("")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--check", action="store_true", help="Fail if the catalogue is stale")
    args = parser.parse_args()

    root = args.root.resolve()
    output = root / OUTPUT_PATH
    try:
        expected = build_caption_catalogue(root).encode("utf-8")
    except (KeyError, OSError, ValueError) as exc:
        print(f"CAPTION CATALOGUE ERROR: {exc}", file=sys.stderr)
        return 2

    if args.check:
        if not output.is_file() or output.read_bytes() != expected:
            print(f"CAPTION CATALOGUE IS STALE: run {Path(__file__).name}")
            return 1
        print("CAPTION CATALOGUE IS CURRENT")
        return 0

    output.write_bytes(expected)
    print(f"Wrote {output.relative_to(root)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
