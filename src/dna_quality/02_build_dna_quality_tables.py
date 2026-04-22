#!/usr/bin/env python3
from __future__ import annotations

# ---------------------------------------------------------------------------
# Standard library imports
# ---------------------------------------------------------------------------

import argparse
import csv
import datetime as dt
import io
import math
import os
import posixpath
import re
import shutil
import subprocess
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET


# ---------------------------------------------------------------------------
# Default external roots and output schemas
# ---------------------------------------------------------------------------

DEFAULT_SURESELECT_ROOT = (
    "/mnt/c/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/"
    "CJD PRNP/Experiments/SureSelect-sequencing/Experiments"
)
DEFAULT_DDPCR_ROOT = (
    "/mnt/c/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/"
    "CJD PRNP/Experiments/ddPCR"
)
DEFAULT_SAMPLES_ROOT = (
    "/mnt/c/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/"
    "CJD PRNP/Samples"
)


FILE_INVENTORY_COLUMNS = [
    "source_family",
    "run_id",
    "stage",
    "instrument",
    "file_type",
    "parser_role",
    "parse_strategy",
    "path",
]

LIBRARY_QC_COLUMNS = [
    "run_id",
    "stage",
    "instrument",
    "parser_source",
    "source_path",
    "native_source_path",
    "pdf_source_path",
    "well",
    "sample_description",
    "sample_id",
    "sample_id_source",
    "batch_id",
    "reported_concentration",
    "reported_concentration_unit",
    "peak_row_count",
    "dominant_peak_size_bp",
    "dominant_peak_calibrated_concentration",
    "dominant_peak_molarity",
    "dominant_peak_area_pct",
    "region_from_bp",
    "region_to_bp",
    "region_average_size_bp",
    "region_concentration",
    "region_concentration_unit",
    "region_molarity",
    "region_molarity_unit",
    "region_percent_total",
    "alert",
    "observations",
]

PREP_METADATA_COLUMNS = [
    "source_family",
    "source_run_id",
    "source_path",
    "sheet_name",
    "row_index",
    "raw_sample_label",
    "sample_id",
    "code",
    "group",
    "region",
    "well_or_barcode",
    "dna_input_ng_ul",
    "clean_dna_ng_ul",
    "qubit_ng_ul",
    "pre_capture_pcr_ng_ul",
    "molarity_nm",
    "elution_volume_ul",
    "notes",
]

SAMPLE_ALIASES_COLUMNS = [
    "alias",
    "alias_key",
    "sample_id",
    "source",
    "source_run_id",
    "notes",
]

SAMPLE_QUALITY_MASTER_COLUMNS = [
    "sample_id",
    "sample_id_source",
    "batch_id",
    "run_id",
    "stage",
    "instrument",
    "parser_source",
    "source_path",
    "well",
    "sample_description",
    "reported_concentration",
    "reported_concentration_unit",
    "peak_row_count",
    "dominant_peak_size_bp",
    "dominant_peak_area_pct",
    "region_average_size_bp",
    "region_concentration",
    "region_concentration_unit",
    "region_molarity",
    "region_percent_total",
    "library_qc_heuristic_score",
    "dna_input_ng_ul",
    "clean_dna_ng_ul",
    "qubit_ng_ul",
    "pre_capture_pcr_ng_ul",
    "prep_molarity_nm",
    "prep_source_families",
    "prep_regions",
    "sequencing_metrics_path",
    "sequencing_on_target_fraction",
    "sequencing_target_mean_depth",
    "sequencing_target_fold80",
    "sequencing_pct_duplication",
    "sequencing_estimated_library_size",
    "sequencing_outcome_band",
    "alert",
    "observations",
]

SCORECARD_COLUMNS = [
    "sample_id",
    "batch_id",
    "library_run_count",
    "stages",
    "mean_reported_concentration",
    "mean_region_average_size_bp",
    "mean_region_percent_total",
    "max_clean_dna_ng_ul",
    "max_qubit_ng_ul",
    "max_pre_capture_pcr_ng_ul",
    "mean_library_qc_heuristic_score",
    "sequencing_on_target_fraction",
    "sequencing_target_mean_depth",
    "sequencing_target_fold80",
    "sequencing_pct_duplication",
    "sequencing_outcome_band",
]

ANALYSIS_SUMMARY_COLUMNS = [
    "summary_type",
    "group_name",
    "metric",
    "n",
    "mean",
    "median",
    "minimum",
    "maximum",
]


# ---------------------------------------------------------------------------
# CLI and string-normalisation helpers
# ---------------------------------------------------------------------------

# Parse the small CLI surface that mirrors the old PowerShell entrypoint.
def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    parser = argparse.ArgumentParser(
        description=(
            "Build descriptive DNA-quality proxy tables by harmonising Tapestation, "
            "prep metadata, ddPCR quantity records, and sequencing QC."
        )
    )
    parser.add_argument("--sureselect-root", default=DEFAULT_SURESELECT_ROOT)
    parser.add_argument("--ddpcr-root", default=DEFAULT_DDPCR_ROOT)
    parser.add_argument("--samples-root", default=DEFAULT_SAMPLES_ROOT)
    parser.add_argument("--output-run", default="")
    parser.add_argument("--repo-root", default=str(repo_root))
    parser.add_argument("--skip-pdf-extraction", action="store_true")
    return parser.parse_args()


# Convert Windows-style paths when the workflow is launched from WSL, while
# still accepting native Linux and native Windows path forms.
def normalize_cli_path(value: str) -> str:
    text = value.strip()
    if os.name != "nt" and re.match(r"^[A-Za-z]:[\\/]", text):
        drive = text[0].lower()
        tail = text[2:].replace("\\", "/")
        return f"/mnt/{drive}{tail}"
    if os.name == "nt" and re.match(r"^/mnt/[A-Za-z]/", text):
        drive = text[5].upper()
        tail = text[6:].replace("/", "\\")
        return f"{drive}:{tail}"
    return text


# Collapse repeated whitespace and non-breaking spaces so joins and regexes
# behave consistently across spreadsheets, CSVs, and PDF text.
def normalize_whitespace(value: object) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", " ", str(value).replace("\xa0", " ")).strip()


# Normalise labels into a simple lowercase lookup key that tolerates common
# encoding noise from spreadsheet exports.
def normalize_key(value: object) -> str:
    text = normalize_whitespace(value)
    text = text.replace("\u00b5", "u").replace("\u03bc", "u").replace("\u00c2", "")
    return re.sub(r"[^a-z0-9]+", "", text.lower())


# Normalise aliases into an uppercase alphanumeric key for sample matching.
def get_alias_key(value: object) -> str:
    return re.sub(r"[^A-Z0-9]+", "", normalize_whitespace(value).upper())


# Pull canonical sample IDs directly from free-text labels when possible.
def get_sample_id_hint(text: object) -> str:
    normalized = normalize_whitespace(text)
    if not normalized:
        return ""
    for pattern in (
        r"(?i)(CJD\d+)(?=$|[^0-9A-Z]|_)",
        r"(?i)(CTRL\d+)(?=$|[^0-9A-Z]|_)",
        r"(?i)(NA[0-9A-Z_]+)(?=$|[^0-9A-Z]|_)",
    ):
        match = re.search(pattern, normalized)
        if match:
            return match.group(1).upper()
    return ""


# Normalise sample IDs for case-insensitive joins across manifests and QC files.
def get_sample_id_key(value: object) -> str:
    return normalize_whitespace(value).upper()


# Return the first non-empty value whose column name matches one of the
# expected aliases used across the project spreadsheets.
def get_first_property_value(row: dict[str, object], candidates: list[str]) -> str:
    for candidate in candidates:
        target = normalize_key(candidate)
        for key, value in row.items():
            if normalize_key(key) == target:
                text = normalize_whitespace(value)
                if text:
                    return text
    return ""


# Keep object-field access centralized so downstream helpers return blanks
# instead of raising on missing optional fields.
def get_object_property_value(row: dict[str, object] | None, field: str) -> object:
    if row is None:
        return ""
    return row.get(field, "")


# ---------------------------------------------------------------------------
# Run and stage inference helpers
# ---------------------------------------------------------------------------

# Tokenize run IDs into informative fragments so ambiguous aliases can be
# broken using nearby run-context clues.
def get_run_tokens(run_id: str) -> list[str]:
    if not run_id:
        return []
    stop_words = {
        "SS",
        "QC",
        "SEQUENCE",
        "SAMPLE",
        "SAMPLES",
        "PART",
        "RUN",
        "FINAL",
        "FIRST",
        "SECOND",
        "THIRD",
        "DAY",
        "BATCH",
        "OF",
        "THE",
        "MIRKA",
    }
    tokens = []
    for token in run_id.upper().split("_"):
        if token and token not in stop_words and token not in tokens:
            tokens.append(token)
    return tokens


# Score the overlap between two run IDs so alias matches prefer contextually
# similar prep and QC runs when names collide.
def get_run_overlap_score(left: str, right: str) -> int:
    if not left or not right:
        return 0
    right_lookup = set(get_run_tokens(right))
    score = 0
    for token in get_run_tokens(left):
        if token in right_lookup:
            if re.fullmatch(r"\d{4}", token):
                score += 3
            elif re.fullmatch(r"\d+", token):
                score += 2
            else:
                score += 1
    return score


# Map directory names into the three high-level library QC stages used in the
# descriptive summaries.
def get_stage_for_path(path: Path) -> str:
    text = str(path).lower()
    if "1st qc" in text or "first day" in text:
        return "pre_capture"
    if "post-capture" in text:
        return "post_capture"
    if "qc before sample submission" in text or "sample submission" in text:
        return "submission_qc"
    if " second sequencing run " in text or text.endswith("second sequencing run"):
        return "submission_qc"
    return "unknown"


# Build a stable run identifier from the file's location under the external
# source root, falling back to the filename when needed.
def get_run_id(path: Path, root: Path) -> str:
    parent = path.parent
    try:
        relative = parent.relative_to(root)
        if str(relative) != ".":
            return re.sub(r"[^A-Za-z0-9]+", "_", str(relative)).strip("_")
    except ValueError:
        pass
    return re.sub(r"[^A-Za-z0-9]+", "_", path.stem).strip("_")


# ---------------------------------------------------------------------------
# XLSX parsing helpers
# ---------------------------------------------------------------------------

# Read one XML member from an XLSX archive, returning blank when that member
# is absent in a given workbook.
def get_zip_entry_text(archive: zipfile.ZipFile, entry_name: str) -> str:
    try:
        with archive.open(entry_name) as handle:
            return handle.read().decode("utf-8")
    except KeyError:
        return ""


# Strip XML namespaces so workbook parsing stays readable.
def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


# Decode the workbook shared-string table so sheet cells can be materialised
# into plain text values.
def get_spreadsheet_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    text = get_zip_entry_text(archive, "xl/sharedStrings.xml")
    if not text:
        return []
    root = ET.fromstring(text)
    values: list[str] = []
    for node in root.iter():
        if local_name(node.tag) != "si":
            continue
        parts = [part.text or "" for part in node.iter() if local_name(part.tag) == "t"]
        values.append(normalize_whitespace("".join(parts)))
    return values


# Resolve visible sheet names to the ZIP member that stores each worksheet XML.
def get_spreadsheet_sheets(archive: zipfile.ZipFile) -> list[dict[str, str]]:
    workbook_text = get_zip_entry_text(archive, "xl/workbook.xml")
    rels_text = get_zip_entry_text(archive, "xl/_rels/workbook.xml.rels")
    if not workbook_text or not rels_text:
        return []

    workbook = ET.fromstring(workbook_text)
    rels = ET.fromstring(rels_text)
    rel_map: dict[str, str] = {}
    for relationship in rels.iter():
        if local_name(relationship.tag) != "Relationship":
            continue
        rel_id = relationship.attrib.get("Id", "")
        target = relationship.attrib.get("Target", "")
        if rel_id and target:
            rel_map[rel_id] = target

    sheets: list[dict[str, str]] = []
    for sheet in workbook.iter():
        if local_name(sheet.tag) != "sheet":
            continue
        rid = ""
        for key, value in sheet.attrib.items():
            if key == "id" or key.endswith("}id"):
                rid = value
                break
        target = rel_map.get(rid, "")
        if not target:
            continue
        entry_path = target.lstrip("/")
        if not entry_path.startswith("xl/"):
            entry_path = posixpath.normpath(posixpath.join("xl", entry_path))
        sheets.append(
            {
                "Name": normalize_whitespace(sheet.attrib.get("name", "")),
                "EntryPath": entry_path,
            }
        )
    return sheets


# Convert Excel column letters like "AA" into one-based numeric positions.
def get_cell_column_index(reference: str) -> int:
    match = re.match(r"^([A-Z]+)", reference or "")
    if not match:
        return 0
    index = 0
    for letter in match.group(1):
        index = (index * 26) + (ord(letter) - ord("A") + 1)
    return index


# Decode one worksheet cell, including inline strings and shared-string cells.
def get_cell_text(cell: ET.Element, shared_strings: list[str]) -> str:
    cell_type = cell.attrib.get("t", "")
    if cell_type == "inlineStr":
        parts = [part.text or "" for part in cell.iter() if local_name(part.tag) == "t"]
        return normalize_whitespace("".join(parts))

    value_node = next((node for node in cell if local_name(node.tag) == "v"), None)
    if value_node is None:
        return ""
    raw = normalize_whitespace(value_node.text or "")
    if cell_type == "s" and re.fullmatch(r"\d+", raw):
        index = int(raw)
        if index < len(shared_strings):
            return shared_strings[index]
    return raw


# Materialise a simple row-oriented view of an XLSX workbook without adding a
# pandas dependency to the WSL workflow.
def get_xlsx_rows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with zipfile.ZipFile(path) as archive:
        shared_strings = get_spreadsheet_shared_strings(archive)
        sheets = get_spreadsheet_sheets(archive)
        for sheet in sheets:
            sheet_text = get_zip_entry_text(archive, sheet["EntryPath"])
            if not sheet_text:
                continue
            sheet_xml = ET.fromstring(sheet_text)
            header_names: list[str] | None = None
            for row in sheet_xml.iter():
                if local_name(row.tag) != "row":
                    continue
                cell_map: dict[int, str] = {}
                for cell in row:
                    if local_name(cell.tag) != "c":
                        continue
                    column_index = get_cell_column_index(cell.attrib.get("r", ""))
                    if column_index <= 0:
                        continue
                    cell_map[column_index] = get_cell_text(cell, shared_strings)
                if not cell_map:
                    continue

                max_column = max(cell_map)
                values = [cell_map.get(column, "") for column in range(1, max_column + 1)]

                if header_names is None:
                    header_names = []
                    seen: set[str] = set()
                    for index, value in enumerate(values, start=1):
                        header = normalize_whitespace(value) or f"column_{index}"
                        candidate = header
                        suffix = 2
                        while candidate.lower() in seen:
                            candidate = f"{header}_{suffix}"
                            suffix += 1
                        seen.add(candidate.lower())
                        header_names.append(candidate)
                    continue

                if not any(values):
                    continue

                record: dict[str, object] = {
                    "source_path": str(path),
                    "sheet_name": sheet["Name"],
                    "row_index": int(row.attrib.get("r", "0") or 0),
                }
                for index, header in enumerate(header_names):
                    record[header] = values[index] if index < len(values) else ""
                rows.append(record)
    return rows


# ---------------------------------------------------------------------------
# CSV and Tapestation parsing helpers
# ---------------------------------------------------------------------------

# Try the encodings seen in the project exports so WSL runs can read the same
# raw files as the Windows workflow.
def read_delimited_rows(path: Path, delimiter: str) -> list[dict[str, str]]:
    last_error: UnicodeDecodeError | None = None
    for encoding in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            with path.open("r", encoding=encoding, newline="") as handle:
                text = handle.read()
            reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
            return [dict(row) for row in reader]
        except UnicodeDecodeError as exc:
            last_error = exc
    if last_error is not None:
        raise last_error
    return []


# Parse numeric text fields while leaving blanks untouched for downstream TSVs.
def convert_to_number_or_blank(value: object) -> float | str:
    text = normalize_whitespace(value)
    if not text:
        return ""
    text = text.replace(",", ".")
    try:
        return float(text)
    except ValueError:
        return ""


# Return the first sorted match so repeated glob lookups stay deterministic.
def find_first(globbed: list[Path]) -> Path | None:
    return sorted(globbed, key=lambda item: str(item))[0] if globbed else None


# Capture the paired native Tapestation file for provenance reporting.
def get_run_native_pair(files: list[Path]) -> str:
    for file_path in sorted(files, key=lambda item: str(item)):
        if file_path.suffix.lower() in {".d1000", ".hsd1000"}:
            return str(file_path)
    return ""


# Parse one Tapestation CSV bundle into row-level library QC observations.
def parse_tapestation_csv_bundle(directory: Path, run_id: str, stage: str) -> list[dict[str, object]]:
    sample_csv = find_first(list(directory.glob("*_sampleTable.csv")))
    peak_csv = find_first(list(directory.glob("*_compactPeakTable.csv")))
    region_csv = find_first(list(directory.glob("*_compactRegionTable.csv")))
    pdf_pair = find_first(list(directory.glob("*.pdf")))
    native_pair = get_run_native_pair([path for path in directory.iterdir() if path.is_file()])
    if sample_csv is None:
        return []

    samples = read_delimited_rows(sample_csv, ",")
    peaks = read_delimited_rows(peak_csv, ",") if peak_csv else []
    regions = read_delimited_rows(region_csv, ",") if region_csv else []

    peak_map: dict[str, list[dict[str, str]]] = defaultdict(list)
    for peak in peaks:
        well = normalize_whitespace(peak.get("Well", ""))
        peak_comment = normalize_whitespace(peak.get("Peak Comment", ""))
        if well == "EL1" or re.search("Marker", peak_comment, flags=re.IGNORECASE):
            continue
        peak_map[well].append(peak)

    region_map = {
        normalize_whitespace(region.get("WellId", "")): region
        for region in regions
        if normalize_whitespace(region.get("WellId", ""))
    }

    rows: list[dict[str, object]] = []
    for sample in samples:
        well = normalize_whitespace(sample.get("Well", ""))
        if well == "EL1":
            continue
        peak_rows = peak_map.get(well, [])
        dominant_peak = None
        if peak_rows:
            dominant_peak = max(
                enumerate(peak_rows),
                key=lambda item: (
                    convert_to_number_or_blank(item[1].get("% Integrated Area", "")) or 0,
                    item[0],
                ),
            )[1]
        region = region_map.get(well)
        rows.append(
            {
                "run_id": run_id,
                "stage": stage,
                "instrument": "D1000",
                "parser_source": "csv",
                "source_path": str(sample_csv),
                "native_source_path": native_pair,
                "pdf_source_path": str(pdf_pair) if pdf_pair else "",
                "well": well,
                "sample_description": normalize_whitespace(sample.get("Sample Description", "")),
                "reported_concentration": convert_to_number_or_blank(
                    get_first_property_value(sample, ["Conc. [ng/ul]"])
                ),
                "reported_concentration_unit": "ng/ul",
                "peak_row_count": len(peak_rows),
                "dominant_peak_size_bp": (
                    convert_to_number_or_blank(dominant_peak.get("Size [bp]", ""))
                    if dominant_peak
                    else ""
                ),
                "dominant_peak_calibrated_concentration": (
                    convert_to_number_or_blank(
                        get_first_property_value(dominant_peak, ["Calibrated Conc. [ng/ul]"])
                    )
                    if dominant_peak
                    else ""
                ),
                "dominant_peak_molarity": (
                    convert_to_number_or_blank(dominant_peak.get("Peak Molarity [nmol/l]", ""))
                    if dominant_peak
                    else ""
                ),
                "dominant_peak_area_pct": (
                    convert_to_number_or_blank(dominant_peak.get("% Integrated Area", ""))
                    if dominant_peak
                    else ""
                ),
                "region_from_bp": (
                    convert_to_number_or_blank(region.get("From [bp]", "")) if region else ""
                ),
                "region_to_bp": (
                    convert_to_number_or_blank(region.get("To [bp]", "")) if region else ""
                ),
                "region_average_size_bp": (
                    convert_to_number_or_blank(region.get("Average Size [bp]", "")) if region else ""
                ),
                "region_concentration": (
                    convert_to_number_or_blank(
                        get_first_property_value(region, ["Conc. [ng/ul]"])
                    )
                    if region
                    else ""
                ),
                "region_concentration_unit": "ng/ul",
                "region_molarity": (
                    convert_to_number_or_blank(region.get("Region Molarity [nmol/l]", ""))
                    if region
                    else ""
                ),
                "region_molarity_unit": "nmol/l",
                "region_percent_total": (
                    convert_to_number_or_blank(region.get("% of Total", "")) if region else ""
                ),
                "alert": normalize_whitespace(sample.get("Alert", "")),
                "observations": normalize_whitespace(sample.get("Observations", "")),
            }
        )
    return rows


# Extract line-oriented text from PDFs, preferring Ghostscript so the output
# matches the original PowerShell workflow as closely as possible.
def get_pdf_lines(path: Path, skip_pdf_extraction: bool) -> list[str]:
    if skip_pdf_extraction:
        return []
    if shutil.which("gs"):
        command = [
            "gs",
            "-q",
            "-dNOPAUSE",
            "-dBATCH",
            "-sDEVICE=txtwrite",
            "-sOutputFile=-",
            str(path),
        ]
    elif shutil.which("pdftotext"):
        command = ["pdftotext", "-layout", str(path), "-"]
    else:
        raise RuntimeError(
            "PDF extraction requires either `gs` or `pdftotext` on PATH. "
            "Use --skip-pdf-extraction to skip PDF parsing."
        )
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return result.stdout.splitlines()


# Parse a Tapestation PDF report when no CSV sidecars are available for that run.
def parse_tapestation_pdf(
    path: Path, run_id: str, stage: str, skip_pdf_extraction: bool
) -> list[dict[str, object]]:
    lines = get_pdf_lines(path, skip_pdf_extraction)
    if not lines:
        return []
    instrument = "HSD1000" if any("High Sensitivity D1000" in line for line in lines) else "D1000"
    conc_unit = "pg/ul" if instrument == "HSD1000" else "ng/ul"
    molarity_unit = "pmol/l" if instrument == "HSD1000" else "nmol/l"
    native_pair = get_run_native_pair([item for item in path.parent.iterdir() if item.is_file()])
    rows: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    section = ""

    def flush_current() -> None:
        nonlocal current
        if current and current.get("well") != "EL1":
            rows.append(current)
        current = None

    for line in lines:
        trimmed = normalize_whitespace(line)
        if not trimmed:
            continue

        match = re.match(r"^(EL1|[A-H]\d{1,2}):\s*(.+)$", trimmed)
        if match:
            flush_current()
            current = {
                "run_id": run_id,
                "stage": stage,
                "instrument": instrument,
                "parser_source": "pdf",
                "source_path": str(path),
                "native_source_path": native_pair,
                "pdf_source_path": str(path),
                "well": match.group(1),
                "sample_description": normalize_whitespace(match.group(2)),
                "reported_concentration": "",
                "reported_concentration_unit": conc_unit,
                "peak_row_count": 0,
                "dominant_peak_size_bp": "",
                "dominant_peak_calibrated_concentration": "",
                "dominant_peak_molarity": "",
                "dominant_peak_area_pct": "",
                "region_from_bp": "",
                "region_to_bp": "",
                "region_average_size_bp": "",
                "region_concentration": "",
                "region_concentration_unit": conc_unit,
                "region_molarity": "",
                "region_molarity_unit": molarity_unit,
                "region_percent_total": "",
                "alert": "",
                "observations": "",
            }
            section = ""
            continue

        if current is None:
            continue
        if trimmed == "Sample Table":
            section = "sample"
            continue
        if trimmed == "Peak Table":
            section = "peak"
            continue
        if trimmed == "Region Table":
            section = "region"
            continue
        if (
            trimmed.startswith("TapeStation Analysis Software")
            or "ScreenTape" in trimmed
            or trimmed.startswith("Filename:")
            or trimmed == "Sample Info"
            or trimmed == "Default image (Contrast 100%)"
            or trimmed == "Default image (Contrast 50%), Image is Scaled to Sample"
        ):
            continue

        if section == "sample":
            sample_match = re.match(r"^(EL1|[A-H]\d{1,2})\s+(.*)$", trimmed)
            if sample_match and sample_match.group(1) == current["well"]:
                tail = normalize_whitespace(sample_match.group(2))
                leading_number = re.match(r"^([0-9]+(?:\.[0-9]+)?)\s+", tail)
                if leading_number:
                    current["reported_concentration"] = convert_to_number_or_blank(
                        leading_number.group(1)
                    )
                caution_match = re.search(r"(Caution!.*)$", tail)
                if caution_match:
                    current["observations"] = normalize_whitespace(caution_match.group(1))
            continue

        if section == "peak":
            if re.search("Marker", trimmed, flags=re.IGNORECASE):
                continue
            tokens = re.split(r"\s+", trimmed)
            if (
                instrument == "D1000"
                and len(tokens) >= 5
                and re.fullmatch(r"\d+(?:\.\d+)?", tokens[0])
                and re.fullmatch(r"\d+(?:\.\d+)?", tokens[4])
            ):
                area = convert_to_number_or_blank(tokens[4])
                current["peak_row_count"] = int(current["peak_row_count"]) + 1
                existing = current["dominant_peak_area_pct"]
                if existing == "" or float(area) > float(existing):
                    current["dominant_peak_size_bp"] = convert_to_number_or_blank(tokens[0])
                    current["dominant_peak_calibrated_concentration"] = convert_to_number_or_blank(
                        tokens[1]
                    )
                    current["dominant_peak_molarity"] = convert_to_number_or_blank(tokens[3])
                    current["dominant_peak_area_pct"] = area
            elif (
                instrument == "HSD1000"
                and len(tokens) >= 4
                and re.fullmatch(r"\d+(?:\.\d+)?", tokens[0])
                and re.fullmatch(r"\d+(?:\.\d+)?", tokens[3])
            ):
                area = convert_to_number_or_blank(tokens[3])
                current["peak_row_count"] = int(current["peak_row_count"]) + 1
                existing = current["dominant_peak_area_pct"]
                if existing == "" or float(area) > float(existing):
                    current["dominant_peak_calibrated_concentration"] = convert_to_number_or_blank(
                        tokens[0]
                    )
                    current["dominant_peak_molarity"] = convert_to_number_or_blank(tokens[2])
                    current["dominant_peak_area_pct"] = area
            continue

        if section == "region":
            region_match = re.match(
                r"^(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+"
                r"(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)",
                trimmed,
            )
            if region_match:
                current["region_from_bp"] = convert_to_number_or_blank(region_match.group(1))
                current["region_to_bp"] = convert_to_number_or_blank(region_match.group(2))
                current["region_average_size_bp"] = convert_to_number_or_blank(region_match.group(3))
                current["region_concentration"] = convert_to_number_or_blank(region_match.group(4))
                current["region_molarity"] = convert_to_number_or_blank(region_match.group(5))
                current["region_percent_total"] = convert_to_number_or_blank(region_match.group(6))

    flush_current()
    return rows


# ---------------------------------------------------------------------------
# Spreadsheet metadata extraction and summary helpers
# ---------------------------------------------------------------------------

# Pull supporting prep and quantity metadata out of the SureSelect, ddPCR,
# and sample-manifest workbooks.
def get_metadata_rows(
    files: list[Path],
    source_family: str,
    sureselect_root: Path,
    ddpcr_root: Path,
    samples_root: Path,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for file_path in files:
        if source_family == "sureselect_prep":
            source_run_id = get_run_id(file_path, sureselect_root)
        elif source_family == "ddpcr_quantity":
            source_run_id = get_run_id(file_path, ddpcr_root)
        elif source_family == "sample_manifest":
            source_run_id = get_run_id(file_path, samples_root)
        else:
            source_run_id = ""

        workbook_rows = get_xlsx_rows(file_path)
        for row in workbook_rows:
            sample_pieces = [
                get_first_property_value(
                    row,
                    [
                        "new_name",
                        "new name",
                        "sample name",
                        "sample",
                        "# sample / code",
                        "sample code",
                        "name",
                        "patient",
                    ],
                ),
                get_first_property_value(row, ["code"]),
                get_first_property_value(row, ["sample description"]),
            ]
            sample_pieces = [piece for piece in sample_pieces if piece]

            sample_id = ""
            for piece in sample_pieces:
                candidate = get_sample_id_hint(piece)
                if candidate:
                    sample_id = candidate
                    break
            if not sample_id:
                for value in row.values():
                    candidate = get_sample_id_hint(value)
                    if candidate:
                        sample_id = candidate
                        break

            record = {
                "source_family": source_family,
                "source_run_id": source_run_id,
                "source_path": str(file_path),
                "sheet_name": row["sheet_name"],
                "row_index": row["row_index"],
                "raw_sample_label": get_first_property_value(
                    row,
                    [
                        "new_name",
                        "new name",
                        "sample name",
                        "sample",
                        "# sample / code",
                        "sample code",
                        "name",
                        "patient",
                    ],
                ),
                "sample_id": sample_id,
                "code": get_first_property_value(row, ["code"]),
                "group": get_first_property_value(row, ["group"]),
                "region": get_first_property_value(row, ["brain region", "region"]),
                "well_or_barcode": get_first_property_value(
                    row, ["sureselect barcode", "barcode", "index", "well"]
                ),
                "dna_input_ng_ul": convert_to_number_or_blank(
                    get_first_property_value(
                        row,
                        ["DNA fr (ng/ul)", "DNA (ng/ul)", "DNA ng/ul"],
                    )
                ),
                "clean_dna_ng_ul": convert_to_number_or_blank(
                    get_first_property_value(row, ["clean DNA (ng/ul)", "clean DNA"])
                ),
                "qubit_ng_ul": convert_to_number_or_blank(
                    get_first_property_value(
                        row,
                        ["qubit", "qubit conc. (ng/ul)", "qubit conc", "DNA conc qubit"],
                    )
                ),
                "pre_capture_pcr_ng_ul": convert_to_number_or_blank(
                    get_first_property_value(
                        row,
                        [
                            "pre-capt. PCR product conc. (ng/ul)",
                            "pre capt PCR product conc. (ng/ul)",
                        ],
                    )
                ),
                "molarity_nm": convert_to_number_or_blank(
                    get_first_property_value(
                        row,
                        ["final sample conc. / molarity [nm]", "molarity", "molarity (nm)"],
                    )
                ),
                "elution_volume_ul": convert_to_number_or_blank(
                    get_first_property_value(row, ["eluted in (ul)"])
                ),
                "notes": get_first_property_value(row, ["comments", "comment", "notes", "observations"]),
            }

            if any(
                record[field]
                for field in (
                    "raw_sample_label",
                    "sample_id",
                    "well_or_barcode",
                    "dna_input_ng_ul",
                    "clean_dna_ng_ul",
                    "qubit_ng_ul",
                    "pre_capture_pcr_ng_ul",
                    "molarity_nm",
                )
            ):
                rows.append(record)
    return rows


# Collect one numeric field from a row set while skipping blanks.
def get_numeric_values(rows: list[dict[str, object]], field: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = get_object_property_value(row, field)
        if value == "":
            continue
        values.append(float(value))
    return values


# Surface the maximum numeric value for one field across a row group.
def get_max_numeric_value(rows: list[dict[str, object]], field: str) -> float | str:
    values = get_numeric_values(rows, field)
    return max(values) if values else ""


# Surface the mean numeric value for one field across a row group.
def get_mean_numeric_value(rows: list[dict[str, object]], field: str) -> float | str:
    values = get_numeric_values(rows, field)
    if not values:
        return ""
    return round(sum(values) / len(values), 3)


# Surface the median numeric value for one field across a row group.
def get_median_numeric_value(rows: list[dict[str, object]], field: str) -> float | str:
    values = sorted(get_numeric_values(rows, field))
    if not values:
        return ""
    count = len(values)
    if count % 2 == 1:
        return values[count // 2]
    return round((values[(count // 2) - 1] + values[count // 2]) / 2, 3)


# Join distinct non-empty values into the semicolon-delimited format used by
# the existing project summaries.
def get_distinct_joined_value(rows: list[dict[str, object]], field: str) -> str:
    values = sorted(
        {
            normalize_whitespace(get_object_property_value(row, field))
            for row in rows
            if normalize_whitespace(get_object_property_value(row, field))
        }
    )
    return "; ".join(values)


# Resolve ambiguous aliases by preferring exact sample-ID agreement first and
# run-context agreement second.
def resolve_alias_match(
    alias_entries: list[dict[str, object]], candidate_values: list[object], run_id: str
) -> dict[str, str] | None:
    for candidate_value in candidate_values:
        key = get_alias_key(candidate_value)
        if not key:
            continue
        matches = [entry for entry in alias_entries if entry["alias_key"] == key]
        if not matches:
            continue

        unique_sample_ids: dict[str, str] = {}
        for match in matches:
            sample_id_text = normalize_whitespace(match["sample_id"])
            unique_sample_ids.setdefault(get_sample_id_key(sample_id_text), sample_id_text)
        if len(unique_sample_ids) == 1:
            return {
                "sample_id": next(iter(unique_sample_ids.values())),
                "source": str(matches[0]["source"]),
            }

        best_match = None
        best_score = -1
        best_sample_id_key = ""
        tie = False
        for match in matches:
            score = get_run_overlap_score(run_id, str(match["source_run_id"]))
            sample_id_key = get_sample_id_key(match["sample_id"])
            if score > best_score:
                best_score = score
                best_match = match
                best_sample_id_key = sample_id_key
                tie = False
            elif score == best_score and sample_id_key != best_sample_id_key:
                tie = True

        if best_match and not tie and best_score > 0:
            return {
                "sample_id": str(best_match["sample_id"]),
                "source": f"{best_match['source']}_context",
            }
    return None


# Match filenames using the compact case-insensitive patterns used in the
# original workflow inventory filters.
def path_matches(path: Path, pattern: str) -> bool:
    return re.search(pattern, path.name, flags=re.IGNORECASE) is not None


# Walk one external root and return the files that belong to this workflow.
def list_files(root: Path, predicate) -> list[Path]:
    if not root.exists():
        raise FileNotFoundError(f"Required input root does not exist: {root}")
    return sorted((path for path in root.rglob("*") if path.is_file() and predicate(path)), key=str)


# Format scalar values for TSV output while preserving blanks and stable float
# rendering across reruns.
def format_scalar(value: object) -> str:
    if value is None or value == "":
        return ""
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return ""
        return f"{value:.15g}"
    return str(value)


# Write a quoted TSV with a fixed column order so downstream tables remain
# diff-friendly and deterministic.
def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            quoting=csv.QUOTE_ALL,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: format_scalar(get_object_property_value(row, column)) for column in columns})


# Write the small Markdown report that summarises run scope and output paths.
def write_report(
    output_path: Path,
    output_run: str,
    sureselect_files: list[Path],
    ddpcr_files: list[Path],
    sample_files: list[Path],
    resolved_library_rows: list[dict[str, object]],
    resolved_prep_metadata: list[dict[str, object]],
    scorecard_rows: list[dict[str, object]],
) -> None:
    unresolved_library_count = sum(1 for row in resolved_library_rows if not row["sample_id"])
    report_lines = [
        "# DNA Quality Analysis",
        "",
        f"Output run: {output_run}",
        "",
        "## Scope",
        f"- SureSelect files inventoried: {len(sureselect_files)}",
        f"- ddPCR files inventoried: {len(ddpcr_files)}",
        f"- Samples metadata files inventoried: {len(sample_files)}",
        f"- Library QC rows parsed: {len(resolved_library_rows)}",
        f"- Prep metadata rows extracted: {len(resolved_prep_metadata)}",
        f"- Library QC rows still unmatched to `sample_id`: {unresolved_library_count}",
        f"- Samples with scorecards: {len(scorecard_rows)}",
        "",
        "## Notes",
        "- CSV sidecars were preferred over PDFs when both existed in the same Tapestation run directory.",
        "- Native `.D1000` / `.HSD1000` files are inventoried as provenance in `file_inventory.tsv` but are not parsed directly in this v1 workflow.",
        "- `library_qc_heuristic_score` is a transparent helper score, not a validated assay threshold.",
        "",
        "## Outputs",
        "- `file_inventory.tsv`: all useful external files discovered for this analysis.",
        "- `library_qc.tsv`: per-run Tapestation-derived QC rows.",
        "- `prep_metadata.tsv`: extracted supporting metadata from SureSelect, ddPCR, and Samples spreadsheets.",
        "- `input_dna_quantity.tsv`: subset of metadata rows with upstream quantity-related values.",
        "- `sample_aliases.tsv`: alias map used for sample-resolution joins.",
        "- `sample_quality_master.tsv`: joined library QC + upstream quantity + sequencing outcome table.",
        "- `dna_quality_scorecard.tsv`: per-sample summary for manuscript/plots.",
        "- `analysis_summary.tsv`: grouped descriptive summaries by sequencing outcome and QC stage.",
        "",
    ]
    output_path.write_text("\n".join(report_lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

# Run the full DNA-quality proxy workflow from source discovery through
# harmonised outputs and the summary report.
def main() -> int:
    args = parse_args()

    # Resolve user-supplied or default roots up front so later path handling is
    # consistent across WSL and native Windows contexts.
    output_run = args.output_run or dt.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    repo_root = Path(normalize_cli_path(args.repo_root)).resolve()
    sureselect_root = Path(normalize_cli_path(args.sureselect_root)).resolve()
    ddpcr_root = Path(normalize_cli_path(args.ddpcr_root)).resolve()
    samples_root = Path(normalize_cli_path(args.samples_root)).resolve()
    output_dir = repo_root / "results" / "dna_quality" / output_run
    output_dir.mkdir(parents=True, exist_ok=True)

    # Inventory the three external source families that feed the harmonised QC
    # tables: SureSelect runs, ddPCR quantity workbooks, and sample manifests.
    sureselect_files = list_files(
        sureselect_root,
        lambda path: path.suffix.lower() in {".pdf", ".csv", ".xlsx", ".d1000", ".hsd1000"},
    )
    sureselect_xlsx = [path for path in sureselect_files if path.suffix.lower() == ".xlsx"]
    csv_bundle_directories = sorted(
        {path.parent for path in sureselect_files if path.name.endswith("_sampleTable.csv")},
        key=str,
    )
    csv_bundle_directory_set = {str(path) for path in csv_bundle_directories}
    pdf_files_for_parsing = [
        path
        for path in sureselect_files
        if path.suffix.lower() == ".pdf" and str(path.parent) not in csv_bundle_directory_set
    ]

    ddpcr_files = list_files(
        ddpcr_root,
        lambda path: path.suffix.lower() in {".xlsx", ".docx"}
        and path_matches(
            path,
            r"Covaris and Tapestation|summary table patients ddPCR|qubit|clean DNA|qPCR - samples",
        ),
    )
    ddpcr_xlsx = [path for path in ddpcr_files if path.suffix.lower() == ".xlsx"]
    sample_files = list_files(samples_root, lambda path: path.suffix.lower() == ".xlsx")

    # Record one provenance row per discovered external file before any parsing
    # decisions are applied.
    file_inventory: list[dict[str, object]] = []
    for file_path in sureselect_files:
        suffix = file_path.suffix.lower()
        file_inventory.append(
            {
                "source_family": "sureselect",
                "run_id": get_run_id(file_path, sureselect_root),
                "stage": get_stage_for_path(file_path),
                "instrument": (
                    "HSD1000"
                    if suffix == ".hsd1000"
                    else "D1000"
                    if suffix == ".d1000" or "D1000" in file_path.name
                    else ""
                ),
                "file_type": file_path.suffix.lstrip("."),
                "parser_role": (
                    "sample_table"
                    if file_path.name.endswith("_sampleTable.csv")
                    else "peak_table"
                    if file_path.name.endswith("_compactPeakTable.csv")
                    else "region_table"
                    if file_path.name.endswith("_compactRegionTable.csv")
                    else "report"
                    if suffix == ".pdf"
                    else "native"
                    if suffix in {".d1000", ".hsd1000"}
                    else "metadata"
                ),
                "parse_strategy": (
                    "pdf_text"
                    if suffix == ".pdf" and str(file_path.parent) not in csv_bundle_directory_set
                    else "csv_bundle"
                    if suffix == ".csv"
                    else "inventory_only"
                    if suffix in {".d1000", ".hsd1000"}
                    else "spreadsheet"
                ),
                "path": str(file_path),
            }
        )
    for file_path in ddpcr_files:
        suffix = file_path.suffix.lower()
        file_inventory.append(
            {
                "source_family": "ddpcr",
                "run_id": get_run_id(file_path, ddpcr_root),
                "stage": "upstream_quantity",
                "instrument": "documentation" if suffix == ".docx" else "spreadsheet",
                "file_type": file_path.suffix.lstrip("."),
                "parser_role": "protocol_note" if suffix == ".docx" else "metadata",
                "parse_strategy": "spreadsheet" if suffix == ".xlsx" else "inventory_only",
                "path": str(file_path),
            }
        )
    for file_path in sample_files:
        file_inventory.append(
            {
                "source_family": "samples",
                "run_id": get_run_id(file_path, samples_root),
                "stage": "sample_manifest",
                "instrument": "spreadsheet",
                "file_type": file_path.suffix.lstrip("."),
                "parser_role": "metadata",
                "parse_strategy": "spreadsheet",
                "path": str(file_path),
            }
        )

    # Parse library QC observations from CSV bundles first, then fall back to
    # PDF parsing only for runs that do not ship the CSV sidecars.
    library_rows: list[dict[str, object]] = []
    for directory in csv_bundle_directories:
        run_id = get_run_id(directory, sureselect_root)
        stage = get_stage_for_path(directory)
        library_rows.extend(parse_tapestation_csv_bundle(directory, run_id, stage))
    for pdf_path in pdf_files_for_parsing:
        run_id = get_run_id(pdf_path, sureselect_root)
        stage = get_stage_for_path(pdf_path)
        library_rows.extend(parse_tapestation_pdf(pdf_path, run_id, stage, args.skip_pdf_extraction))

    # Extract prep and quantity metadata from all workbook sources into one
    # shared row model.
    prep_metadata: list[dict[str, object]] = []
    prep_metadata.extend(
        get_metadata_rows(sureselect_xlsx, "sureselect_prep", sureselect_root, ddpcr_root, samples_root)
    )
    prep_metadata.extend(
        get_metadata_rows(ddpcr_xlsx, "ddpcr_quantity", sureselect_root, ddpcr_root, samples_root)
    )
    prep_metadata.extend(
        get_metadata_rows(sample_files, "sample_manifest", sureselect_root, ddpcr_root, samples_root)
    )

    alias_entries: list[dict[str, object]] = []
    alias_map: dict[str, dict[str, object]] = {}

    # Build the alias table once so both metadata rows and library QC rows can
    # reuse the same sample-resolution rules.
    def add_alias_entry(
        alias: object, sample_id: object, source: str, source_run_id: str, notes: object
    ) -> None:
        key = get_alias_key(alias)
        sample_id_text = normalize_whitespace(sample_id)
        if not key or not sample_id_text:
            return
        if key in alias_map:
            if alias_map[key]["sample_id"] != sample_id_text:
                alias_map[key]["ambiguous"] = True
        else:
            alias_map[key] = {
                "sample_id": sample_id_text,
                "source": source,
                "ambiguous": False,
            }
        alias_entries.append(
            {
                "alias": normalize_whitespace(alias),
                "alias_key": key,
                "sample_id": sample_id_text,
                "source": source,
                "source_run_id": source_run_id,
                "notes": normalize_whitespace(notes),
            }
        )

    # Seed the alias table from parsed metadata first, then extend it with the
    # preprocessing manifest and any manual override file.
    for row in prep_metadata:
        if not row["sample_id"]:
            continue
        add_alias_entry(row["sample_id"], row["sample_id"], str(row["source_family"]), str(row["source_run_id"]), row["source_path"])
        add_alias_entry(
            row["raw_sample_label"],
            row["sample_id"],
            str(row["source_family"]),
            str(row["source_run_id"]),
            row["source_path"],
        )
        add_alias_entry(
            row["code"],
            row["sample_id"],
            str(row["source_family"]),
            str(row["source_run_id"]),
            row["source_path"],
        )
        add_alias_entry(
            row["well_or_barcode"],
            row["sample_id"],
            str(row["source_family"]),
            str(row["source_run_id"]),
            row["source_path"],
        )

    preprocessing_manifest_path = repo_root / "config" / "preprocessing_samples.tsv"
    preprocessing_rows = (
        read_delimited_rows(preprocessing_manifest_path, "\t") if preprocessing_manifest_path.exists() else []
    )
    batch_map: dict[str, str] = {}
    for row in preprocessing_rows:
        sample_id = normalize_whitespace(row.get("sample_id", ""))
        batch_id = normalize_whitespace(row.get("batch_id", ""))
        if sample_id:
            batch_map[get_sample_id_key(sample_id)] = batch_id
            add_alias_entry(sample_id, sample_id, "config_preprocessing", "", batch_id)

    override_path = repo_root / "config" / "dna_quality_sample_alias_overrides.tsv"
    if override_path.exists():
        for row in read_delimited_rows(override_path, "\t"):
            add_alias_entry(row.get("alias", ""), row.get("sample_id", ""), "alias_override", "", row.get("notes", ""))

    # Resolve sample IDs onto prep metadata rows using direct ID hints first and
    # alias matching only when the raw workbook labels are ambiguous.
    resolved_prep_metadata: list[dict[str, object]] = []
    for row in prep_metadata:
        sample_id = normalize_whitespace(row["sample_id"])
        if not sample_id:
            for candidate in (row["raw_sample_label"], row["code"], row["well_or_barcode"]):
                direct_sample_id = get_sample_id_hint(candidate)
                if direct_sample_id:
                    sample_id = direct_sample_id
                    break
        if not sample_id:
            alias_match = resolve_alias_match(
                alias_entries,
                [row["raw_sample_label"], row["code"], row["well_or_barcode"]],
                str(row["source_run_id"]),
            )
            if alias_match:
                sample_id = alias_match["sample_id"]

        resolved_row = dict(row)
        resolved_row["sample_id"] = sample_id
        resolved_prep_metadata.append(resolved_row)

    # Resolve sample IDs onto library QC rows using the same alias table and
    # carry batch information forward from preprocessing_samples.tsv.
    resolved_library_rows: list[dict[str, object]] = []
    for row in library_rows:
        sample_id = get_sample_id_hint(row["sample_description"])
        sample_id_source = ""
        if sample_id:
            sample_id_source = "inline_label"
        else:
            alias_match = resolve_alias_match(
                alias_entries,
                [row["sample_description"], row["well"]],
                str(row["run_id"]),
            )
            if alias_match:
                sample_id = alias_match["sample_id"]
                sample_id_source = alias_match["source"]

        resolved_library_rows.append(
            {
                "run_id": row["run_id"],
                "stage": (
                    "submission_qc"
                    if row["stage"] == "unknown" and row["instrument"] == "HSD1000"
                    else row["stage"]
                ),
                "instrument": row["instrument"],
                "parser_source": row["parser_source"],
                "source_path": row["source_path"],
                "native_source_path": row["native_source_path"],
                "pdf_source_path": row["pdf_source_path"],
                "well": row["well"],
                "sample_description": row["sample_description"],
                "sample_id": sample_id,
                "sample_id_source": sample_id_source,
                "batch_id": batch_map.get(get_sample_id_key(sample_id), "") if sample_id else "",
                "reported_concentration": row["reported_concentration"],
                "reported_concentration_unit": row["reported_concentration_unit"],
                "peak_row_count": row["peak_row_count"],
                "dominant_peak_size_bp": row["dominant_peak_size_bp"],
                "dominant_peak_calibrated_concentration": row["dominant_peak_calibrated_concentration"],
                "dominant_peak_molarity": row["dominant_peak_molarity"],
                "dominant_peak_area_pct": row["dominant_peak_area_pct"],
                "region_from_bp": row["region_from_bp"],
                "region_to_bp": row["region_to_bp"],
                "region_average_size_bp": row["region_average_size_bp"],
                "region_concentration": row["region_concentration"],
                "region_concentration_unit": row["region_concentration_unit"],
                "region_molarity": row["region_molarity"],
                "region_molarity_unit": row["region_molarity_unit"],
                "region_percent_total": row["region_percent_total"],
                "alert": row["alert"],
                "observations": row["observations"],
            }
        )

    # Index prep metadata and downstream sequencing metrics by sample so the
    # final joined table can be assembled in one pass over library QC rows.
    prep_by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in resolved_prep_metadata:
        if row["sample_id"]:
            prep_by_sample[get_sample_id_key(row["sample_id"])].append(row)

    seq_files = sorted(
        (path for path in (repo_root / "results").rglob("sequencing_metrics_per_sample.tsv") if path.is_file()),
        key=lambda item: item.stat().st_mtime,
        reverse=True,
    )
    preprocessing_sample_ids = {get_sample_id_key(row.get("sample_id", "")) for row in preprocessing_rows}
    preprocessing_sample_ids.discard("")
    sequencing_by_sample: dict[str, dict[str, object]] = {}
    for file_path in seq_files:
        for row in read_delimited_rows(file_path, "\t"):
            sample_id = normalize_whitespace(row.get("sample_id", ""))
            sample_id_key = get_sample_id_key(sample_id)
            if sample_id_key not in preprocessing_sample_ids or sample_id_key in sequencing_by_sample:
                continue
            enriched = dict(row)
            enriched["source_path"] = str(file_path)
            sequencing_by_sample[sample_id_key] = enriched

    # Keep a dedicated upstream-quantity export for the manuscript notes and
    # ad hoc spot checks.
    input_dna_quantity = [
        row
        for row in resolved_prep_metadata
        if any(
            row[field] != ""
            for field in (
                "dna_input_ng_ul",
                "clean_dna_ng_ul",
                "qubit_ng_ul",
                "pre_capture_pcr_ng_ul",
                "molarity_nm",
            )
        )
    ]

    # Join library QC, upstream quantity, and sequencing outcomes into the
    # main per-observation table used for later summaries.
    sample_quality_master: list[dict[str, object]] = []
    for row in resolved_library_rows:
        prep_rows = (
            prep_by_sample.get(get_sample_id_key(row["sample_id"]), [])
            if row["sample_id"]
            else []
        )
        seq = (
            sequencing_by_sample.get(get_sample_id_key(row["sample_id"]), None)
            if row["sample_id"]
            else None
        )

        library_score = 0
        if row["reported_concentration"] != "" or row["region_concentration"] != "":
            library_score += 1
        if (
            row["region_average_size_bp"] != ""
            and 300 <= float(row["region_average_size_bp"]) <= 500
        ):
            library_score += 1
        if row["region_percent_total"] != "" and float(row["region_percent_total"]) >= 85:
            library_score += 1
        if row["peak_row_count"] != "" and int(row["peak_row_count"]) <= 2:
            library_score += 1

        sequencing_band = ""
        if seq is not None:
            mean_depth = convert_to_number_or_blank(seq.get("target_mean_depth", ""))
            on_target = convert_to_number_or_blank(seq.get("on_target_fraction", ""))
            fold80 = convert_to_number_or_blank(seq.get("target_fold80", ""))
            if (
                mean_depth != ""
                and on_target != ""
                and fold80 != ""
                and float(mean_depth) >= 3000
                and float(on_target) >= 0.70
                and float(fold80) <= 1.50
            ):
                sequencing_band = "high"
            elif (
                mean_depth != ""
                and on_target != ""
                and float(mean_depth) >= 2000
                and float(on_target) >= 0.60
            ):
                sequencing_band = "moderate"
            else:
                sequencing_band = "low"

        sample_quality_master.append(
            {
                "sample_id": row["sample_id"],
                "sample_id_source": row["sample_id_source"],
                "batch_id": row["batch_id"],
                "run_id": row["run_id"],
                "stage": row["stage"],
                "instrument": row["instrument"],
                "parser_source": row["parser_source"],
                "source_path": row["source_path"],
                "well": row["well"],
                "sample_description": row["sample_description"],
                "reported_concentration": row["reported_concentration"],
                "reported_concentration_unit": row["reported_concentration_unit"],
                "peak_row_count": row["peak_row_count"],
                "dominant_peak_size_bp": row["dominant_peak_size_bp"],
                "dominant_peak_area_pct": row["dominant_peak_area_pct"],
                "region_average_size_bp": row["region_average_size_bp"],
                "region_concentration": row["region_concentration"],
                "region_concentration_unit": row["region_concentration_unit"],
                "region_molarity": row["region_molarity"],
                "region_percent_total": row["region_percent_total"],
                "library_qc_heuristic_score": library_score,
                "dna_input_ng_ul": get_max_numeric_value(prep_rows, "dna_input_ng_ul"),
                "clean_dna_ng_ul": get_max_numeric_value(prep_rows, "clean_dna_ng_ul"),
                "qubit_ng_ul": get_max_numeric_value(prep_rows, "qubit_ng_ul"),
                "pre_capture_pcr_ng_ul": get_max_numeric_value(prep_rows, "pre_capture_pcr_ng_ul"),
                "prep_molarity_nm": get_max_numeric_value(prep_rows, "molarity_nm"),
                "prep_source_families": get_distinct_joined_value(prep_rows, "source_family"),
                "prep_regions": get_distinct_joined_value(prep_rows, "region"),
                "sequencing_metrics_path": seq["source_path"] if seq else "",
                "sequencing_on_target_fraction": (
                    convert_to_number_or_blank(seq.get("on_target_fraction", "")) if seq else ""
                ),
                "sequencing_target_mean_depth": (
                    convert_to_number_or_blank(seq.get("target_mean_depth", "")) if seq else ""
                ),
                "sequencing_target_fold80": (
                    convert_to_number_or_blank(seq.get("target_fold80", "")) if seq else ""
                ),
                "sequencing_pct_duplication": (
                    convert_to_number_or_blank(seq.get("pct_duplication", "")) if seq else ""
                ),
                "sequencing_estimated_library_size": (
                    convert_to_number_or_blank(seq.get("estimated_library_size", "")) if seq else ""
                ),
                "sequencing_outcome_band": sequencing_band,
                "alert": row["alert"],
                "observations": row["observations"],
            }
        )

    # Collapse the per-observation table into one scorecard row per sample.
    scorecard_rows: list[dict[str, object]] = []
    for sample_id in sorted({str(row["sample_id"]) for row in sample_quality_master if row["sample_id"]}):
        rows = [row for row in sample_quality_master if row["sample_id"] == sample_id]
        scorecard_rows.append(
            {
                "sample_id": sample_id,
                "batch_id": get_distinct_joined_value(rows, "batch_id"),
                "library_run_count": len(rows),
                "stages": get_distinct_joined_value(rows, "stage"),
                "mean_reported_concentration": get_mean_numeric_value(rows, "reported_concentration"),
                "mean_region_average_size_bp": get_mean_numeric_value(rows, "region_average_size_bp"),
                "mean_region_percent_total": get_mean_numeric_value(rows, "region_percent_total"),
                "max_clean_dna_ng_ul": get_max_numeric_value(rows, "clean_dna_ng_ul"),
                "max_qubit_ng_ul": get_max_numeric_value(rows, "qubit_ng_ul"),
                "max_pre_capture_pcr_ng_ul": get_max_numeric_value(rows, "pre_capture_pcr_ng_ul"),
                "mean_library_qc_heuristic_score": get_mean_numeric_value(
                    rows, "library_qc_heuristic_score"
                ),
                "sequencing_on_target_fraction": get_max_numeric_value(
                    rows, "sequencing_on_target_fraction"
                ),
                "sequencing_target_mean_depth": get_max_numeric_value(
                    rows, "sequencing_target_mean_depth"
                ),
                "sequencing_target_fold80": get_max_numeric_value(rows, "sequencing_target_fold80"),
                "sequencing_pct_duplication": get_max_numeric_value(
                    rows, "sequencing_pct_duplication"
                ),
                "sequencing_outcome_band": get_distinct_joined_value(rows, "sequencing_outcome_band"),
            }
        )

    # Build grouped descriptive summaries for manuscript drafting and quick
    # sanity checks across outcome bands and library-QC stages.
    analysis_summary_rows: list[dict[str, object]] = []
    scorecard_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in scorecard_rows:
        group_name = str(row["sequencing_outcome_band"])
        if not group_name:
            continue
        scorecard_groups[group_name].append(row)
    for group_name in sorted(scorecard_groups):
        rows = scorecard_groups[group_name]
        for metric in (
            "mean_region_average_size_bp",
            "mean_region_percent_total",
            "max_clean_dna_ng_ul",
            "mean_library_qc_heuristic_score",
            "sequencing_on_target_fraction",
            "sequencing_target_mean_depth",
        ):
            numeric_values = sorted(get_numeric_values(rows, metric))
            analysis_summary_rows.append(
                {
                    "summary_type": "scorecard_by_outcome",
                    "group_name": group_name,
                    "metric": metric,
                    "n": len(numeric_values),
                    "mean": get_mean_numeric_value(rows, metric),
                    "median": get_median_numeric_value(rows, metric),
                    "minimum": numeric_values[0] if numeric_values else "",
                    "maximum": numeric_values[-1] if numeric_values else "",
                }
            )

    stage_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in sample_quality_master:
        stage_groups[str(row["stage"])].append(row)
    for group_name in sorted(stage_groups):
        rows = stage_groups[group_name]
        for metric in (
            "reported_concentration",
            "region_average_size_bp",
            "region_percent_total",
            "library_qc_heuristic_score",
        ):
            numeric_values = sorted(get_numeric_values(rows, metric))
            analysis_summary_rows.append(
                {
                    "summary_type": "library_rows_by_stage",
                    "group_name": group_name,
                    "metric": metric,
                    "n": len(numeric_values),
                    "mean": get_mean_numeric_value(rows, metric),
                    "median": get_median_numeric_value(rows, metric),
                    "minimum": numeric_values[0] if numeric_values else "",
                    "maximum": numeric_values[-1] if numeric_values else "",
                }
            )

    # Write all workflow outputs only after every table has been built so the
    # run directory is either complete or absent.
    write_tsv(output_dir / "file_inventory.tsv", file_inventory, FILE_INVENTORY_COLUMNS)
    write_tsv(output_dir / "library_qc.tsv", resolved_library_rows, LIBRARY_QC_COLUMNS)
    write_tsv(output_dir / "prep_metadata.tsv", resolved_prep_metadata, PREP_METADATA_COLUMNS)
    write_tsv(output_dir / "input_dna_quantity.tsv", input_dna_quantity, PREP_METADATA_COLUMNS)
    write_tsv(output_dir / "sample_aliases.tsv", alias_entries, SAMPLE_ALIASES_COLUMNS)
    write_tsv(
        output_dir / "sample_quality_master.tsv",
        sample_quality_master,
        SAMPLE_QUALITY_MASTER_COLUMNS,
    )
    write_tsv(output_dir / "dna_quality_scorecard.tsv", scorecard_rows, SCORECARD_COLUMNS)
    write_tsv(output_dir / "analysis_summary.tsv", analysis_summary_rows, ANALYSIS_SUMMARY_COLUMNS)
    write_report(
        output_dir / "report.md",
        output_run,
        sureselect_files,
        ddpcr_files,
        sample_files,
        resolved_library_rows,
        resolved_prep_metadata,
        scorecard_rows,
    )

    print(f"DNA quality outputs written to {output_dir}")
    return 0


# Keep the file import-safe so helpers can be exercised from ad hoc checks.
if __name__ == "__main__":
    sys.exit(main())
