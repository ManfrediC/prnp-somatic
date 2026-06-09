#!/usr/bin/env python3
"""Build the local ddPCR raw-file database and manifests."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import posixpath
import re
import shutil
import sys
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET


KEY_ENV_PATH = Path(__file__).resolve().parents[2] / "env" / "ddpcr_key.env"


def load_ddpcr_key() -> bytes:
    value = os.environ.get("DDPCR_KEY")
    if not value and KEY_ENV_PATH.exists():
        for line in KEY_ENV_PATH.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            name, separator, raw_value = line.partition("=")
            if separator and name.strip() == "DDPCR_KEY":
                value = raw_value.strip().strip('"').strip("'")
                break
    if not value:
        raise RuntimeError("DDPCR_KEY is required; set it in env/ddpcr_key.env or export it.")
    return value.encode("utf-8")


PASSWORD = load_ddpcr_key()
ASSAYS = ("D178N", "E200K", "P102L")
ASSAY_ORDER = {assay: index for index, assay in enumerate(ASSAYS)}
DATE_RE = re.compile(r"(20\d{2}-\d{2}-\d{2})")
CSV_RE = re.compile(r"^(20\d{2}-\d{2}-\d{2})\s*_SNV_(D1789?N|E200K|P102L)(?:_(.+))?\.csv$", re.IGNORECASE)
ASSAY_RE = re.compile(r"(D1789?N|E200K|P102L)", re.IGNORECASE)
XLSX_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
CONTROL_SAMPLE_RE = re.compile(r"^(NTC|WT|WT[-_ ]?HIGH|WT[-_ ]?COMMERCIAL|D178N[-_ ]?MUT|E200K[-_ ]?MUT|P102L[-_ ]?MUT)$", re.IGNORECASE)
LAYOUT_FROM_DDPCR_METADATA = {
    "2021-01-28": "Reviewed XLSX candidates were rejected; use .ddpcr plate metadata as sample-layout provenance.",
}

RUN_FIELDS = [
    "run_id", "run_date", "experiment", "assay", "status", "csv_source_kind",
    "csv_original_path", "csv_database_path", "csv_database_relative_path",
    "original_ddpcr_path", "renamed_ddpcr_path", "renamed_ddpcr_relative_path",
    "layout_original_path", "layout_database_path", "layout_database_relative_path",
    "archive_contents_dir", "archive_contents_relative_dir",
    "csv_rows", "csv_wells", "archive_wells_with_samples", "archive_wells_excluded",
    "layout_match_status", "layout_match_score", "notes",
]
FILE_FIELDS = [
    "run_id", "run_date", "experiment", "assay", "file_kind", "status",
    "original_path", "database_path", "database_relative_path",
    "original_name", "database_name", "sha256", "size_bytes", "modified_time",
    "duplicate_of", "notes",
]
SAMPLE_FIELDS = [
    "run_id", "run_date", "experiment", "assay", "well", "well_row", "well_col",
    "well_order", "target_order", "sample", "target", "target_clean", "target_role",
    "target_type", "dye", "channel", "source_csv", "source_ddpcr", "source_layout",
    "accepted_droplets", "positives", "negatives", "ch1_ch2_pos",
    "ch1_pos_ch2_neg", "ch1_neg_ch2_pos", "ch1_ch2_neg", "fractional_abundance",
    "poisson_fractional_abundance_min", "poisson_fractional_abundance_max",
]
ARCHIVE_CONTENT_FIELDS = [
    "run_id", "run_date", "experiment", "assay", "archive_member_path",
    "database_path", "database_relative_path", "content_kind", "well",
    "target_or_channel", "sha256", "size_bytes",
]
EXCLUDED_WELL_FIELDS = [
    "run_id", "run_date", "experiment", "assay", "well", "well_row", "well_col",
    "well_order", "sample_ids", "panel_targets", "reason", "source_ddpcr",
    "database_ddpcr", "notes",
]
SUPERSEDED_FIELDS = ["file_kind", "source_path", "superseded_by", "canonical_key", "reason", "sha256", "size_bytes"]
EXCLUDED_MATERIAL_FIELDS = ["file_kind", "category", "reason", "source_path", "inferred_date", "inferred_assay", "sha256", "size_bytes", "notes"]
LAYOUT_MATCH_FIELDS = [
    "run_date", "rank", "chosen", "status", "score", "source_path",
    "matched_informative_samples", "matched_common_samples",
    "missing_informative_samples", "assay_hits", "date_path_hit", "notes",
]


@dataclass(frozen=True)
class CsvSource:
    key: tuple[str, str]
    path: Path
    source_kind: str
    original_name: str
    canonical_name: str


@dataclass(frozen=True)
class ArchiveInfo:
    path: Path
    relative_path: str
    sha256: str
    size_bytes: int
    modified_time: str
    inferred_date: str
    inferred_assay: str
    assay_source: str
    plate_targets: tuple[str, ...]
    wells_with_samples: int
    archive_entries: int


@dataclass(frozen=True)
class LayoutChoice:
    run_date: str
    source_path: Path | None
    canonical_name: str
    score: float
    status: str
    matched_samples: tuple[str, ...]
    missing_samples: tuple[str, ...]
    notes: str


def normalise_assay(value: str | None) -> str:
    if not value:
        return ""
    match = ASSAY_RE.search(value)
    if not match:
        return ""
    return match.group(1).upper().replace("D1789N", "D178N")


def clean_target(value: str | None) -> str:
    if not value:
        return ""
    out = str(value).strip()
    out = re.sub(r"(-mut|_FAM1|_VIC2)$", "", out, flags=re.IGNORECASE)
    out = out.upper().replace("D1789N", "D178N")
    if out == "PRNP":
        return "WT"
    return out


def canonical_csv_name(run_date: str, assay: str) -> str:
    return f"{run_date}_SNV_{assay}.csv"


def canonical_ddpcr_name(run_date: str, assay: str) -> str:
    return f"{run_date}_SNV_{assay}.ddpcr"


def canonical_layout_name(path: Path, run_date: str) -> str:
    clean_stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", path.stem).strip("._")
    clean_stem = re.sub(r"-WS\d+(?:-\d+)?$", "", clean_stem)
    clean_stem = clean_stem or "layout"
    if clean_stem.startswith(run_date):
        return f"{clean_stem}{path.suffix.lower()}"
    return f"{run_date}_{clean_stem}{path.suffix.lower()}"


def layout_family_name(path: Path) -> str:
    stem = path.stem
    stem = re.sub(r"-WS\d+(?:-\d+)?$", "", stem, flags=re.IGNORECASE)
    stem = re.sub(r"\s*\(Automatisch gespeichert\)\s*", "", stem, flags=re.IGNORECASE)
    stem = re.sub(r"\s+", " ", stem).strip().lower()
    return stem


def parse_csv_key(path: Path) -> tuple[str, str, str] | None:
    match = CSV_RE.match(path.name)
    if not match:
        return None
    run_date, assay, suffix = match.groups()
    return run_date, normalise_assay(assay), suffix or ""


def infer_date(path: Path) -> str:
    match = DATE_RE.search(path.as_posix())
    return match.group(1) if match else ""


def date_variants(run_date: str) -> tuple[str, ...]:
    year, month, day = run_date.split("-")
    return (
        run_date,
        f"{day}-{month}-{year}",
        f"{day}_{month}_{year}",
        f"{day}.{month}.{year}",
        f"{year}-{day}-{month}",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def windows_path(path: Path | None) -> str:
    if path is None:
        return ""
    resolved = path.resolve()
    text = resolved.as_posix()
    if text.startswith("/mnt/c/"):
        return "C:\\" + text[len("/mnt/c/") :].replace("/", "\\")
    return str(resolved)


def path_from_manifest_value(value: str) -> Path:
    if value.startswith("C:\\") and Path("/mnt/c").exists():
        return Path("/mnt/c") / value[3:].replace("\\", "/")
    return Path(value)


def rel_db_path(path: Path | None, raw_root: Path) -> str:
    if path is None:
        return ""
    try:
        return path.resolve().relative_to(raw_root.resolve()).as_posix()
    except ValueError:
        return ""


def path_is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
    except ValueError:
        return False
    return True


def modified_time(path: Path) -> str:
    return str(int(path.stat().st_mtime))


def safe_member_path(destination: Path, member_name: str) -> Path:
    normalised = posixpath.normpath(member_name.replace("\\", "/"))
    if normalised.startswith("../") or normalised == ".." or posixpath.isabs(normalised):
        raise ValueError(f"Unsafe archive member path: {member_name}")
    target = (destination / Path(*normalised.split("/"))).resolve()
    destination_resolved = destination.resolve()
    if destination_resolved not in target.parents and target != destination_resolved:
        raise ValueError(f"Archive member escapes destination: {member_name}")
    return target


def read_ddplt_payload(archive_path: Path) -> dict:
    with zipfile.ZipFile(archive_path) as archive:
        ddplt_entries = [name for name in archive.namelist() if name.lower().endswith(".ddplt")]
        if not ddplt_entries:
            return {}
        payload = archive.read(ddplt_entries[0], pwd=PASSWORD)
    return json.loads(payload.decode("utf-8"))


def target_names_from_plate(plate: dict) -> tuple[str, ...]:
    targets: set[str] = set()
    for well in plate.get("WellSamples", []) or []:
        for target in (well.get("Panel") or {}).get("Targets", []) or []:
            cleaned = clean_target(target.get("TargetName"))
            if cleaned and cleaned != "WT":
                targets.add(cleaned)
    return tuple(sorted(targets))


def well_label_from_index(index: int) -> str:
    row_index = index // 12
    column_index = index % 12 + 1
    return f"{chr(ord('A') + row_index)}{column_index:02d}"


def well_parts(well: str) -> tuple[str, int, int]:
    match = re.match(r"^([A-H])(\d{1,2})$", well)
    if not match:
        return "", 0, 9999
    row = match.group(1)
    column = int(match.group(2))
    order = (ord(row) - ord("A")) * 12 + column
    return row, column, order


def sample_ids_from_well(well: dict) -> tuple[str, ...]:
    sample_ids = well.get("SampleIds") or []
    return tuple(str(item) for item in sample_ids if item is not None)


def archive_well_rows(archive_path: Path) -> list[dict[str, str]]:
    plate = read_ddplt_payload(archive_path)
    rows: list[dict[str, str]] = []
    for well in plate.get("WellSamples", []) or []:
        sample_ids = sample_ids_from_well(well)
        if not sample_ids:
            continue
        well_index = int(well.get("WellIndex", -1))
        well_label = well_label_from_index(well_index) if well_index >= 0 else ""
        targets = []
        channels = []
        dyes = []
        for target in (well.get("Panel") or {}).get("Targets", []) or []:
            targets.append(clean_target(target.get("TargetName")))
            dye = target.get("Dye") or {}
            channels.append(str(dye.get("Channel", "")))
            dyes.append(str(dye.get("DyeName", "")))
        row_letter, column, order = well_parts(well_label)
        rows.append(
            {
                "well": well_label,
                "well_row": row_letter,
                "well_col": str(column),
                "well_order": str(order),
                "sample_ids": ";".join(sample_ids),
                "panel_targets": ";".join(targets),
                "panel_channels": ";".join(channels),
                "panel_dyes": ";".join(dyes),
            }
        )
    return rows


def infer_archive_info(path: Path, archive_root: Path) -> ArchiveInfo:
    sha = sha256_file(path)
    relative_path = path.relative_to(archive_root).as_posix()
    inferred_date = infer_date(Path(relative_path)) or infer_date(path)
    filename_assay = normalise_assay(path.name)
    plate_targets: tuple[str, ...] = ()
    wells_with_samples = 0
    archive_entries = 0
    try:
        plate = read_ddplt_payload(path)
        plate_targets = target_names_from_plate(plate)
        wells_with_samples = sum(1 for well in plate.get("WellSamples", []) or [] if sample_ids_from_well(well))
        with zipfile.ZipFile(path) as archive:
            archive_entries = len(archive.namelist())
    except Exception:
        pass

    if filename_assay:
        inferred_assay = filename_assay
        assay_source = "filename"
    elif len(plate_targets) == 1 and plate_targets[0] in ASSAYS:
        inferred_assay = plate_targets[0]
        assay_source = "plate_targets"
    else:
        inferred_assay = ""
        assay_source = ""

    return ArchiveInfo(
        path=path,
        relative_path=relative_path,
        sha256=sha,
        size_bytes=path.stat().st_size,
        modified_time=modified_time(path),
        inferred_date=inferred_date,
        inferred_assay=inferred_assay,
        assay_source=assay_source,
        plate_targets=plate_targets,
        wells_with_samples=wells_with_samples,
        archive_entries=archive_entries,
    )


def existing_csv_source_kinds(runs_path: Path) -> dict[tuple[str, str], str]:
    rows = read_manifest_rows(runs_path)
    source_kinds: dict[tuple[str, str], str] = {}
    for row in rows:
        run_date = row.get("run_date", "").strip()
        assay = normalise_assay(row.get("assay", ""))
        source_kind = row.get("csv_source_kind", "").strip()
        if run_date and assay and source_kind:
            source_kinds[(run_date, assay)] = source_kind
    return source_kinds


def collect_active_csvs(
    csv_source_root: Path,
    corrected_csv: Path | None,
    source_kinds: dict[tuple[str, str], str] | None = None,
) -> tuple[dict[tuple[str, str], CsvSource], list[dict[str, str]]]:
    active: dict[tuple[str, str], CsvSource] = {}
    superseded: list[dict[str, str]] = []
    source_kinds = source_kinds or {}

    for path in sorted(csv_source_root.rglob("*.csv")):
        parsed = parse_csv_key(path)
        if not parsed:
            continue
        run_date, assay, suffix = parsed
        key = (run_date, assay)
        if suffix:
            continue
        active[key] = CsvSource(
            key=key,
            path=path,
            source_kind=source_kinds.get(key, "canonical_csv_export"),
            original_name=path.name,
            canonical_name=canonical_csv_name(run_date, assay),
        )

    corrected_key = ("2021-01-26", "E200K")
    if corrected_csv is not None and corrected_csv.exists():
        old_source = active.get(corrected_key)
        if old_source is not None:
            superseded.append(
                {
                    "file_kind": "csv_export",
                    "source_path": windows_path(old_source.path),
                    "superseded_by": windows_path(corrected_csv),
                    "canonical_key": "_".join(corrected_key),
                    "reason": "replaced_by_corrected_2021-01-26_E200K_export",
                    "sha256": sha256_file(old_source.path),
                    "size_bytes": str(old_source.path.stat().st_size),
                }
            )
        active[corrected_key] = CsvSource(
            key=corrected_key,
            path=corrected_csv,
            source_kind=source_kinds.get(corrected_key, "corrected_csv_override"),
            original_name=corrected_csv.name,
            canonical_name=canonical_csv_name(*corrected_key),
        )

    corrected_csv_search_root = corrected_csv.parent if corrected_csv is not None and corrected_csv.exists() else None
    if corrected_csv_search_root and corrected_csv_search_root.exists():
        for path in sorted(corrected_csv_search_root.glob("*.csv")):
            parsed = parse_csv_key(path)
            if not parsed:
                continue
            run_date, assay, suffix = parsed
            if suffix and "old" in suffix.lower():
                superseded.append(
                    {
                        "file_kind": "csv_export",
                        "source_path": windows_path(path),
                        "superseded_by": windows_path(corrected_csv) if (run_date, assay) == corrected_key else "",
                        "canonical_key": f"{run_date}_{assay}",
                        "reason": "filename_marks_superseded_old_export",
                        "sha256": sha256_file(path),
                        "size_bytes": str(path.stat().st_size),
                    }
                )

    return active, superseded


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [row for row in csv.DictReader(handle) if row.get("Well")]


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_manifest_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_run_checkpoint(
    checkpoint_root: Path,
    run_id: str,
    run_rows: list[dict[str, str]],
    file_rows: list[dict[str, str]],
    sample_rows: list[dict[str, str]],
    archive_content_rows: list[dict[str, str]],
    excluded_well_rows: list[dict[str, str]],
) -> None:
    run_dir = checkpoint_root / run_id
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    write_csv(run_dir / "runs.csv", run_rows, RUN_FIELDS)
    write_csv(run_dir / "files.csv", file_rows, FILE_FIELDS)
    write_csv(run_dir / "sample_manifest.csv", sample_rows, SAMPLE_FIELDS)
    write_csv(run_dir / "archive_contents.csv", archive_content_rows, ARCHIVE_CONTENT_FIELDS)
    write_csv(run_dir / "excluded_archive_wells.csv", excluded_well_rows, EXCLUDED_WELL_FIELDS)
    (run_dir / "_complete.txt").write_text("complete\n", encoding="utf-8")


def checkpointed_run_ids(checkpoint_root: Path) -> set[str]:
    if not checkpoint_root.exists():
        return set()
    return {
        child.name
        for child in checkpoint_root.iterdir()
        if child.is_dir() and (child / "_complete.txt").exists()
    }


def read_checkpoint_fragments(checkpoint_root: Path, filename: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not checkpoint_root.exists():
        return rows
    for run_dir in sorted(child for child in checkpoint_root.iterdir() if child.is_dir()):
        if not (run_dir / "_complete.txt").exists():
            continue
        rows.extend(read_manifest_rows(run_dir / filename))
    return rows


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() == destination.resolve():
        return
    shutil.copy2(source, destination)


def extract_archive(source: Path, destination: Path) -> list[dict[str, str]]:
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)

    extracted: list[dict[str, str]] = []
    with zipfile.ZipFile(source) as archive:
        for member in archive.infolist():
            target = safe_member_path(destination, member.filename)
            if member.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(member, pwd=PASSWORD) as opened, target.open("wb") as output:
                shutil.copyfileobj(opened, output)
            extracted.append({"member": member.filename, "path": str(target)})
    return extracted


def content_kind(member_name: str) -> tuple[str, str, str]:
    member_path = Path(member_name.replace("\\", "/"))
    suffix = member_path.suffix.lower()
    stem = member_path.stem
    well = ""
    target_or_channel = ""
    if suffix == ".ddplt":
        return "plate_definition", well, target_or_channel
    if suffix == ".ddpeakjson":
        well = stem.split("_")[0]
        target_or_channel = "_".join(stem.split("_")[1:])
        return "droplet_amplitudes", well, target_or_channel
    if suffix == ".ddmetajson":
        well = stem
        return "droplet_classification_and_gates", well, target_or_channel
    if suffix == ".json":
        return "json", well, target_or_channel
    return suffix.lstrip(".") or "file", well, target_or_channel


def normalise_text(value: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "", str(value).upper())


def xlsx_strings(path: Path) -> set[str]:
    values: set[str] = set()
    try:
        with zipfile.ZipFile(path) as archive:
            shared_strings: list[str] = []
            try:
                root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
                for item in root.findall(XLSX_NS + "si"):
                    text = "".join(node.text or "" for node in item.iter(XLSX_NS + "t"))
                    shared_strings.append(text)
            except KeyError:
                pass

            worksheet_names = [
                name for name in archive.namelist()
                if name.startswith("xl/worksheets/") and name.endswith(".xml")
            ]
            for worksheet_name in worksheet_names:
                root = ET.fromstring(archive.read(worksheet_name))
                for cell in root.iter(XLSX_NS + "c"):
                    cell_type = cell.attrib.get("t")
                    raw_value = cell.find(XLSX_NS + "v")
                    value = ""
                    if cell_type == "s" and raw_value is not None and raw_value.text is not None:
                        index = int(raw_value.text)
                        if index < len(shared_strings):
                            value = shared_strings[index]
                    elif cell_type == "inlineStr":
                        value = "".join(node.text or "" for node in cell.iter(XLSX_NS + "t"))
                    elif raw_value is not None and raw_value.text is not None:
                        value = raw_value.text
                    if value:
                        values.add(normalise_text(value))
    except (KeyError, zipfile.BadZipFile, ET.ParseError):
        return set()
    return values


def collect_layout_candidates(usz_ddpcr_root: Path) -> list[Path]:
    runs_root = usz_ddpcr_root / "Runs"
    if not runs_root.exists():
        return []
    candidates = [
        path for path in runs_root.rglob("*")
        if path.is_file() and path.suffix.lower() in {".xlsx", ".xlsm"}
    ]
    return sorted(candidates)


def collect_local_layout_candidates(layout_root: Path) -> list[Path]:
    if not layout_root.exists():
        return []
    return sorted(
        path for path in layout_root.rglob("*")
        if path.is_file() and path.suffix.lower() in {".xlsx", ".xlsm"}
    )


def score_layouts(
    run_date: str,
    assays: set[str],
    samples: set[str],
    candidates: list[Path],
    workbook_string_cache: dict[Path, set[str]],
) -> tuple[LayoutChoice, list[dict[str, str]]]:
    informative_samples = sorted(
        sample for sample in samples
        if sample and not CONTROL_SAMPLE_RE.match(sample) and normalise_text(sample) not in {"", "WT"}
    )
    all_samples = sorted(sample for sample in samples if sample)
    scored_all: list[dict[str, object]] = []
    run_date_variants = date_variants(run_date)

    for candidate in candidates:
        if candidate not in workbook_string_cache:
            workbook_string_cache[candidate] = xlsx_strings(candidate)
        strings = workbook_string_cache[candidate]
        if not strings:
            continue
        matched = [
            sample for sample in informative_samples
            if normalise_text(sample) in strings
        ]
        matched_common = [
            sample for sample in all_samples
            if sample not in matched and normalise_text(sample) in strings
        ]
        assay_hits = sum(1 for assay in assays if normalise_text(assay) in strings)
        path_text = candidate.as_posix()
        date_path_hit = 1 if any(variant in path_text for variant in run_date_variants) else 0
        name_penalty = 0.0
        if re.search(r"-WS\d+", candidate.name, flags=re.IGNORECASE):
            name_penalty += 0.25
        if "automatisch gespeichert" in candidate.name.lower():
            name_penalty += 2.0
        if re.search(r"\bkopie\b|\bcopy\b", candidate.name, flags=re.IGNORECASE):
            name_penalty += 2.0
        score = (
            12.0 * len(matched)
            + 2.0 * len(matched_common)
            + 5.0 * assay_hits
            + 75.0 * date_path_hit
            - name_penalty
        )
        if score <= 0:
            continue
        scored_all.append(
            {
                "source_path": candidate,
                "score": score,
                "matched": tuple(matched),
                "matched_common": tuple(matched_common),
                "assay_hits": assay_hits,
                "date_path_hit": date_path_hit,
                "name_penalty": name_penalty,
            }
        )

    scored = [item for item in scored_all if int(item["date_path_hit"]) == 1] or scored_all
    scored.sort(key=lambda item: (-float(item["score"]), windows_path(item["source_path"])))
    top_rows: list[dict[str, str]] = []
    for rank, item in enumerate(scored[:10], start=1):
        source_path = item["source_path"]
        assert isinstance(source_path, Path)
        matched = tuple(item["matched"])
        matched_common = tuple(item["matched_common"])
        missing = tuple(sample for sample in informative_samples if sample not in matched)
        top_rows.append(
            {
                "run_date": run_date,
                "rank": str(rank),
                "chosen": "false",
                "status": "candidate",
                "score": f"{float(item['score']):.2f}",
                "source_path": windows_path(source_path),
                "matched_informative_samples": ";".join(matched),
                "matched_common_samples": ";".join(matched_common),
                "missing_informative_samples": ";".join(missing[:40]),
                "assay_hits": str(item["assay_hits"]),
                "date_path_hit": str(item["date_path_hit"]),
                "notes": "",
            }
        )

    if not scored:
        return (
            LayoutChoice(
                run_date=run_date,
                source_path=None,
                canonical_name="",
                score=0.0,
                status="no_match",
                matched_samples=(),
                missing_samples=tuple(informative_samples),
                notes="No XLSX candidate had assay/sample/date evidence.",
            ),
            top_rows,
        )

    best = scored[0]
    second_score = float(scored[1]["score"]) if len(scored) > 1 else 0.0
    best_path = best["source_path"]
    assert isinstance(best_path, Path)
    best_matched = tuple(best["matched"])
    best_missing = tuple(sample for sample in informative_samples if sample not in best_matched)
    status = "matched"
    notes = ""
    if len(best_matched) == 0 and not bool(best["date_path_hit"]):
        status = "weak_match"
        notes = "Chosen by weak evidence; review before relying on layout."
    elif (
        second_score
        and float(best["score"]) - second_score < 5
        and layout_family_name(best_path) != layout_family_name(scored[1]["source_path"])
    ):
        status = "ambiguous"
        notes = "Top two layout candidates have similar scores."

    if top_rows:
        top_rows[0]["chosen"] = "true"
        top_rows[0]["status"] = status
        top_rows[0]["notes"] = notes

    return (
        LayoutChoice(
            run_date=run_date,
            source_path=best_path,
            canonical_name=canonical_layout_name(best_path, run_date),
            score=float(best["score"]),
            status=status,
            matched_samples=best_matched,
            missing_samples=best_missing,
            notes=notes,
        ),
        top_rows,
    )


def classify_excluded(path: Path) -> str:
    text = path.as_posix().lower()
    if "analysis_testing" in text or "analysis test" in text:
        return "analysis_test"
    if "failed" in text:
        return "failed"
    if "cnv" in text:
        return "cnv"
    if "lod" in text:
        return "lod"
    if "gradient" in text or "temperature" in text:
        return "gradient"
    if "old" in text or "superseded" in text:
        return "old_or_superseded"
    if "snv" in text:
        return "inactive_snv"
    return "other_ddpcr_related"


def choose_archives(active_keys: set[tuple[str, str]], archive_infos: list[ArchiveInfo]) -> tuple[dict[tuple[str, str], ArchiveInfo], list[dict[str, str]]]:
    by_key: dict[tuple[str, str], list[ArchiveInfo]] = defaultdict(list)
    for info in archive_infos:
        if info.inferred_date and info.inferred_assay:
            by_key[(info.inferred_date, info.inferred_assay)].append(info)

    chosen: dict[tuple[str, str], ArchiveInfo] = {}
    excluded: list[dict[str, str]] = []
    selected_paths: set[Path] = set()

    for key in sorted(active_keys, key=lambda item: (item[0], ASSAY_ORDER.get(item[1], 99))):
        candidates = by_key.get(key, [])
        if not candidates:
            continue
        candidates = sorted(
            candidates,
            key=lambda info: (
                0 if info.assay_source == "filename" else 1,
                0 if key[0] in info.relative_path else 1,
                info.relative_path,
            ),
        )
        chosen[key] = candidates[0]
        selected_paths.add(candidates[0].path)
        for item in candidates[1:]:
            excluded.append(
                {
                    "file_kind": "ddpcr_archive",
                    "category": "unselected_duplicate_or_alternative",
                    "reason": "another_archive_selected_for_active_run_key",
                    "source_path": windows_path(item.path),
                    "inferred_date": item.inferred_date,
                    "inferred_assay": item.inferred_assay,
                    "sha256": item.sha256,
                    "size_bytes": str(item.size_bytes),
                    "notes": f"selected={windows_path(candidates[0].path)}",
                }
            )

    for info in archive_infos:
        if info.path in selected_paths:
            continue
        key = (info.inferred_date, info.inferred_assay)
        if key in active_keys:
            continue
        reason = "no_active_snv_csv_for_archive_key" if info.inferred_assay in ASSAYS else "not_an_active_snv_assay_archive"
        excluded.append(
            {
                "file_kind": "ddpcr_archive",
                "category": classify_excluded(info.path),
                "reason": reason,
                "source_path": windows_path(info.path),
                "inferred_date": info.inferred_date,
                "inferred_assay": info.inferred_assay,
                "sha256": info.sha256,
                "size_bytes": str(info.size_bytes),
                "notes": f"plate_targets={';'.join(info.plate_targets)}",
            }
        )

    return chosen, excluded


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument(
        "--use-local-archives",
        action="store_true",
        help=(
            "Use raw/ddpcr/ddpcr_archive as the .ddpcr source and preserve "
            "the local archive, CSV, and layout directories while rebuilding "
            "archive_contents and manifests."
        ),
    )
    parser.add_argument(
        "--usz-ddpcr-root",
        type=Path,
        default=Path("/mnt/c/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/CJD PRNP/Experiments/ddPCR"),
    )
    parser.add_argument(
        "--csv-source-root",
        type=Path,
        default=None,
        help="Directory containing reviewed CSV exports; defaults to raw/ddpcr/csv_export.",
    )
    parser.add_argument(
        "--corrected-csv",
        type=Path,
        default=None,
        help="Optional corrected CSV override for 2021-01-26 E200K.",
    )
    parser.add_argument("--raw-root", type=Path, default=None)
    parser.add_argument("--resume", action="store_true", help="Keep completed per-run checkpoints and skip them.")
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    raw_root = (args.raw_root or repo_root / "raw/ddpcr").resolve()
    usz_ddpcr_root = args.usz_ddpcr_root.resolve()

    expected_raw_root = (repo_root / "raw/ddpcr").resolve()
    if raw_root != expected_raw_root:
        print(f"Refusing to write outside raw/ddpcr: {raw_root}", file=sys.stderr)
        return 2

    manifests_root = raw_root / "manifests"
    ddpcr_archive_root = raw_root / "ddpcr_archive"
    archive_contents_root = raw_root / "archive_contents"
    csv_export_root = raw_root / "csv_export"
    layout_root = raw_root / "layout_xlsx"
    checkpoint_root = manifests_root / "_by_run"
    csv_source_root = (args.csv_source_root or csv_export_root).resolve()
    corrected_csv = args.corrected_csv.resolve() if args.corrected_csv else None
    csv_source_kinds = existing_csv_source_kinds(manifests_root / "runs.csv")
    archive_source_root = ddpcr_archive_root if args.use_local_archives else usz_ddpcr_root / "ddPCR-files"

    if not args.resume:
        managed_dirs = [manifests_root, archive_contents_root]
        if not args.use_local_archives:
            managed_dirs.extend([ddpcr_archive_root, layout_root])
        if not path_is_relative_to(csv_source_root, csv_export_root):
            managed_dirs.append(csv_export_root)
        for managed_dir in managed_dirs:
            if managed_dir.exists():
                shutil.rmtree(managed_dir)
    completed_run_ids = checkpointed_run_ids(checkpoint_root) if args.resume else set()

    print("Collecting active CSV exports...", flush=True)
    active_csvs, superseded_rows = collect_active_csvs(csv_source_root, corrected_csv, csv_source_kinds)
    active_keys = set(active_csvs)

    print("Inspecting .ddpcr archives...", flush=True)
    archive_infos = [
        infer_archive_info(path, archive_source_root)
        for path in sorted(archive_source_root.rglob("*.ddpcr"))
    ]
    chosen_archives, excluded_rows = choose_archives(active_keys, archive_infos)

    csv_rows_by_key: dict[tuple[str, str], list[dict[str, str]]] = {
        key: read_csv_rows(source.path)
        for key, source in active_csvs.items()
    }

    date_samples: dict[str, set[str]] = defaultdict(set)
    date_assays: dict[str, set[str]] = defaultdict(set)
    for (run_date, assay), rows in csv_rows_by_key.items():
        date_assays[run_date].add(assay)
        for row in rows:
            sample = row.get("Sample", "").strip()
            if sample:
                date_samples[run_date].add(sample)

    print("Scoring sample-layout workbooks...", flush=True)
    layout_candidates = (
        collect_local_layout_candidates(layout_root)
        if args.use_local_archives
        else collect_layout_candidates(usz_ddpcr_root)
    )
    layout_choices: dict[str, LayoutChoice] = {}
    layout_match_rows: list[dict[str, str]] = []
    workbook_string_cache: dict[Path, set[str]] = {}
    for run_date in sorted(date_samples):
        choice, rows = score_layouts(
            run_date,
            date_assays[run_date],
            date_samples[run_date],
            layout_candidates,
            workbook_string_cache,
        )
        if run_date in LAYOUT_FROM_DDPCR_METADATA:
            rows.insert(
                0,
                {
                    "run_date": run_date,
                    "rank": "1",
                    "chosen": "true",
                    "status": "resolved_from_ddpcr_metadata",
                    "score": "0.00",
                    "source_path": "",
                    "matched_informative_samples": "",
                    "matched_common_samples": "",
                    "missing_informative_samples": "",
                    "assay_hits": "0",
                    "date_path_hit": "0",
                    "notes": LAYOUT_FROM_DDPCR_METADATA[run_date],
                },
            )
            for rank, row in enumerate(rows[1:], start=2):
                row["rank"] = str(rank)
                row["chosen"] = "false"
                row["status"] = "rejected_candidate"
                row["notes"] = "Reviewed and rejected as active sample-layout provenance."
            choice = LayoutChoice(
                run_date=run_date,
                source_path=None,
                canonical_name="",
                score=0.0,
                status="resolved_from_ddpcr_metadata",
                matched_samples=(),
                missing_samples=(),
                notes=LAYOUT_FROM_DDPCR_METADATA[run_date],
            )
        layout_choices[run_date] = choice
        layout_match_rows.extend(rows)

    runs_rows: list[dict[str, str]] = []
    file_rows: list[dict[str, str]] = []
    sample_rows: list[dict[str, str]] = []
    archive_content_rows: list[dict[str, str]] = []
    excluded_well_rows: list[dict[str, str]] = []

    copied_layouts: dict[Path, Path] = {}

    for key in sorted(active_csvs, key=lambda item: (item[0], ASSAY_ORDER.get(item[1], 99))):
        run_date, assay = key
        run_id = f"{run_date}_SNV_{assay}"
        if run_id in completed_run_ids:
            print(f"Skipping checkpointed {run_id}...", flush=True)
            continue
        print(f"Building {run_id}...", flush=True)
        file_start = len(file_rows)
        sample_start = len(sample_rows)
        archive_content_start = len(archive_content_rows)
        excluded_well_start = len(excluded_well_rows)
        csv_source = active_csvs[key]
        csv_destination = csv_export_root / run_date / csv_source.canonical_name
        copy_file(csv_source.path, csv_destination)

        archive_info = chosen_archives.get(key)
        archive_destination: Path | None = None
        contents_destination: Path | None = None
        raw_well_rows: list[dict[str, str]] = []
        if archive_info is not None:
            archive_destination = ddpcr_archive_root / run_date / canonical_ddpcr_name(run_date, assay)
            copy_file(archive_info.path, archive_destination)
            contents_destination = archive_contents_root / run_date / f"SNV_{assay}"
            extracted = extract_archive(archive_info.path, contents_destination)
            raw_well_rows = archive_well_rows(archive_info.path)
            for item in extracted:
                extracted_path = Path(item["path"])
                kind, well, target_or_channel = content_kind(item["member"])
                archive_content_rows.append(
                    {
                        "run_id": run_id,
                        "run_date": run_date,
                        "experiment": "SNV",
                        "assay": assay,
                        "archive_member_path": item["member"],
                        "database_path": windows_path(extracted_path),
                        "database_relative_path": rel_db_path(extracted_path, raw_root),
                        "content_kind": kind,
                        "well": well,
                        "target_or_channel": target_or_channel,
                        "sha256": sha256_file(extracted_path),
                        "size_bytes": str(extracted_path.stat().st_size),
                    }
                )

        layout_choice = layout_choices.get(run_date)
        layout_destination: Path | None = None
        if layout_choice and layout_choice.source_path is not None:
            layout_destination = layout_root / run_date / layout_choice.canonical_name
            if layout_choice.source_path not in copied_layouts:
                copy_file(layout_choice.source_path, layout_destination)
                copied_layouts[layout_choice.source_path] = layout_destination
            else:
                layout_destination = copied_layouts[layout_choice.source_path]

        active_rows = csv_rows_by_key[key]
        active_wells = {row.get("Well", "") for row in active_rows if row.get("Well")}
        raw_wells = {row["well"] for row in raw_well_rows if row.get("well")}
        excluded_raw_wells = sorted(raw_wells - active_wells, key=lambda well: well_parts(well)[2])
        raw_well_by_name = {row["well"]: row for row in raw_well_rows}
        for well in excluded_raw_wells:
            well_row = raw_well_by_name[well]
            excluded_well_rows.append(
                {
                    "run_id": run_id,
                    "run_date": run_date,
                    "experiment": "SNV",
                    "assay": assay,
                    "well": well,
                    "well_row": well_row["well_row"],
                    "well_col": well_row["well_col"],
                    "well_order": well_row["well_order"],
                    "sample_ids": well_row["sample_ids"],
                    "panel_targets": well_row["panel_targets"],
                    "reason": "present_in_archive_not_active_csv_export",
                    "source_ddpcr": windows_path(archive_info.path) if archive_info else "",
                    "database_ddpcr": windows_path(archive_destination),
                    "notes": "Excluded from sample_manifest to avoid importing failed or non-exported wells.",
                }
            )

        runs_rows.append(
            {
                "run_id": run_id,
                "run_date": run_date,
                "experiment": "SNV",
                "assay": assay,
                "status": "active",
                "csv_source_kind": csv_source.source_kind,
                "csv_original_path": windows_path(csv_source.path),
                "csv_database_path": windows_path(csv_destination),
                "csv_database_relative_path": rel_db_path(csv_destination, raw_root),
                "original_ddpcr_path": windows_path(archive_info.path) if archive_info else "",
                "renamed_ddpcr_path": windows_path(archive_destination),
                "renamed_ddpcr_relative_path": rel_db_path(archive_destination, raw_root),
                "layout_original_path": windows_path(layout_choice.source_path) if layout_choice else "",
                "layout_database_path": windows_path(layout_destination),
                "layout_database_relative_path": rel_db_path(layout_destination, raw_root),
                "archive_contents_dir": windows_path(contents_destination),
                "archive_contents_relative_dir": rel_db_path(contents_destination, raw_root),
                "csv_rows": str(len(active_rows)),
                "csv_wells": str(len(active_wells)),
                "archive_wells_with_samples": str(len(raw_wells)),
                "archive_wells_excluded": str(len(excluded_raw_wells)),
                "layout_match_status": layout_choice.status if layout_choice else "no_match",
                "layout_match_score": f"{layout_choice.score:.2f}" if layout_choice else "0.00",
                "notes": "" if archive_info else "No USZ .ddpcr archive matched this active CSV key.",
            }
        )

        if archive_destination is not None:
            file_rows.append(
                {
                    "run_id": run_id,
                    "run_date": run_date,
                    "experiment": "SNV",
                    "assay": assay,
                    "file_kind": "ddpcr_archive",
                    "status": "active",
                    "original_path": windows_path(archive_info.path) if archive_info else "",
                    "database_path": windows_path(archive_destination),
                    "database_relative_path": rel_db_path(archive_destination, raw_root),
                    "original_name": archive_info.path.name if archive_info else "",
                    "database_name": archive_destination.name,
                    "sha256": sha256_file(archive_destination),
                    "size_bytes": str(archive_destination.stat().st_size),
                    "modified_time": modified_time(archive_destination),
                    "duplicate_of": "",
                    "notes": "",
                }
            )

        file_rows.append(
            {
                "run_id": run_id,
                "run_date": run_date,
                "experiment": "SNV",
                "assay": assay,
                "file_kind": "csv_export",
                "status": "active",
                "original_path": windows_path(csv_source.path),
                "database_path": windows_path(csv_destination),
                "database_relative_path": rel_db_path(csv_destination, raw_root),
                "original_name": csv_source.original_name,
                "database_name": csv_destination.name,
                "sha256": sha256_file(csv_destination),
                "size_bytes": str(csv_destination.stat().st_size),
                "modified_time": modified_time(csv_destination),
                "duplicate_of": "",
                "notes": csv_source.source_kind,
            }
        )

        if layout_destination is not None and layout_choice and layout_choice.source_path is not None:
            file_rows.append(
                {
                    "run_id": run_id,
                    "run_date": run_date,
                    "experiment": "SNV",
                    "assay": assay,
                    "file_kind": "sample_layout_xlsx",
                    "status": layout_choice.status,
                    "original_path": windows_path(layout_choice.source_path),
                    "database_path": windows_path(layout_destination),
                    "database_relative_path": rel_db_path(layout_destination, raw_root),
                    "original_name": layout_choice.source_path.name,
                    "database_name": layout_destination.name,
                    "sha256": sha256_file(layout_destination),
                    "size_bytes": str(layout_destination.stat().st_size),
                    "modified_time": modified_time(layout_destination),
                    "duplicate_of": "",
                    "notes": layout_choice.notes,
                }
            )

        for row_index, row in enumerate(active_rows, start=1):
            well = row.get("Well", "")
            well_row_letter, well_col, well_order = well_parts(well)
            target = row.get("Target", "")
            target_clean = clean_target(target)
            target_role = "reference" if target_clean == "WT" else "mutant" if target_clean == assay else "other"
            dye = row.get("DyeName(s)", "")
            channel = "Ch1" if dye.upper() == "FAM" else "Ch2" if dye.upper() in {"VIC", "HEX"} else ""
            sample_rows.append(
                {
                    "run_id": run_id,
                    "run_date": run_date,
                    "experiment": "SNV",
                    "assay": assay,
                    "well": well,
                    "well_row": well_row_letter,
                    "well_col": str(well_col),
                    "well_order": str(well_order),
                    "target_order": str(row_index),
                    "sample": row.get("Sample", ""),
                    "target": target,
                    "target_clean": target_clean,
                    "target_role": target_role,
                    "target_type": row.get("TargetType", ""),
                    "dye": dye,
                    "channel": channel,
                    "source_csv": windows_path(csv_destination),
                    "source_ddpcr": windows_path(archive_destination),
                    "source_layout": windows_path(layout_destination),
                    "accepted_droplets": row.get("Accepted Droplets", ""),
                    "positives": row.get("Positives", ""),
                    "negatives": row.get("Negatives", ""),
                    "ch1_ch2_pos": row.get("Ch1+Ch2+", ""),
                    "ch1_pos_ch2_neg": row.get("Ch1+Ch2-", ""),
                    "ch1_neg_ch2_pos": row.get("Ch1-Ch2+", ""),
                    "ch1_ch2_neg": row.get("Ch1-Ch2-", ""),
                    "fractional_abundance": row.get("Fractional Abundance", ""),
                    "poisson_fractional_abundance_min": row.get("PoissonFractionalAbundanceMin", ""),
                    "poisson_fractional_abundance_max": row.get("PoissonFractionalAbundanceMax", ""),
                }
            )

        write_run_checkpoint(
            checkpoint_root,
            run_id,
            runs_rows[-1:],
            file_rows[file_start:],
            sample_rows[sample_start:],
            archive_content_rows[archive_content_start:],
            excluded_well_rows[excluded_well_start:],
        )

    selected_layout_sources = {choice.source_path for choice in layout_choices.values() if choice.source_path is not None}
    for candidate in layout_candidates:
        if candidate in selected_layout_sources:
            continue
        excluded_rows.append(
            {
                "file_kind": "sample_layout_xlsx",
                "category": classify_excluded(candidate),
                "reason": "not_selected_as_active_snv_layout",
                "source_path": windows_path(candidate),
                "inferred_date": infer_date(candidate),
                "inferred_assay": normalise_assay(candidate.name),
                "sha256": sha256_file(candidate),
                "size_bytes": str(candidate.stat().st_size),
                "notes": "",
            }
        )

    legacy_csv_search_root = repo_root / "legacy/scratch/ddpcr/2026-06-finalisation/ddpcr_raw_import/ddpcr-files"
    active_csv_paths = {source.path.resolve() for source in active_csvs.values()}
    superseded_paths = {path_from_manifest_value(row["source_path"]).resolve() for row in superseded_rows if row["source_path"]}
    for csv_root in [csv_source_root, legacy_csv_search_root]:
        if not csv_root.exists():
            continue
        for path in sorted(csv_root.rglob("*.csv")):
            resolved = path.resolve()
            if resolved in active_csv_paths or resolved in superseded_paths:
                continue
            parsed = parse_csv_key(path)
            if parsed:
                run_date, assay, suffix = parsed
            else:
                run_date, assay, suffix = infer_date(path), normalise_assay(path.name), ""
            excluded_rows.append(
                {
                    "file_kind": "csv_export",
                    "category": classify_excluded(path),
                    "reason": "not_active_snv_workflow_csv_export",
                    "source_path": windows_path(path),
                    "inferred_date": run_date,
                    "inferred_assay": assay,
                    "sha256": sha256_file(path),
                    "size_bytes": str(path.stat().st_size),
                    "notes": f"suffix={suffix}" if suffix else "",
                }
            )

    checkpointed_rows = read_checkpoint_fragments(checkpoint_root, "runs.csv")
    if checkpointed_rows:
        runs_rows = checkpointed_rows
        file_rows = read_checkpoint_fragments(checkpoint_root, "files.csv")
        sample_rows = read_checkpoint_fragments(checkpoint_root, "sample_manifest.csv")
        archive_content_rows = read_checkpoint_fragments(checkpoint_root, "archive_contents.csv")
        excluded_well_rows = read_checkpoint_fragments(checkpoint_root, "excluded_archive_wells.csv")

    sample_rows.sort(
        key=lambda row: (
            row["run_date"],
            ASSAY_ORDER.get(row["assay"], 99),
            int(row["well_order"] or "9999"),
            int(row["target_order"] or "9999"),
        )
    )
    archive_content_rows.sort(key=lambda row: (row["run_date"], ASSAY_ORDER.get(row["assay"], 99), row["archive_member_path"]))
    file_rows.sort(key=lambda row: (row["run_date"], ASSAY_ORDER.get(row["assay"], 99), row["file_kind"], row["database_name"]))
    runs_rows.sort(key=lambda row: (row["run_date"], ASSAY_ORDER.get(row["assay"], 99)))
    excluded_rows.sort(key=lambda row: (row["category"], row["file_kind"], row["source_path"]))
    excluded_well_rows.sort(key=lambda row: (row["run_date"], ASSAY_ORDER.get(row["assay"], 99), int(row["well_order"] or "9999")))
    superseded_rows.sort(key=lambda row: row["source_path"])

    write_csv(
        manifests_root / "runs.csv",
        runs_rows,
        [
            "run_id", "run_date", "experiment", "assay", "status", "csv_source_kind",
            "csv_original_path", "csv_database_path", "csv_database_relative_path",
            "original_ddpcr_path", "renamed_ddpcr_path", "renamed_ddpcr_relative_path",
            "layout_original_path", "layout_database_path", "layout_database_relative_path",
            "archive_contents_dir", "archive_contents_relative_dir",
            "csv_rows", "csv_wells", "archive_wells_with_samples", "archive_wells_excluded",
            "layout_match_status", "layout_match_score", "notes",
        ],
    )
    write_csv(
        manifests_root / "files.csv",
        file_rows,
        [
            "run_id", "run_date", "experiment", "assay", "file_kind", "status",
            "original_path", "database_path", "database_relative_path",
            "original_name", "database_name", "sha256", "size_bytes", "modified_time",
            "duplicate_of", "notes",
        ],
    )
    write_csv(
        manifests_root / "sample_manifest.csv",
        sample_rows,
        [
            "run_id", "run_date", "experiment", "assay", "well", "well_row", "well_col",
            "well_order", "target_order", "sample", "target", "target_clean", "target_role",
            "target_type", "dye", "channel", "source_csv", "source_ddpcr", "source_layout",
            "accepted_droplets", "positives", "negatives", "ch1_ch2_pos",
            "ch1_pos_ch2_neg", "ch1_neg_ch2_pos", "ch1_ch2_neg", "fractional_abundance",
            "poisson_fractional_abundance_min", "poisson_fractional_abundance_max",
        ],
    )
    write_csv(
        manifests_root / "archive_contents.csv",
        archive_content_rows,
        [
            "run_id", "run_date", "experiment", "assay", "archive_member_path",
            "database_path", "database_relative_path", "content_kind", "well",
            "target_or_channel", "sha256", "size_bytes",
        ],
    )
    write_csv(
        manifests_root / "superseded_files.csv",
        superseded_rows,
        ["file_kind", "source_path", "superseded_by", "canonical_key", "reason", "sha256", "size_bytes"],
    )
    write_csv(
        manifests_root / "excluded_material.csv",
        excluded_rows,
        ["file_kind", "category", "reason", "source_path", "inferred_date", "inferred_assay", "sha256", "size_bytes", "notes"],
    )
    write_csv(
        manifests_root / "excluded_archive_wells.csv",
        excluded_well_rows,
        [
            "run_id", "run_date", "experiment", "assay", "well", "well_row", "well_col",
            "well_order", "sample_ids", "panel_targets", "reason", "source_ddpcr",
            "database_ddpcr", "notes",
        ],
    )
    write_csv(
        manifests_root / "layout_matching.csv",
        layout_match_rows,
        [
            "run_date", "rank", "chosen", "status", "score", "source_path",
            "matched_informative_samples", "matched_common_samples",
            "missing_informative_samples", "assay_hits", "date_path_hit", "notes",
        ],
    )

    active_missing_archives = [
        f"{run_date}_SNV_{assay}"
        for run_date, assay in sorted(active_keys)
        if (run_date, assay) not in chosen_archives
    ]
    layout_status_counts: defaultdict[str, int] = defaultdict(int)
    for row in runs_rows:
        layout_status_counts[row["layout_match_status"]] += 1
    summary_lines = [
        "ddPCR raw database build summary",
        f"raw_root: {windows_path(raw_root)}",
        f"active_runs: {len(runs_rows)}",
        f"active_csv_exports: {len(active_csvs)}",
        f"active_ddpcr_archives: {sum(1 for row in runs_rows if row['original_ddpcr_path'])}",
        f"sample_manifest_rows: {len(sample_rows)}",
        f"archive_content_rows: {len(archive_content_rows)}",
        f"excluded_archive_wells: {len(excluded_well_rows)}",
        f"superseded_files: {len(superseded_rows)}",
        f"excluded_material_files: {len(excluded_rows)}",
        "layout_status_counts: " + "; ".join(f"{key}={value}" for key, value in sorted(layout_status_counts.items())),
        "missing_active_archives: " + ("; ".join(active_missing_archives) if active_missing_archives else "none"),
    ]
    (manifests_root / "run_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("\n".join(summary_lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
