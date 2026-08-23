#!/usr/bin/env python3
"""Check the manuscript figure/table crosswalk, wrappers and bundle inputs."""

from __future__ import annotations

import argparse
import csv
import contextlib
import io
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from generate_caption_list import OUTPUT_PATH as CAPTION_CATALOGUE_PATH
from generate_caption_list import build_caption_catalogue


FIGURES = [f"Figure {i}" for i in range(1, 9)] + [f"Supplementary Figure S{i}" for i in range(1, 5)]
TABLES = [f"Table {i}" for i in range(1, 8)] + [f"Supplementary Table S{i}" for i in range(1, 6)]

MAIN_FIGURE_SOURCES = [
    "manuscript/figures/main/figure_mosaicism_models.tex",
    "manuscript/figures/main/figure_dna_purification_protocol.tex",
    "manuscript/figures/main/figure_snv_detection_workflow.tex",
    "manuscript/figures/main/figure_rt_quic_decontamination.tex",
    "manuscript/figures/main/figure_ddpcr_mutant_allele_frequencies.tex",
    "manuscript/figures/main/figure_ddpcr_e200k_positive_wells.tex",
    "manuscript/figures/main/figure_ddpcr_e200k_positive_well_vaf.tex",
    "manuscript/figures/main/figure_snv_lollipop_variants.tex",
]
SUPPLEMENT_FIGURE_SOURCES = [
    "manuscript/figures/supplement/figure_control4_gdna_screentape.tex",
    "manuscript/figures/supplement/figure_ddpcr_gating_strategy.tex",
    "manuscript/figures/supplement/figure_ddpcr_pooled_participant_analysis.tex",
    "manuscript/figures/supplement/figure_ddpcr_lod_summary.tex",
]
MAIN_TABLE_SOURCES = [
    "manuscript/tables/main/table_patient_cohort_overview.tex",
    "manuscript/tables/main/table_ddpcr_lod_quantification.tex",
    "manuscript/tables/ddpcr_sample_number/ddpcr_sample_number.tex",
    "manuscript/tables/main/table_a117v_spike_in_lod.tex",
    "manuscript/tables/main/table_sequencing_metrics.tex",
    "manuscript/tables/main/table_prnp_somatic_snv_summary.tex",
    "manuscript/tables/main/table_prnp_junction_metrics.tex",
]
SUPPLEMENT_TABLE_SOURCES = [
    "manuscript/tables/supplement/table_dna_quality_evidence.tex",
    "manuscript/tables/supplement/table_ddpcr_participant_haploid_genomes.tex",
    "manuscript/tables/supplement/table_prnp_orr_repeat_tool_summary.tex",
    "manuscript/tables/supplement/table_ddpcr_e200k_positive_well_results.tex",
]
EXTERNAL_SUPPLEMENT_TABLE_ASSETS = {
    "Supplementary Table S3": "manuscript/tables/supplement/ddPCR_results_by_region_S3.xlsx",
}
SUPPLEMENT_TABLE_BUNDLE_LABELS = ["Table S1", "Table S2", "Table S4", "Table S5"]

EXPECTED_WRAPPERS = {
    "manuscript/figures/all_figures_main.tex": MAIN_FIGURE_SOURCES,
    "manuscript/figures/all_figures_supplement.tex": SUPPLEMENT_FIGURE_SOURCES,
    "manuscript/tables/all_tables_main.tex": MAIN_TABLE_SOURCES,
    "manuscript/tables/all_tables_supplement.tex": SUPPLEMENT_TABLE_SOURCES,
    "manuscript/figures/all_figures_with_legends/figures_with_legends.tex": [
        "manuscript/figures/all_figures_main.tex",
        "manuscript/figures/all_figures_supplement.tex",
    ],
    "manuscript/tables/all_tables_with_legends/all_tables.tex": [
        "manuscript/tables/all_tables_main.tex",
        "manuscript/tables/all_tables_supplement.tex",
    ],
}

DOCX_W_NS = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
DOCX_DEL = DOCX_W_NS + "del"
DOCX_TEXT = DOCX_W_NS + "t"
PROVENANCE_STATUSES = {
    "manual",
    "reproducible",
    "reproducible_with_manual_review",
    "needs_provenance_recovery",
}
PROVENANCE_MARKERS = {"none", "not_recorded", "manual", "manual/legacy"}


def visible_docx_text(path: Path) -> tuple[str, str]:
    with zipfile.ZipFile(path) as archive:
        xml = ET.fromstring(archive.read("word/document.xml"))

    visible: list[str] = []
    deleted: list[str] = []

    def visit(node: ET.Element, in_deleted: bool = False) -> None:
        now_deleted = in_deleted or node.tag == DOCX_DEL
        if node.tag == DOCX_TEXT and node.text:
            (deleted if now_deleted else visible).append(node.text)
        for child in node:
            visit(child, now_deleted)

    visit(xml)
    return " ".join(visible), " ".join(deleted)


def expand_reference_token(token: str) -> list[str]:
    token = token.replace("–", "-").replace("—", "-").strip()
    if "-" in token and re.fullmatch(r"S?\d+-S?\d+", token, flags=re.IGNORECASE):
        left, right = token.split("-", 1)
        prefix = "S" if left.upper().startswith("S") else ""
        start = int(re.search(r"\d+", left).group())
        end = int(re.search(r"\d+", right).group())
        return [f"{prefix}{i}" for i in range(start, end + 1)]
    return [token.upper() if token.upper().startswith("S") else str(int(token))]


def docx_sequence(text: str, kind: str) -> list[str]:
    pattern = re.compile(
        rf"\b(Supplementary\s+)?{kind}s?\s+([S]?\d+(?:\s*(?:and|,|&)\s*[S]?\d+|\s*[-–—]\s*[S]?\d+)*)",
        flags=re.IGNORECASE,
    )
    seen: list[str] = []
    for match in pattern.finditer(text):
        supplementary = match.group(1) is not None or match.group(2).lstrip().upper().startswith("S")
        body = match.group(2)
        chunks = re.split(r"\s*(?:and|,|&)\s*", body, flags=re.IGNORECASE)
        expanded: list[str] = []
        for chunk in chunks:
            expanded.extend(expand_reference_token(chunk))
        prefix = "Supplementary " if supplementary else ""
        for item in expanded:
            label = f"{prefix}{kind} {item}"
            if label not in seen:
                seen.append(label)
    return seen


def input_paths(path: Path) -> list[str]:
    return re.findall(r"\\input\{([^}]+)\}", path.read_text(encoding="utf-8"))


def normalise_repo_path(root: Path, raw: str) -> Path:
    raw = raw.strip().replace("\\", "/")
    root_text = root.as_posix().rstrip("/").lower()
    if raw.lower().startswith(root_text + "/"):
        return root / raw[len(root_text) + 1 :]
    candidate = Path(raw)
    return candidate if candidate.is_absolute() else root / raw


def split_paths(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def is_provenance_marker(value: str) -> bool:
    value = value.strip().lower()
    return value in PROVENANCE_MARKERS or value.startswith("manual (")


def check_provenance(root: Path, crosswalk_rows: list[dict[str, str]], errors: list[str]) -> None:
    """Validate creation metadata and its one-to-many input dependency table."""

    expected_refs = FIGURES + TABLES
    required_crosswalk_columns = {
        "provenance_status",
        "generator_script",
        "transform_script",
        "asset_output",
        "generated_tex",
        "verification",
    }
    if crosswalk_rows and not required_crosswalk_columns.issubset(crosswalk_rows[0]):
        missing = sorted(required_crosswalk_columns - set(crosswalk_rows[0]))
        errors.append(f"Crosswalk is missing provenance columns: {', '.join(missing)}")
        return

    by_ref = {row.get("docx_reference", ""): row for row in crosswalk_rows}
    for reference in expected_refs:
        row = by_ref.get(reference)
        if row is None:
            continue
        provenance_status = row["provenance_status"].strip()
        if provenance_status not in PROVENANCE_STATUSES:
            errors.append(f"Unknown provenance status for {reference}: {provenance_status!r}")
        if reference in {"Figure 1", "Figure 2", "Figure 3", "Table 1"} and provenance_status != "manual":
            errors.append(f"{reference} must have provenance status 'manual'")
        if row["generated_tex"].strip() not in {"manual", "generated"}:
            errors.append(f"Invalid generated_tex value for {reference}: {row['generated_tex']!r}")
        for column in ("generator_script", "transform_script", "asset_output"):
            for raw_path in split_paths(row[column]):
                if is_provenance_marker(raw_path):
                    continue
                path = normalise_repo_path(root, raw_path)
                if not path.is_file() or path.stat().st_size == 0:
                    errors.append(
                        f"Missing or empty provenance {column} for {reference}: {raw_path}"
                    )

    inputs_path = root / "manuscript/config/figure_table_provenance_inputs.tsv"
    if not inputs_path.is_file() or inputs_path.stat().st_size == 0:
        errors.append(f"Missing or empty provenance input table: {inputs_path.relative_to(root)}")
        return
    with inputs_path.open(newline="", encoding="utf-8") as handle:
        input_rows = list(csv.DictReader(handle, delimiter="\t"))
    required_input_columns = {"docx_reference", "input_role", "path", "required", "produced_by", "notes"}
    if not input_rows or not required_input_columns.issubset(input_rows[0]):
        missing = sorted(required_input_columns - set(input_rows[0] if input_rows else {}))
        errors.append(f"Provenance input table is missing columns: {', '.join(missing)}")
        return
    input_refs = {row["docx_reference"] for row in input_rows}
    if input_refs != set(expected_refs):
        missing = sorted(set(expected_refs) - input_refs)
        extra = sorted(input_refs - set(expected_refs))
        if missing:
            errors.append(f"Provenance input table is missing references: {', '.join(missing)}")
        if extra:
            errors.append(f"Provenance input table has unexpected references: {', '.join(extra)}")
    for index, row in enumerate(input_rows, start=2):
        reference = row["docx_reference"]
        if reference not in by_ref:
            errors.append(f"Provenance input row {index} has unknown reference: {reference}")
        if row["required"].strip().lower() not in {"yes", "no"}:
            errors.append(f"Provenance input row {index} has invalid required flag")
        raw_path = row["path"].strip()
        if not is_provenance_marker(raw_path):
            path = normalise_repo_path(root, raw_path)
            if not path.is_file() or path.stat().st_size == 0:
                errors.append(f"Missing or empty provenance input row {index}: {raw_path}")
    for reference in expected_refs:
        rows = [row for row in input_rows if row["docx_reference"] == reference]
        if not rows:
            continue
        status = by_ref.get(reference, {}).get("provenance_status", "")
        if status in {"reproducible", "reproducible_with_manual_review"}:
            unresolved = [row["path"] for row in rows if is_provenance_marker(row["path"])]
            if unresolved:
                errors.append(
                    f"Reproducible artefact {reference} has unresolved provenance inputs: "
                    + ", ".join(unresolved)
                )


def check_external_supplement_assets(root: Path, errors: list[str]) -> None:
    """Check supplementary tables delivered outside the combined TeX bundle."""
    for reference, raw_path in EXTERNAL_SUPPLEMENT_TABLE_ASSETS.items():
        path = normalise_repo_path(root, raw_path)
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"Missing or empty external asset for {reference}: {raw_path}")


def check_caption_catalogue(root: Path, errors: list[str]) -> None:
    """Check that the plain-text catalogue matches the canonical captions."""
    path = root / CAPTION_CATALOGUE_PATH
    try:
        expected = build_caption_catalogue(root).encode("utf-8")
    except (KeyError, OSError, ValueError) as exc:
        errors.append(f"Could not generate caption catalogue: {exc}")
        return
    if not path.is_file() or path.read_bytes() != expected:
        errors.append(
            "Caption catalogue is stale; run "
            "src/tools/manuscript/generate_caption_list.py"
        )


def referenced_assets(root: Path, entry_points: list[str]) -> tuple[set[Path], list[str]]:
    tex_seen: set[Path] = set()
    assets: set[Path] = set()
    errors: list[str] = []
    pending = [root / path for path in entry_points]
    while pending:
        tex = pending.pop()
        if tex in tex_seen:
            continue
        tex_seen.add(tex)
        if not tex.is_file() or tex.stat().st_size == 0:
            errors.append(f"Missing or empty TeX source: {tex.relative_to(root)}")
            continue
        for raw in input_paths(tex):
            child = normalise_repo_path(root, raw)
            pending.append(child)
        for raw in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", tex.read_text(encoding="utf-8")):
            asset = normalise_repo_path(root, raw)
            if not asset.suffix:
                matches = [asset.with_suffix(ext) for ext in (".pdf", ".svg", ".png", ".jpg", ".jpeg") if asset.with_suffix(ext).is_file()]
                if len(matches) == 1:
                    asset = matches[0]
                elif not matches:
                    errors.append(f"Missing graphics asset: {raw} (referenced by {tex.relative_to(root)})")
                    continue
            assets.add(asset)
            if not asset.is_file() or asset.stat().st_size == 0:
                errors.append(f"Missing or empty graphics asset: {asset}")
    for path in assets:
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"Missing or empty referenced asset: {path}")
    return assets, errors


def extract_pdf_labels(root: Path, path: Path) -> list[str]:
    candidates = [
        r"C:/Users/Manfredi/Documents/Codex/tools/miktex-portable/texmfs/install/miktex/bin/x64/pdftotext.exe",
        str(Path(r"C:/Users/Manfredi/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftotext.exe")),
        shutil.which("pdftotext"),
        r"C:/Users/Manfredi/AppData/Local/Programs/MiKTeX/miktex/bin/x64/pdftotext.exe",
    ]
    executable = next((candidate for candidate in candidates if candidate and Path(candidate).exists()), None)
    if executable is not None:
        result = subprocess.run(
            [executable, "-layout", str(path), "-"],
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="ignore",
            timeout=120,
        )
        text = result.stdout or ""
    else:
        # pypdf is a portable fallback, but it is much slower on the large
        # embedded-asset figure bundle than Poppler's text extractor.
        try:
            from pypdf import PdfReader

            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                text = "\n".join(page.extract_text() or "" for page in PdfReader(str(path), strict=False).pages)
        except Exception as exc:
            raise RuntimeError("pdftotext or pypdf is required to verify bundle caption order") from exc
    labels: list[str] = []
    # Caption labels are followed by a colon; references to another artefact
    # inside a caption (for example, “Table S1.”) are not bundle labels.
    for match in re.finditer(r"\b(Figure|Table)\s+(S?\d+)\s*:", text):
        label = f"{match.group(1)} {match.group(2)}"
        if label not in labels:
            labels.append(label)
    return labels


def check(root: Path) -> list[str]:
    errors: list[str] = []
    docx = root / "manuscript/v2/Carta_PRNP_v2.docx"
    if not docx.is_file() or docx.stat().st_size == 0:
        errors.append(f"Missing authoritative DOCX: {docx}")
    else:
        visible, _deleted = visible_docx_text(docx)
        figure_refs = docx_sequence(visible, "Figure")
        table_refs = docx_sequence(visible, "Table")
        main_figures = [label for label in figure_refs if not label.startswith("Supplementary ")]
        supplementary_figures = [label for label in figure_refs if label.startswith("Supplementary ")]
        main_tables = [label for label in table_refs if not label.startswith("Supplementary ")]
        supplementary_tables = [label for label in table_refs if label.startswith("Supplementary ")]
        if main_figures != FIGURES[:8] or sorted(supplementary_figures) != sorted(FIGURES[8:]):
            errors.append(f"DOCX figure sequences mismatch: main={main_figures}, supplementary={supplementary_figures}")
        if main_tables != TABLES[:7] or sorted(supplementary_tables) != sorted(TABLES[7:]):
            errors.append(f"DOCX table sequences mismatch: main={main_tables}, supplementary={supplementary_tables}")

    crosswalk = root / "manuscript/config/figure_table_crosswalk.tsv"
    with crosswalk.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_rows = FIGURES + TABLES
    if [row["docx_reference"] for row in rows] != expected_rows:
        errors.append("Crosswalk references are not exact, ordered and gap-free")
    if len({row["tex_file"] for row in rows}) != len(rows):
        errors.append("Crosswalk contains duplicate canonical source paths")
    for row in rows:
        tex = root / row["tex_file"]
        if not tex.is_file() or tex.stat().st_size == 0:
            errors.append(f"Crosswalk source is missing or empty: {row['tex_file']}")
        if "Table 3" == row["docx_reference"] and "ddpcr_sample_number.R" not in row["generator_script"]:
            errors.append("Table 3 crosswalk row does not name its generator")

    check_provenance(root, rows, errors)
    check_external_supplement_assets(root, errors)
    check_caption_catalogue(root, errors)

    for wrapper, expected in EXPECTED_WRAPPERS.items():
        actual = [normalise_repo_path(root, raw).as_posix().replace(root.as_posix() + "/", "") for raw in input_paths(root / wrapper)]
        if actual != expected:
            errors.append(f"Wrapper order mismatch for {wrapper}: {actual} != {expected}")

    _, asset_errors = referenced_assets(root, list(EXPECTED_WRAPPERS))
    errors.extend(asset_errors)

    lod = root / "manuscript/tables/lod_calculations/LoD_data.csv"
    if lod.is_file():
        with lod.open(newline="", encoding="utf-8-sig") as handle:
            for index, row in enumerate(csv.DictReader(handle), start=2):
                try:
                    low, high = float(row["ci_low"]), float(row["ci_high"])
                except (KeyError, TypeError, ValueError) as exc:
                    errors.append(f"LoD row {index} has invalid confidence bounds: {exc}")
                    continue
                if low > high:
                    errors.append(f"LoD row {index} has descending confidence bounds: {low} > {high}")
    else:
        errors.append(f"Missing LoD data: {lod}")

    stale_patterns = [
        re.compile(r"Figure 1--7"),
        re.compile(r"Table 1--3"),
        re.compile(r"0\.8\s*%"),
        re.compile(r"Detection of somatic intronic variants in prion disease brain tissue", re.IGNORECASE),
    ]
    scan_paths = [root / path for path in EXPECTED_WRAPPERS]
    scan_paths += [root / "manuscript/README.md", root / "manuscript/config/README.md", root / "manuscript/scripts/README.md", crosswalk]
    for path in scan_paths:
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        for pattern in stale_patterns:
            if pattern.search(text):
                errors.append(f"Stale string {pattern.pattern!r} found in {path.relative_to(root)}")

    bundle_checks = [
        (root / "manuscript/figures/all_figures_with_legends/figures_with_legends.pdf", [f"Figure {i}" for i in range(1, 9)] + [f"Figure S{i}" for i in range(1, 5)]),
        (root / "manuscript/tables/all_tables_with_legends/all_tables.pdf", [f"Table {i}" for i in range(1, 8)] + SUPPLEMENT_TABLE_BUNDLE_LABELS),
    ]
    for pdf, expected in bundle_checks:
        if not pdf.is_file() or pdf.stat().st_size == 0:
            errors.append(f"Missing or empty combined PDF: {pdf.relative_to(root)}")
            continue
        try:
            actual = extract_pdf_labels(root, pdf)
            if actual != expected:
                errors.append(f"PDF caption order mismatch for {pdf.relative_to(root)}: {actual} != {expected}")
        except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
            errors.append(f"Could not verify PDF caption order for {pdf.relative_to(root)}: {exc}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[3])
    args = parser.parse_args()
    errors = check(args.root.resolve())
    if errors:
        print("MANUSCRIPT CONSISTENCY CHECK FAILED")
        for error in errors:
            print(f"- {error}")
        return 1
    print("MANUSCRIPT CONSISTENCY CHECK PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
