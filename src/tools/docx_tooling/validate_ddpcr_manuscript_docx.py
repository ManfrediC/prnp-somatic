#!/usr/bin/env python3
"""Validate the corrected ddPCR manuscript DOCX artefacts."""

from __future__ import annotations

import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import create_ddpcr_corrected_manuscript_docx as manuscript


PROJECT_ROOT = Path(__file__).resolve().parents[3]
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript" / "v2"
CLEAN_DOCX = MANUSCRIPT_DIR / "2026_05_11_Carta_PRNP_v2.ddpcr_corrected_clean_word.docx"
TRACKED_DOCX = MANUSCRIPT_DIR / "2026_05_11_Carta_PRNP_v2.ddpcr_corrected_tracked.docx"

W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
M_NS = "http://schemas.openxmlformats.org/officeDocument/2006/math"


def qn(ns: str, tag: str) -> str:
    return f"{{{ns}}}{tag}"


def extract_text(xml: str, include_deletions: bool = True) -> str:
    root = ET.fromstring(xml)
    text_parts: list[str] = []
    for node in root.iter():
        text_tags = {qn(W_NS, "t"), qn(M_NS, "t")}
        if include_deletions:
            text_tags.add(qn(W_NS, "delText"))
        if node.tag in text_tags:
            text_parts.append(node.text or "")
    return "".join(text_parts)


def read_document_xml(path: Path) -> str:
    with zipfile.ZipFile(path, "r") as zf:
        bad_member = zf.testzip()
        if bad_member is not None:
            raise RuntimeError(f"Bad compressed member in {path}: {bad_member}")
        return zf.read("word/document.xml").decode("utf-8")


def main() -> int:
    if not TRACKED_DOCX.exists():
        print(f"Missing tracked DOCX: {TRACKED_DOCX}", file=sys.stderr)
        return 1
    if not CLEAN_DOCX.exists():
        print(f"Missing clean revised DOCX: {CLEAN_DOCX}", file=sys.stderr)
        return 1

    try:
        document_xml = read_document_xml(TRACKED_DOCX)
        clean_xml = read_document_xml(CLEAN_DOCX)
    except Exception as exc:
        print(f"Could not read DOCX artefacts: {exc}", file=sys.stderr)
        return 1

    tracked_text = extract_text(document_xml)
    clean_text = extract_text(clean_xml, include_deletions=False)
    summary = manuscript.DDPCR_SUMMARY
    checks = {
        "tracked insertions": document_xml.count("<w:ins") > 0,
        "tracked deletions": document_xml.count("<w:del") > 0,
        "technical pooled E200K calls": str(summary["technical_pooled_pass_description"]) in clean_text,
        "zero non-pooled call": "no non-pooled non-CJD30 sample-region measurement passed both" in clean_text,
        "participant pooled positive": str(summary["participant_pass_description"]) in clean_text,
        "participant pooled follow-up": str(summary["pooled_follow_up_description"]) in clean_text,
        "LoB-only count": f"further {summary['lob_only_count']} sample-region measurements" in clean_text,
        "count-scale LoB": "count-scale LoB" in clean_text,
        "Poisson-corrected FA": "Poisson-corrected molecule-occupancy" in clean_text,
        "LoB formula": "LoB_count" in tracked_text,
    }

    print(f"Clean revised DOCX: {CLEAN_DOCX.as_posix()}")
    print(f"Tracked DOCX: {TRACKED_DOCX.as_posix()}")
    print(f"w:ins count: {document_xml.count('<w:ins')}")
    print(f"w:del count: {document_xml.count('<w:del')}")
    for label, ok in checks.items():
        print(f"{label}: {'PASS' if ok else 'FAIL'}")

    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
