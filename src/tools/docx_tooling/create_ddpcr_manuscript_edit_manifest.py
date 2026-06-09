#!/usr/bin/env python3
"""Write the ddPCR manuscript replacement manifest used by Word automation."""

from __future__ import annotations

import json
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import create_ddpcr_corrected_manuscript_docx as manuscript


PROJECT_ROOT = Path(__file__).resolve().parents[3]
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript" / "v2"
MANIFEST_PATH = MANUSCRIPT_DIR / "ddpcr_manuscript_replacements.json"


def compute_word_counts() -> dict[str, int]:
    """Recreate the edited manuscript text and count abstract/main words."""
    with zipfile.ZipFile(manuscript.ORIGINAL_DOCX, "r") as zf:
        root = ET.fromstring(zf.read("word/document.xml"))

    body = root.find(manuscript.qn(manuscript.W_NS, "body"))
    if body is None:
        raise RuntimeError("Could not find Word document body")

    paragraphs = manuscript.all_paragraphs(root)
    anchors = {}
    for prefix, new_text in manuscript.REPLACEMENTS.items():
        anchors[prefix] = manuscript.replace_paragraph_by_prefix(paragraphs, prefix, new_text)
        paragraphs = manuscript.all_paragraphs(root)

    data_processing = anchors["Data processing: Raw droplet fluorescence data were exported"]
    lambda_p = manuscript.insert_after(
        body,
        data_processing,
        manuscript.formula_paragraph(manuscript.FORMULAS["lambda"]),
    )
    manuscript.insert_after(body, lambda_p, manuscript.formula_paragraph(manuscript.FORMULAS["fa"]))

    lob = anchors["Limit of blank: To rule out false positives"]
    manuscript.insert_after(body, lob, manuscript.formula_paragraph(manuscript.FORMULAS["lob"]))

    gtotal_old = [p for p in manuscript.all_paragraphs(root) if manuscript.paragraph_text(p).startswith("Gtotal =")]
    if len(gtotal_old) != 1:
        raise RuntimeError(f"Expected one Gtotal paragraph, found {len(gtotal_old)}")
    manuscript.replace_paragraph_with_formula(gtotal_old[0], manuscript.FORMULAS["gtotal"])

    texts = [manuscript.paragraph_text(p) for p in manuscript.all_paragraphs(root) if manuscript.paragraph_text(p)]
    abstract_start = texts.index("Abstract")
    intro_start = texts.index("Introduction")
    acknowledgements_start = texts.index("Acknowledgements")

    return {
        "abstract": manuscript.count_words(" ".join(texts[abstract_start + 1:intro_start])),
        "main": manuscript.count_words(" ".join(texts[intro_start + 1:acknowledgements_start])),
    }


def main() -> None:
    replacements = [
        {"prefix": prefix, "text": text}
        for prefix, text in manuscript.REPLACEMENTS.items()
    ]
    formulas = manuscript.FORMULAS

    manifest = {
        "replacements": replacements,
        "formula_insertions": [
            {
                "after_prefix": manuscript.REPLACEMENTS[
                    "Data processing: Raw droplet fluorescence data were exported"
                ][:90],
                "text": formulas["lambda"],
            },
            {"after_prefix": formulas["lambda"], "text": formulas["fa"]},
            {
                "after_prefix": manuscript.REPLACEMENTS[
                    "Limit of blank: To rule out false positives"
                ][:90],
                "text": formulas["lob"],
            },
        ],
        "formula_replacements": [
            {"prefix": "Gtotal =", "text": formulas["gtotal"]},
        ],
        "word_counts": compute_word_counts(),
    }

    MANIFEST_PATH.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(MANIFEST_PATH.as_posix())
    print(f"Abstract word count: {manifest['word_counts']['abstract']}")
    print(f"Main text word count: {manifest['word_counts']['main']}")


if __name__ == "__main__":
    main()
