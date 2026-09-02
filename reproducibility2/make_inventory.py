#!/usr/bin/env python3
"""Create the isolated sequencing2 script inventory."""

from __future__ import annotations

import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOTS = (
    REPO_ROOT / "src2",
    REPO_ROOT / "manuscript2",
    REPO_ROOT / "tests2",
)
OUTPUT = REPO_ROOT / "reproducibility2" / "inventory.tsv"
EXTENSION_LANGUAGE = {".py": "python", ".sh": "bash", ".R": "r", ".r": "r"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    rows: list[dict[str, str]] = []
    for root in SCRIPT_ROOTS:
        for path in sorted(root.rglob("*")):
            language = EXTENSION_LANGUAGE.get(path.suffix)
            if not path.is_file() or language is None:
                continue
            stat = path.stat()
            rows.append(
                {
                    "relpath": path.relative_to(REPO_ROOT).as_posix(),
                    "language": language,
                    "size_bytes": str(stat.st_size),
                    "mtime_utc": datetime.fromtimestamp(
                        stat.st_mtime, tz=timezone.utc
                    ).isoformat(),
                    "sha256": sha256_file(path),
                    "purpose": "sequencing2 isolated workflow or verification",
                    "inputs": "see manuscript2/config2/figure_table_provenance_inputs.tsv",
                    "outputs": "see reproducibility2/final_outputs_manifest.tsv",
                    "deps/tools": "project environment",
                    "notes": "not part of the original pipeline",
                }
            )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {OUTPUT}")


if __name__ == "__main__":
    main()
