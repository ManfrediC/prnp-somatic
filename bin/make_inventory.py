#!/usr/bin/env python3

"""
Create a script inventory for active (non-legacy) code under src/.

Output:
- doc/inventory.tsv

Current scope:
- include script-like files in src/
- exclude src/legacy/
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from datetime import datetime, timezone

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
LEGACY_DIR = SRC_DIR / "legacy"
OUT_PATH = REPO_ROOT / "doc" / "inventory.tsv"

EXT_TO_LANG = {
    ".py": "python",
    ".sh": "bash",
    ".R": "r",
    ".r": "r",
    ".pl": "perl",
    ".awk": "awk",
}

def sha256_file(path: Path) -> str:
    # File hash is used to track content-level changes across snapshots.
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def main() -> None:
    rows = []
    for p in sorted(SRC_DIR.rglob("*")):
        if not p.is_file():
            continue
        # Inventory active scripts only: skip archived legacy tree.
        if LEGACY_DIR in p.parents:
            continue
        lang = EXT_TO_LANG.get(p.suffix, "")
        if not lang:
            continue  # Ignore non-script files for now.

        stat = p.stat()
        mtime_utc = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat()

        # Keep a stable row schema so downstream curation can fill blanks in-place.
        rows.append({
            "relpath": str(p.relative_to(REPO_ROOT)),
            "language": lang,
            "size_bytes": str(stat.st_size),
            "mtime_utc": mtime_utc,
            "sha256": sha256_file(p),
            "purpose": "",
            "inputs": "",
            "outputs": "",
            "deps/tools": "",
            "notes": "",
        })

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w", newline="", encoding="utf-8") as f:
        if not rows:
            raise SystemExit(f"No scripts discovered under {SRC_DIR} (excluding {LEGACY_DIR})")
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUT_PATH}")

if __name__ == "__main__":
    if not SRC_DIR.exists():
        raise SystemExit(f"Missing src dir: {SRC_DIR}")
    main()
