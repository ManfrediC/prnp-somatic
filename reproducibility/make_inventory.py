#!/usr/bin/env python3

"""
Create a script inventory for active (non-legacy) workflow and manuscript code.

Output:
- reproducibility/inventory.tsv

Current scope:
- include script-like files in src/;
- include canonical manuscript figure/legend scripts under manuscript/figures/;
- exclude legacy and historical manuscript figure directories.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from datetime import datetime, timezone

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOTS = (REPO_ROOT / "src", REPO_ROOT / "manuscript" / "figures")
EXCLUDED_DIR_NAMES = {"legacy", "old_versions", "height_4_5"}
OUT_PATH = REPO_ROOT / "reproducibility" / "inventory.tsv"

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
    for script_root in SCRIPT_ROOTS:
        if not script_root.exists():
            continue
        for p in sorted(script_root.rglob("*")):
            if not p.is_file():
                continue
            # Inventory active scripts only: skip archived or historical trees.
            if any(part in EXCLUDED_DIR_NAMES for part in p.relative_to(script_root).parts):
                continue
            lang = EXT_TO_LANG.get(p.suffix, "")
            if not lang:
                continue  # Ignore non-script files for now.

            stat = p.stat()
            mtime_utc = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat()

            # Keep a stable row schema; use NA placeholders to avoid trailing TSV whitespace.
            rows.append({
                "relpath": p.relative_to(REPO_ROOT).as_posix(),
                "language": lang,
                "size_bytes": str(stat.st_size),
                "mtime_utc": mtime_utc,
                "sha256": sha256_file(p),
                "purpose": "NA",
                "inputs": "NA",
                "outputs": "NA",
                "deps/tools": "NA",
                "notes": "NA",
            })

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w", newline="", encoding="utf-8") as f:
        if not rows:
            raise SystemExit(f"No active scripts discovered under {SCRIPT_ROOTS}")
        w = csv.DictWriter(
            f,
            fieldnames=list(rows[0].keys()),
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUT_PATH}")

if __name__ == "__main__":
    missing_roots = [str(root) for root in SCRIPT_ROOTS if not root.exists()]
    if len(missing_roots) == len(SCRIPT_ROOTS):
        raise SystemExit(f"Missing script roots: {', '.join(missing_roots)}")
    main()
