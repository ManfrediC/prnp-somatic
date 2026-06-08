#!/usr/bin/env python3
"""Build the local ddPCR LoD raw-file database and provenance manifests."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import posixpath
import re
import shutil
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path


KEY_ENV_PATH = Path(__file__).resolve().parents[2] / "env" / "ddpcr_key.env"
FORBIDDEN_NAME_RE = re.compile(r"(Manfredi|Mirka)", re.IGNORECASE)
WELL_RE = re.compile(r"^[A-H](?:0[1-9]|1[0-2])$")


@dataclass(frozen=True)
class LodRun:
    slug: str
    assay: str
    run_date: str
    experiment_title: str
    ddpcr_source: str
    raw_csv_source: str
    analysis_csv_source: str
    layout_source: str
    excluded: tuple[tuple[str, str, str], ...] = ()

    @property
    def run_id(self) -> str:
        return f"{self.run_date}_SNV_{self.assay}"


LOD_RUNS = (
    LodRun(
        slug="d178n",
        assay="D178N",
        run_date="2020-09-11",
        experiment_title="2020-09-11 D178N LOD",
        ddpcr_source=r"Runs\2020-09-11 D178N LOD\Manfredi_D178N_LOD.ddpcr",
        raw_csv_source=(
            r"Runs\2020-09-11 D178N LOD"
            r"\Manfredi_D178N_LOD_2020-09-11-18-27\Manfredi_D178N_LOD.csv"
        ),
        analysis_csv_source=r"LOD graphs\D178N\Manfredi_D178N_LOD.csv",
        layout_source=r"Runs\2020-09-11 D178N LOD\2020-09-11 SNV LOD D178N.xlsx",
        excluded=(
            (
                "ddpcr_archive",
                r"Runs\2020_08_05 ddPCR first test"
                r"\Manfredi_PRNP_RPPH1_2020-08-05-10-44\Manfredi_D178N_LOD.ddpcr",
                "early_test_not_final_lod_source",
            ),
        ),
    ),
    LodRun(
        slug="e200k",
        assay="E200K",
        run_date="2020-10-08",
        experiment_title="2020-10-08 P102L test and E200K LOD",
        ddpcr_source=(
            r"Runs\2020-10-08 P102L test and E200K LOD"
            r"\Manfredi_SNV_E200K_LOD_only.ddpcr"
        ),
        raw_csv_source=(
            r"Runs\2020-10-08 P102L test and E200K LOD"
            r"\Manfredi_SNV_E200K_P102L_2020-10-08-19-38"
            r"\Manfredi_SNV_E200K_P102L.csv"
        ),
        analysis_csv_source=r"LOD graphs\E200K\Manfredi_E200K_LOD.csv",
        layout_source=(
            r"Runs\2020-10-08 P102L test and E200K LOD"
            r"\2020-10-08 P102L test and E200K LOD.xlsx"
        ),
        excluded=(),
    ),
    LodRun(
        slug="p102l",
        assay="P102L",
        run_date="2021-01-25",
        experiment_title="2021-01-25 SNV LOD P102L, CNV LOD",
        ddpcr_source=(
            r"Runs\2021-01-25 SNV LOD P102L, CNV LOD"
            r"\Manfredi-2021-01-25 SNV LOD P102L only.ddpcr"
        ),
        raw_csv_source=(
            r"Runs\2021-01-25 SNV LOD P102L, CNV LOD"
            r"\Manfredi-2021-01-25 SNV LOD P102L only.csv"
        ),
        analysis_csv_source=r"LOD graphs\P102L\Manfredi_P102L_LOD.csv",
        layout_source=(
            r"Runs\2021-01-25 SNV LOD P102L, CNV LOD"
            r"\2021-01-25 SNV LOD P102L, CNV LOD.xlsx"
        ),
        excluded=(
            (
                "csv_export",
                r"Runs\2021-01-22 SNV LOD P102L, CNV LOD"
                r"\Manfredi-2021-01-22 SNV LOD P102L - CNV LOD_2021-01-22-17-44"
                r"\Manfredi-2021-01-22 SNV LOD P102L - CNV LOD.csv",
                "superseded_not_final_lod_source",
            ),
            (
                "csv_export",
                r"Runs\2020-11-24 P102L and rescue samples\Manfredi-2020-11-24_P102L_LOD_only.csv",
                "superseded_not_final_lod_source",
            ),
            (
                "csv_export",
                r"LOD graphs\P102L\old\Manfredi_P102L_LOD_combined.csv",
                "superseded_not_final_lod_source",
            ),
        ),
    ),
)

FILE_FIELDS = [
    "assay",
    "run_id",
    "file_kind",
    "source_path",
    "destination_path",
    "source_name",
    "destination_name",
    "sha256",
    "size_bytes",
    "modified_time",
]
RUN_FIELDS = [
    "assay",
    "run_id",
    "run_date",
    "experiment_title",
    "ddpcr_archive",
    "archive_contents_dir",
    "csv_export",
    "analysis_csv",
    "layout_xlsx",
    "source_ddpcr_archive",
    "source_csv_export",
    "source_analysis_csv",
    "source_layout_xlsx",
    "notes",
]
SAMPLE_FIELDS = [
    "assay",
    "run_id",
    "well",
    "sample",
    "sample_normalised",
    "target",
    "target_normalised",
    "accepted_droplets",
    "positives",
    "negatives",
    "merged_wells",
    "fractional_abundance",
    "source_csv",
]
ARCHIVE_FIELDS = [
    "assay",
    "run_id",
    "archive_member_path",
    "destination_path",
    "content_kind",
    "well",
    "target_or_channel",
    "sha256",
    "size_bytes",
]
EXCLUDED_FIELDS = [
    "assay",
    "file_kind",
    "source_path",
    "reason",
    "sha256",
    "size_bytes",
    "exists",
]


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


def path_from_text(value: str) -> Path:
    if re.match(r"^[A-Za-z]:\\", value):
        drive = value[0].lower()
        return Path(f"/mnt/{drive}") / value[3:].replace("\\", "/")
    return Path(value)


def source_path(root: Path, relative_windows_path: str) -> Path:
    return root / relative_windows_path.replace("\\", "/")


def windows_path(path: Path) -> str:
    resolved = path.resolve()
    text = resolved.as_posix()
    if text.startswith("/mnt/c/"):
        return "C:\\" + text[len("/mnt/c/") :].replace("/", "\\")
    return str(resolved)


def relative_path(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def modified_time(path: Path) -> str:
    return str(int(path.stat().st_mtime))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def copy_file(source: Path, destination: Path, dry_run: bool) -> None:
    if dry_run:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def safe_member_path(destination: Path, member_name: str) -> Path:
    normalised = posixpath.normpath(member_name.replace("\\", "/"))
    if normalised.startswith("../") or normalised == ".." or posixpath.isabs(normalised):
        raise ValueError(f"Unsafe archive member path: {member_name}")
    target = (destination / Path(*normalised.split("/"))).resolve()
    destination_resolved = destination.resolve()
    if destination_resolved not in target.parents and target != destination_resolved:
        raise ValueError(f"Archive member escapes destination: {member_name}")
    return target


def canonical_member_name(member_name: str, run: LodRun) -> str:
    member_path = Path(member_name.replace("\\", "/"))
    if member_path.suffix.lower() == ".ddplt":
        return f"{run.run_id}.ddplt"
    return member_name


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


def extract_archive(source: Path, destination: Path, run: LodRun, dry_run: bool) -> list[dict[str, str]]:
    password = load_ddpcr_key()
    if not dry_run:
        if destination.exists():
            shutil.rmtree(destination)
        destination.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    with zipfile.ZipFile(source) as archive:
        for member in archive.infolist():
            canonical_member = canonical_member_name(member.filename, run)
            target = safe_member_path(destination, canonical_member)
            if member.is_dir():
                if not dry_run:
                    target.mkdir(parents=True, exist_ok=True)
                continue
            if not dry_run:
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(member, pwd=password) as opened, target.open("wb") as output:
                    shutil.copyfileobj(opened, output)
                digest = sha256_file(target)
                size = str(target.stat().st_size)
            else:
                digest = ""
                size = str(member.file_size)
            kind, well, target_or_channel = content_kind(canonical_member)
            if dry_run:
                if canonical_member != member.filename:
                    print(
                        "DRY-RUN rename archive member: "
                        f"{member.filename} -> {canonical_member}"
                    )
                print(f"DRY-RUN export archive member: {member.filename} -> {windows_path(target)}")
            rows.append(
                {
                    "assay": run.assay,
                    "run_id": run.run_id,
                    "archive_member_path": member.filename,
                    "destination_path": windows_path(target),
                    "content_kind": kind,
                    "well": well,
                    "target_or_channel": target_or_channel,
                    "sha256": digest,
                    "size_bytes": size,
                }
            )
    return rows


def normalise_sample(value: str, assay: str) -> str:
    text = value.strip()
    text = re.sub(rf"^{assay}[-_ ]?", "", text, flags=re.IGNORECASE)
    text = re.sub(r"-?mut", "", text, flags=re.IGNORECASE)
    text = re.sub(r"[-_ ]?high$", "", text, flags=re.IGNORECASE)
    text = text.replace("%", "")
    text = text.replace("_", " ").strip()
    if text.upper() in {"NTC", "WT"}:
        return text.upper()
    return text


def normalise_target(value: str) -> str:
    text = value.strip()
    text = re.sub(r"[-_ ]?mut$", "", text, flags=re.IGNORECASE)
    text = text.replace("D1789N", "D178N")
    return text


def get_first(row: dict[str, str], names: tuple[str, ...]) -> str:
    for name in names:
        if name in row:
            return row[name]
    return ""


def parse_sample_manifest(path: Path, run: LodRun) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            target = normalise_target(get_first(row, ("Target", "target")))
            if target not in {run.assay, "WT"}:
                continue
            well = get_first(row, ("Well", "well"))
            if not well.startswith("M"):
                continue
            sample = get_first(row, ("Sample", "sample"))
            rows.append(
                {
                    "assay": run.assay,
                    "run_id": run.run_id,
                    "well": well,
                    "sample": sample,
                    "sample_normalised": normalise_sample(sample, run.assay),
                    "target": target,
                    "target_normalised": target,
                    "accepted_droplets": get_first(row, ("Accepted Droplets", "Accepted.Droplets", "AcceptedDroplets")),
                    "positives": get_first(row, ("Positives", "positives")),
                    "negatives": get_first(row, ("Negatives", "negatives")),
                    "merged_wells": get_first(row, ("MergedWells", "Merged Wells", "merged_wells")),
                    "fractional_abundance": get_first(row, ("FractionalAbundance", "Fractional Abundance")),
                    "source_csv": windows_path(path),
                }
            )
    return rows


def file_manifest_row(run: LodRun, kind: str, source: Path, destination: Path) -> dict[str, str]:
    return {
        "assay": run.assay,
        "run_id": run.run_id,
        "file_kind": kind,
        "source_path": windows_path(source),
        "destination_path": windows_path(destination),
        "source_name": source.name,
        "destination_name": destination.name,
        "sha256": sha256_file(destination),
        "size_bytes": str(destination.stat().st_size),
        "modified_time": modified_time(source),
    }


def excluded_rows(run: LodRun, root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for kind, relative_source, reason in run.excluded:
        path = source_path(root, relative_source)
        exists = path.exists()
        rows.append(
            {
                "assay": run.assay,
                "file_kind": kind,
                "source_path": windows_path(path),
                "reason": reason,
                "sha256": sha256_file(path) if exists and path.is_file() else "",
                "size_bytes": str(path.stat().st_size) if exists and path.is_file() else "",
                "exists": str(exists).lower(),
            }
        )
    return rows


def expected_destinations(raw_root: Path, run: LodRun) -> dict[str, Path]:
    assay_root = raw_root / run.slug
    return {
        "assay_root": assay_root,
        "ddpcr": assay_root / "ddpcr_archive" / f"{run.run_id}.ddpcr",
        "archive_contents": assay_root / "archive_contents",
        "raw_csv": assay_root / "csv_export" / f"{run.run_id}_quantaSoft.csv",
        "analysis_csv": assay_root / "analysis_csv" / f"{run.run_id}_lod_analysis.csv",
        "layout": assay_root / "layout_xlsx" / f"{run.run_id}_layout.xlsx",
        "manifests": assay_root / "manifests",
    }


def validate_identifier_free(raw_root: Path) -> list[str]:
    failures: list[str] = []
    for path in raw_root.rglob("*"):
        if FORBIDDEN_NAME_RE.search(path.name):
            failures.append(relative_path(path, raw_root))
    return failures


def validate_archive_contents(archive_dir: Path) -> list[str]:
    failures: list[str] = []
    if not (archive_dir / "RunInfo.json").exists():
        failures.append("missing RunInfo.json")
    if not (archive_dir / "PeakData").is_dir():
        failures.append("missing PeakData/")
    if not (archive_dir / "PeakMetaData").is_dir():
        failures.append("missing PeakMetaData/")
    return failures


def build_run(run: LodRun, source_root: Path, raw_root: Path, dry_run: bool) -> dict[str, int]:
    dest = expected_destinations(raw_root, run)
    source_files = {
        "ddpcr_archive": source_path(source_root, run.ddpcr_source),
        "csv_export": source_path(source_root, run.raw_csv_source),
        "analysis_csv": source_path(source_root, run.analysis_csv_source),
        "layout_xlsx": source_path(source_root, run.layout_source),
    }
    missing = [f"{kind}: {path}" for kind, path in source_files.items() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing LoD source files:\n" + "\n".join(missing))

    if not dry_run and dest["assay_root"].exists():
        shutil.rmtree(dest["assay_root"])

    planned = [
        (source_files["ddpcr_archive"], dest["ddpcr"], "ddpcr_archive"),
        (source_files["csv_export"], dest["raw_csv"], "csv_export"),
        (source_files["analysis_csv"], dest["analysis_csv"], "analysis_csv"),
        (source_files["layout_xlsx"], dest["layout"], "layout_xlsx"),
    ]
    for source, destination, kind in planned:
        print(f"{'DRY-RUN copy' if dry_run else 'copy'} {kind}: {windows_path(source)} -> {windows_path(destination)}")
        copy_file(source, destination, dry_run)

    print(
        f"{'DRY-RUN extract' if dry_run else 'extract'} archive: "
        f"{windows_path(source_files['ddpcr_archive'])} -> {windows_path(dest['archive_contents'])}"
    )
    archive_rows = extract_archive(source_files["ddpcr_archive"], dest["archive_contents"], run, dry_run)

    if dry_run:
        for row in excluded_rows(run, source_root):
            print(
                "DRY-RUN exclude candidate: "
                f"{row['file_kind']} {row['source_path']} ({row['reason']}; exists={row['exists']})"
            )
        return {"files": len(planned), "archive_contents": len(archive_rows), "sample_rows": 0}

    file_rows = [file_manifest_row(run, kind, source, destination) for source, destination, kind in planned]
    sample_rows = parse_sample_manifest(dest["analysis_csv"], run)
    run_rows = [
        {
            "assay": run.assay,
            "run_id": run.run_id,
            "run_date": run.run_date,
            "experiment_title": run.experiment_title,
            "ddpcr_archive": windows_path(dest["ddpcr"]),
            "archive_contents_dir": windows_path(dest["archive_contents"]),
            "csv_export": windows_path(dest["raw_csv"]),
            "analysis_csv": windows_path(dest["analysis_csv"]),
            "layout_xlsx": windows_path(dest["layout"]),
            "source_ddpcr_archive": windows_path(source_files["ddpcr_archive"]),
            "source_csv_export": windows_path(source_files["csv_export"]),
            "source_analysis_csv": windows_path(source_files["analysis_csv"]),
            "source_layout_xlsx": windows_path(source_files["layout_xlsx"]),
            "notes": "Final ddPCR LoD source curated from USZ experiment folder.",
        }
    ]

    manifest_root = dest["manifests"]
    write_csv(manifest_root / "run.csv", run_rows, RUN_FIELDS)
    write_csv(manifest_root / "files.csv", file_rows, FILE_FIELDS)
    write_csv(manifest_root / "sample_manifest.csv", sample_rows, SAMPLE_FIELDS)
    write_csv(manifest_root / "archive_contents.csv", archive_rows, ARCHIVE_FIELDS)
    write_csv(manifest_root / "excluded_material.csv", excluded_rows(run, source_root), EXCLUDED_FIELDS)

    archive_failures = validate_archive_contents(dest["archive_contents"])
    if archive_failures:
        raise RuntimeError(f"{run.run_id} archive validation failed: {', '.join(archive_failures)}")

    return {"files": len(file_rows), "archive_contents": len(archive_rows), "sample_rows": len(sample_rows)}


def write_manifest_index(raw_root: Path) -> None:
    rows = []
    for run in LOD_RUNS:
        dest = expected_destinations(raw_root, run)
        rows.append(
            {
                "assay": run.assay,
                "run_id": run.run_id,
                "manifest_dir": windows_path(dest["manifests"]),
            }
        )
    write_csv(raw_root / "manifest_index.csv", rows, ["assay", "run_id", "manifest_dir"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument(
        "--usz-ddpcr-root",
        type=Path,
        default=Path(
            "/mnt/c/Users/Manfredi/USZ/Neuropathologie - Carta Manfredi/"
            "CJD PRNP/Experiments/ddPCR"
        ),
    )
    parser.add_argument("--raw-root", type=Path, default=None)
    parser.add_argument("--dry-run", action="store_true", help="List planned work without writing files.")
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    raw_root = (args.raw_root or repo_root / "raw/ddpcr_lod").resolve()
    expected_raw_root = (repo_root / "raw/ddpcr_lod").resolve()
    if raw_root != expected_raw_root:
        print(f"Refusing to write outside raw/ddpcr_lod: {raw_root}", file=sys.stderr)
        return 2

    source_root = path_from_text(str(args.usz_ddpcr_root)).resolve()
    if not source_root.exists():
        print(f"USZ ddPCR root does not exist: {source_root}", file=sys.stderr)
        return 2

    totals = {"files": 0, "archive_contents": 0, "sample_rows": 0}
    for run in LOD_RUNS:
        counts = build_run(run, source_root, raw_root, args.dry_run)
        for key, value in counts.items():
            totals[key] += value

    if not args.dry_run:
        write_manifest_index(raw_root)
        forbidden = validate_identifier_free(raw_root)
        if forbidden:
            print("Identifier-bearing destination names found:", file=sys.stderr)
            for path in forbidden:
                print(f"  {path}", file=sys.stderr)
            return 1

    print("ddPCR LoD raw database build summary")
    print(f"raw_root: {windows_path(raw_root)}")
    print(f"assay_runs: {len(LOD_RUNS)}")
    print(f"copied_files: {totals['files']}")
    print(f"archive_content_rows: {totals['archive_contents']}")
    print(f"sample_manifest_rows: {totals['sample_rows']}")
    print(f"dry_run: {str(args.dry_run).lower()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
