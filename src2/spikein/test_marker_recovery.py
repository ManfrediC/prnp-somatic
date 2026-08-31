#!/usr/bin/env python3
"""Focused first-stage preflight checks; never run a variant caller.

Run with Python from the WSL prnp-spikein environment. Temporary fixtures stay
under results2/spikein/work and are removed after each check. Later stages will
add their own relevant checks here when they are implemented.
"""

import csv
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "src2/spikein/1_sources_haplotypecaller.sh"
MANIFEST = REPO_ROOT / "src2/spikein/samples.tsv"


class SourcePreflightTests(unittest.TestCase):
    def setUp(self):
        work = REPO_ROOT / "results2/spikein/work"
        work.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(prefix="stage1-check-", dir=work)
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.output = self.directory / "discovery"
        with MANIFEST.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.fields = reader.fieldnames
            self.rows = list(reader)

        # A caller invocation is a test failure, even if DRY_RUN regresses.
        guard = self.directory / "bin"
        guard.mkdir()
        gatk = guard / "gatk"
        gatk.write_text("#!/bin/sh\necho CALLER_MUST_NOT_RUN >&2\nexit 97\n")
        gatk.chmod(0o755)
        self.env = os.environ.copy()
        self.env.update(
            PATH=str(guard) + os.pathsep + self.env["PATH"],
            OUT_ROOT=str(self.output),
            SAMPLES_TSV=str(MANIFEST),
            JAVA_MEM_GB="8",
        )
        self.env.pop("DRY_RUN", None)

    def write_manifest(self):
        path = self.directory / "samples.tsv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, self.fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(self.rows)
        self.env["SAMPLES_TSV"] = str(path)

    def invoke(self):
        return subprocess.run(
            ["bash", str(SCRIPT)], cwd=REPO_ROOT, env=self.env,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=60,
        )

    def assert_rejected(self, message):
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn(message, result.stdout)
        self.assertNotIn("CALLER_MUST_NOT_RUN", result.stdout)
        self.assertFalse(self.output.exists())

    def test_default_dry_run_only_calls_pure_samples_and_writes_nothing(self):
        # Selection must follow role names, not manifest row order.
        by_role = {row["role"]: row for row in self.rows}
        self.rows.reverse()
        self.write_manifest()
        before = sorted(str(p.relative_to(self.directory)) for p in self.directory.rglob("*"))
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("DRY RUN COMPLETE", result.stdout)
        self.assertNotIn("CALLER_MUST_NOT_RUN", result.stdout)
        after = sorted(str(p.relative_to(self.directory)) for p in self.directory.rglob("*"))
        self.assertEqual(before, after)
        commands = [shlex.split(line[2:]) for line in result.stdout.splitlines() if line.startswith("+ gatk ")]
        self.assertEqual(len(commands), 4)
        calls = [c for c in commands if "HaplotypeCaller" in c]
        self.assertEqual(len(calls), 2)
        self.assertEqual(
            [c[c.index("-I") + 1] for c in calls],
            [str((REPO_ROOT / by_role[role]["bam"]).resolve()) for role in ("pure_donor", "pure_wt")],
        )
        for call in calls:
            self.assertEqual(call[call.index("-ERC") + 1], "GVCF")
            self.assertEqual(call[call.index("--sample-ploidy") + 1], "2")
            self.assertEqual(call[call.index("--max-reads-per-alignment-start") + 1], "0")
            self.assertEqual(call[call.index("--read-index") + 1], call[call.index("-I") + 1] + ".bai")
        self.assertIn("CombineGVCFs", commands[2])
        self.assertIn("GenotypeGVCFs", commands[3])
        self.assertNotIn("+ samtools index", result.stdout)
        self.assertFalse(self.output.exists())

    def test_missing_role(self):
        self.rows.pop()
        self.write_manifest()
        self.assert_rejected("Manifest must contain exactly")

    def test_duplicate_role(self):
        self.rows[1]["role"] = self.rows[0]["role"]
        self.write_manifest()
        self.assert_rejected("Duplicate sample role")

    def test_unexpected_role(self):
        self.rows[0]["role"] = "other_donor"
        self.write_manifest()
        self.assert_rejected("Unexpected sample role")

    def test_reused_bam(self):
        self.rows[1]["bam"] = self.rows[0]["bam"]
        self.write_manifest()
        self.assert_rejected("BAM file reused across roles")

    def test_mismatched_sample_name(self):
        self.rows[0]["sample_id"] = "wrong_sample"
        self.write_manifest()
        self.assert_rejected("BAM SM does not match manifest")

    def test_missing_bam(self):
        self.rows[0]["bam"] = str(self.directory / "missing.bam")
        self.write_manifest()
        self.assert_rejected("Missing, empty or unreadable BAM")

    def test_missing_bai_is_not_created(self):
        # An empty synthetic BAM tests the missing-index path without changing
        # or copying a project BAM. Only its header is needed for preflight.
        header = subprocess.check_output(
            ["samtools", "view", "-H", str(REPO_ROOT / self.rows[0]["bam"])], text=True,
        )
        bam = self.directory / "no_index.bam"
        subprocess.run(
            ["samtools", "view", "--no-PG", "-b", "-o", str(bam), "-"],
            input=header, text=True, check=True,
        )
        self.rows[0]["bam"] = str(bam)
        self.write_manifest()
        self.assert_rejected("Missing existing BAI")
        self.assertFalse(Path(str(bam) + ".bai").exists())

    def test_existing_directory_is_refused_before_execution(self):
        self.output.mkdir()
        sentinel = self.output / "keep.txt"
        sentinel.write_text("unchanged\n")
        self.env["DRY_RUN"] = "0"
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Output directory already exists", result.stdout)
        self.assertNotIn("CALLER_MUST_NOT_RUN", result.stdout)
        self.assertEqual(sentinel.read_text(), "unchanged\n")

    def test_unsafe_output_path(self):
        self.env["OUT_ROOT"] = str(REPO_ROOT / "results/forbidden-spikein-output")
        self.assert_rejected("OUT_ROOT must be a directory below results2/spikein")

    def test_output_symlink_escape(self):
        link = self.directory / "escape"
        link.symlink_to(REPO_ROOT / "results", target_is_directory=True)
        self.env["OUT_ROOT"] = str(link / "forbidden-spikein-output")
        self.assert_rejected("OUT_ROOT must be a directory below results2/spikein")

    def test_output_whitespace_is_rejected_before_java(self):
        self.env["OUT_ROOT"] = str(self.directory / "space folder" / "discovery")
        self.assert_rejected("OUT_ROOT must not contain whitespace")


if __name__ == "__main__":
    for tool in ("bash", "samtools", "bcftools"):
        if not shutil.which(tool):
            raise SystemExit(f"Missing {tool}; activate the WSL prnp-spikein environment")
    unittest.main(verbosity=2)
