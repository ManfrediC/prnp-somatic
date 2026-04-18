#!/usr/bin/env python3
"""Manual PRNP ORR mosaic review helper for one BAM/CRAM sample.

This tool is intentionally conservative:
- it uses the existing genomic alignment as the primary evidence layer
- it only inspects primary, mapped alignments over PRNP ORR
- it requires left and right anchor support for high-confidence block calls
- it keeps one-sided reads in the review table, but down-ranks them

The synthetic-reference component is implemented as rescoring extracted ORR
segments against a small local allele panel rather than remapping the full BAM.
That keeps the method auditable and fast for manual review.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import pysam


BLOCKS = {
    "R1": "CCTCAGGGCGGTGGTGGCTGGGGGCAG",
    "R2": "CCTCATGGTGGTGGCTGGGGGCAG",
    "R2a": "CCTCATGGTGGTGGCTGGGGACAG",
    "R2c": "CCTCATGGCGGTGGCTGGGGGCAG",
    "R3": "CCCCATGGTGGTGGCTGGGGACAG",
    "R3g": "CCCCATGGTGGTGGCTGGGGGCAG",
    "R4": "CCTCATGGTGGTGGCTGGGGTCAA",
}
REFERENCE_ARCHITECTURE = "R1-R2-R2-R3-R4"
REFERENCE_ORR_SEQUENCE = "".join(BLOCKS[token] for token in REFERENCE_ARCHITECTURE.split("-"))
REFERENCE_ORR_LENGTH = len(REFERENCE_ORR_SEQUENCE)

CIGAR_MATCH = {0, 7, 8}
CIGAR_INSERTION = 1
CIGAR_DELETION = 2
CIGAR_REF_SKIP = 3
CIGAR_SOFT_CLIP = 4


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build per-read and per-sample TSVs for conservative manual review "
            "of PRNP ORR mosaic OPRI/OPRD evidence."
        )
    )
    parser.add_argument("--bam", required=True, type=Path, help="Input BAM or CRAM.")
    parser.add_argument(
        "--reference-fasta",
        required=True,
        type=Path,
        help="Reference FASTA used by the alignment.",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        type=Path,
        help="Prefix for output TSVs, e.g. results/repeats/manual/CJD1",
    )
    parser.add_argument(
        "--sample-id",
        default="",
        help="Sample identifier. Defaults to the BAM/CRAM basename without suffix.",
    )
    parser.add_argument(
        "--region-bed",
        default=Path("resources/prnp_orr.hg38.bed"),
        type=Path,
        help="BED with the full PRNP ORR interval.",
    )
    parser.add_argument(
        "--panel-tsv",
        default=Path("resources/repeats/prnp_orr_manual_panel.tsv"),
        type=Path,
        help="Synthetic PRNP ORR allele panel.",
    )
    parser.add_argument(
        "--min-mapq",
        default=20,
        type=int,
        help="Minimum mapping quality for review reads.",
    )
    parser.add_argument(
        "--anchor-probe",
        default=12,
        type=int,
        help="Number of reference bases at each ORR edge used to anchor a read.",
    )
    parser.add_argument(
        "--min-anchor-bases",
        default=8,
        type=int,
        help="Minimum anchored bases required per side.",
    )
    parser.add_argument(
        "--nearby-window",
        default=6,
        type=int,
        help="Reference window around ORR used to flag nearby indels/soft clips.",
    )
    parser.add_argument(
        "--sample-calls-tsv",
        default=None,
        type=Path,
        help="Optional results/repeats/sample_calls.tsv for context annotations.",
    )
    parser.add_argument(
        "--gangstr-calls-tsv",
        default=None,
        type=Path,
        help="Optional results/repeats/gangstr_calls.tsv for context annotations.",
    )
    parser.add_argument(
        "--synthetic-remap-mode",
        default="candidate_only",
        choices=("off", "candidate_only", "all_two_sided_exact"),
        help=(
            "Whether to remap reads to local synthetic PRNP ORR contigs. "
            "'candidate_only' remaps only exact two-sided reads with nonreference "
            "length or nearby indel/soft-clip evidence."
        ),
    )
    parser.add_argument(
        "--synthetic-flank-size",
        default=40,
        type=int,
        help="Reference flank size to add on each side of synthetic contigs.",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Small sequence / file helpers
# ---------------------------------------------------------------------------

def reverse_complement(sequence: str) -> str:
    # Synthetic rescoring compares both orientations, so keep reverse-complement
    # handling in one shared helper.
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_bed_interval(path: Path) -> tuple[str, int, int, str]:
    # The workflow expects a single PRNP ORR interval in the BED file.
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else f"{chrom}:{start + 1}-{end}"
            return chrom, start, end, name
    raise ValueError(f"No interval found in BED: {path}")


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    # Optional cohort context tables are plain TSVs, so keep loading simple.
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    # Create the output directory on demand so one-sample reruns can target a
    # clean prefix without precreating parent folders.
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def architecture_to_sequence(architecture: str) -> str:
    # Convert panel architecture labels into literal ORR sequence strings once
    # so later scoring can compare reads directly against sequence.
    return "".join(BLOCKS[token] for token in architecture.split("-"))


def load_panel(path: Path) -> list[dict[str, str]]:
    # Enrich the panel rows with concrete sequence and length-shift fields so
    # downstream scoring never has to recalculate them per read.
    rows = read_tsv_rows(path)
    for row in rows:
        row["sequence"] = architecture_to_sequence(row["architecture"])
        row["delta_bp"] = str(len(row["sequence"]) - REFERENCE_ORR_LENGTH)
    return rows


def build_panel_buckets(panel_rows: list[dict[str, str]]) -> dict[int, list[dict[str, str]]]:
    # Bucket panel alleles by net length shift so Hamming comparisons only
    # occur between architectures with the same overall size.
    buckets: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in panel_rows:
        buckets[int(row["delta_bp"])].append(row)
    return buckets


def build_synthetic_contigs(
    panel_rows: list[dict[str, str]],
    left_flank: str,
    right_flank: str,
) -> list[dict[str, str]]:
    # Synthetic contigs reuse the native genomic flanks so the local remap step
    # is judged in approximately the same sequence context as the BAM read.
    contigs: list[dict[str, str]] = []
    for row in panel_rows:
        contigs.append(
            {
                **row,
                "contig_sequence": left_flank + row["sequence"] + right_flank,
            }
        )
    return contigs


def infer_sample_id(args: argparse.Namespace) -> str:
    # Fall back to the BAM basename for quick ad hoc local runs.
    if args.sample_id:
        return args.sample_id
    return args.bam.stem


def normalize_qpos(read: pysam.AlignedSegment, qpos: int) -> int:
    # Normalize query coordinates to the forward orientation so anchor logic
    # can treat plus and minus reads the same way.
    if not read.is_reverse:
        return qpos
    return read.query_length - 1 - qpos


def normalized_query_sequence(read: pysam.AlignedSegment) -> str:
    # Return sequence in the same orientation used by normalized query indices.
    sequence = read.query_sequence or ""
    return reverse_complement(sequence) if read.is_reverse else sequence


def hamming_distance(left: str, right: str) -> int:
    # Hamming scoring is only used after matching sequences by length.
    if len(left) != len(right):
        raise ValueError("Hamming distance requires equal-length sequences.")
    return sum(1 for lhs, rhs in zip(left, right) if lhs != rhs)


def overlaps(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    # Shared interval helper for nearby-indel and nearby-softclip checks.
    return start_a < end_b and start_b < end_a


# ---------------------------------------------------------------------------
# Read-level evidence extraction from the genomic BAM alignment
# ---------------------------------------------------------------------------

def extract_anchor_segment(
    read: pysam.AlignedSegment,
    orr_start: int,
    orr_end: int,
    anchor_probe: int,
    min_anchor_bases: int,
) -> dict[str, object]:
    # Track how much of each ORR edge is covered in normalized query space.
    # That lets the downstream logic distinguish reads that truly span the
    # repeat from reads that only clip into one side of it.
    sequence = normalized_query_sequence(read)
    aligned_pairs = read.get_aligned_pairs(matches_only=False)

    left_anchor_qnorm: list[int] = []
    right_anchor_qnorm: list[int] = []
    left_exact_qnorm: list[int] = []
    right_exact_qnorm: list[int] = []
    for qpos, rpos in aligned_pairs:
        if qpos is None or rpos is None:
            continue
        qnorm = normalize_qpos(read, qpos)
        if orr_start <= rpos < orr_start + anchor_probe:
            left_anchor_qnorm.append(qnorm)
        if orr_end - anchor_probe <= rpos < orr_end:
            right_anchor_qnorm.append(qnorm)
        if rpos == orr_start:
            left_exact_qnorm.append(qnorm)
        if rpos == orr_end - 1:
            right_exact_qnorm.append(qnorm)

    left_anchor_bases = len(set(left_anchor_qnorm))
    right_anchor_bases = len(set(right_anchor_qnorm))
    has_left_anchor = left_anchor_bases >= min_anchor_bases
    has_right_anchor = right_anchor_bases >= min_anchor_bases

    # Prefer exact edge-to-edge extraction when the read touches both ORR
    # boundaries exactly. Fall back to a wider window only when the read spans
    # the locus but does not align cleanly to both exact endpoints.
    anchor_tier = "unanchored"
    middle_sequence = ""
    if has_left_anchor and has_right_anchor:
        if left_exact_qnorm and right_exact_qnorm:
            qstart = min(left_exact_qnorm)
            qend = max(right_exact_qnorm) + 1
            if qend >= qstart:
                middle_sequence = sequence[qstart:qend]
                anchor_tier = "two_sided_exact"
        else:
            qstart = min(left_anchor_qnorm)
            qend = max(right_anchor_qnorm) + 1
            if qend >= qstart:
                middle_sequence = sequence[qstart:qend]
                anchor_tier = "two_sided_window"
    elif has_left_anchor:
        qstart = min(left_anchor_qnorm)
        middle_sequence = sequence[qstart:]
        anchor_tier = "left_only"
    elif has_right_anchor:
        qend = max(right_anchor_qnorm) + 1
        middle_sequence = sequence[:qend]
        anchor_tier = "right_only"

    return {
        "anchor_tier": anchor_tier,
        "left_anchor_bases": left_anchor_bases,
        "right_anchor_bases": right_anchor_bases,
        "middle_sequence": middle_sequence,
    }


def inspect_cigar_near_orr(
    read: pysam.AlignedSegment,
    orr_start: int,
    orr_end: int,
    nearby_window: int,
) -> dict[str, object]:
    if read.cigartuples is None:
        return {
            "has_indel_near_orr": False,
            "has_softclip_near_orr": False,
            "indel_ops": "",
        }

    # Keep a simple local record of indels and edge soft clips near the ORR.
    # These events are weak evidence, but they are useful for deciding which
    # reads deserve synthetic rescoring even when the extracted length is still
    # reference-like.
    ref_pos = read.reference_start
    has_indel_near_orr = False
    has_softclip_near_orr = False
    indel_ops: list[str] = []

    for index, (operation, length) in enumerate(read.cigartuples):
        if operation in CIGAR_MATCH:
            ref_pos += length
            continue

        if operation == CIGAR_INSERTION:
            if orr_start - nearby_window <= ref_pos <= orr_end + nearby_window:
                has_indel_near_orr = True
                indel_ops.append(f"I{length}@{ref_pos + 1}")
            continue

        if operation in {CIGAR_DELETION, CIGAR_REF_SKIP}:
            op_start = ref_pos
            op_end = ref_pos + length
            if overlaps(op_start, op_end, orr_start - nearby_window, orr_end + nearby_window):
                has_indel_near_orr = True
                indel_ops.append(f"D{length}@{op_start + 1}")
            ref_pos = op_end
            continue

        if operation == CIGAR_SOFT_CLIP:
            at_left_edge = index == 0 and abs(read.reference_start - orr_start) <= nearby_window
            at_right_edge = (
                index == len(read.cigartuples) - 1
                and abs(read.reference_end - orr_end) <= nearby_window
            )
            if at_left_edge or at_right_edge:
                has_softclip_near_orr = True

    return {
        "has_indel_near_orr": has_indel_near_orr,
        "has_softclip_near_orr": has_softclip_near_orr,
        "indel_ops": ",".join(indel_ops),
    }


def classify_bam_signal(anchor_tier: str, middle_sequence: str) -> tuple[str, str, str]:
    if anchor_tier != "two_sided_exact" or not middle_sequence:
        return "", "", "partial_or_unanchored"

    # Convert exact extracted sequence length into a coarse genomic signal first.
    # Panel scoring happens later; this step just says whether the observed
    # sequence is reference-length, whole-repeat shifted, or more irregular.
    delta_bp = len(middle_sequence) - REFERENCE_ORR_LENGTH
    delta_bp_str = str(delta_bp)
    if delta_bp == 0:
        return delta_bp_str, "0", "reference_length"
    if delta_bp % 24 == 0:
        return delta_bp_str, str(delta_bp // 24), "expansion_like" if delta_bp > 0 else "contraction_like"
    return delta_bp_str, "", "complex_length"


def score_against_panel(
    middle_sequence: str,
    delta_bp: str,
    panel_buckets: dict[int, list[dict[str, str]]],
) -> dict[str, str]:
    # Restrict comparisons to panel alleles with the same overall length shift.
    # This keeps the Hamming-distance step interpretable instead of letting very
    # different copy-number models dominate purely by edit distance.
    if not middle_sequence or delta_bp == "":
        return {}

    try:
        bucket = panel_buckets[int(delta_bp)]
    except (ValueError, KeyError):
        return {}

    scores: list[tuple[int, str, dict[str, str]]] = []
    for allele in bucket:
        sequence = allele["sequence"]
        if len(sequence) != len(middle_sequence):
            continue
        scores.append((hamming_distance(middle_sequence, sequence), allele["allele_id"], allele))

    if not scores:
        return {}

    scores.sort(key=lambda item: (item[0], item[1]))
    best_distance, _, best = scores[0]
    second_distance = scores[1][0] if len(scores) > 1 else ""
    distance_gap = (
        str(int(second_distance) - int(best_distance))
        if second_distance != ""
        else ""
    )

    confidence = "low"
    if best_distance == 0 and (second_distance == "" or int(distance_gap) >= 1):
        confidence = "high"
    elif best_distance <= 1 and (second_distance == "" or int(distance_gap or "0") >= 1):
        confidence = "medium"
    elif second_distance != "" and int(second_distance) == int(best_distance):
        confidence = "ambiguous"

    return {
        "panel_scored": "yes",
        "best_allele_id": best["allele_id"],
        "best_event_type": best["event_type"],
        "best_architecture": best["architecture"],
        "best_evidence_tier": best["evidence_tier"],
        "best_source_label": best["source_label"],
        "best_source_url": best["source_url"],
        "best_hamming_distance": str(best_distance),
        "second_best_hamming_distance": str(second_distance),
        "best_minus_second_hamming": distance_gap,
        "panel_confidence": confidence,
    }


def fitting_edit_alignment(query: str, target: str) -> dict[str, int]:
    """Fit the full query to the best substring of target by edit distance."""

    # This is a lightweight fitting-alignment implementation for short extracted
    # read segments. It avoids invoking a full external aligner while still
    # allowing the query to land on the best local window inside each synthetic
    # contig.
    m = len(query)
    n = len(target)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    trace = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(1, m + 1):
        dp[i][0] = i
        trace[i][0] = 1  # up
    for j in range(1, n + 1):
        dp[0][j] = 0
        trace[0][j] = 2  # left

    for i in range(1, m + 1):
        qbase = query[i - 1]
        for j in range(1, n + 1):
            cost_diag = dp[i - 1][j - 1] + (0 if qbase == target[j - 1] else 1)
            cost_up = dp[i - 1][j] + 1
            cost_left = dp[i][j - 1] + 1
            best_cost = cost_diag
            best_trace = 0
            if cost_up < best_cost:
                best_cost = cost_up
                best_trace = 1
            if cost_left < best_cost:
                best_cost = cost_left
                best_trace = 2
            dp[i][j] = best_cost
            trace[i][j] = best_trace

    best_end = min(range(n + 1), key=lambda j: (dp[m][j], j))
    edit_distance = dp[m][best_end]

    i = m
    j = best_end
    aligned_columns = 0
    while i > 0:
        direction = trace[i][j]
        if direction == 0:
            i -= 1
            j -= 1
            aligned_columns += 1
        elif direction == 1:
            i -= 1
            aligned_columns += 1
        else:
            j -= 1
            aligned_columns += 1

    start = j
    return {
        "edit_distance": edit_distance,
        "target_start": start,
        "target_end": best_end,
        "aligned_columns": aligned_columns,
    }


def synthetic_remap_read(
    query_sequence: str,
    contigs: list[dict[str, str]],
) -> dict[str, str]:
    # Score both read orientations against every synthetic contig, then keep the
    # best local fit. This is intentionally a conservative rescoring layer, not
    # a replacement for the primary genomic alignment.
    candidates: list[tuple[int, int, str, str, dict[str, str], dict[str, int]]] = []
    for orientation, oriented_query in (("+", query_sequence), ("-", reverse_complement(query_sequence))):
        for contig in contigs:
            alignment = fitting_edit_alignment(oriented_query, contig["contig_sequence"])
            candidates.append(
                (
                    alignment["edit_distance"],
                    -alignment["aligned_columns"],
                    contig["allele_id"],
                    orientation,
                    contig,
                    alignment,
                )
            )

    candidates.sort(key=lambda item: item[:4])
    best = candidates[0]
    second = candidates[1] if len(candidates) > 1 else None
    best_distance = best[0]
    second_distance = second[0] if second is not None else ""
    distance_gap = str(int(second_distance) - int(best_distance)) if second is not None else ""

    confidence = "low"
    if best_distance == 0 and (second is None or int(distance_gap) >= 1):
        confidence = "high"
    elif best_distance <= 2 and (second is None or int(distance_gap or "0") >= 1):
        confidence = "medium"
    elif second is not None and int(second_distance) == int(best_distance):
        confidence = "ambiguous"

    _, _, _, orientation, contig, alignment = best
    return {
        "synthetic_remap_attempted": "yes",
        "synthetic_best_allele_id": contig["allele_id"],
        "synthetic_best_event_type": contig["event_type"],
        "synthetic_best_architecture": contig["architecture"],
        "synthetic_best_evidence_tier": contig["evidence_tier"],
        "synthetic_best_orientation": orientation,
        "synthetic_best_edit_distance": str(alignment["edit_distance"]),
        "synthetic_second_best_edit_distance": str(second_distance),
        "synthetic_edit_distance_gap": distance_gap,
        "synthetic_best_target_start_1based": str(alignment["target_start"] + 1),
        "synthetic_best_target_end_1based": str(alignment["target_end"]),
        "synthetic_aligned_columns": str(alignment["aligned_columns"]),
        "synthetic_remap_confidence": confidence,
    }


def load_sample_context(path: Path | None) -> dict[str, dict[str, str]]:
    # Optional sample-level context tables are keyed once here so the main loop
    # can enrich its summary row without repeated file reads.
    if path is None or not path.exists():
        return {}
    return {row["sample_id"]: row for row in read_tsv_rows(path)}


def should_attempt_synthetic_remap(
    synthetic_remap_mode: str,
    anchor_tier: str,
    bam_signal: str,
    has_indel_near_orr: bool,
    has_softclip_near_orr: bool,
) -> tuple[bool, str]:
    # The default mode avoids remapping every reference-like read. Instead it
    # spends effort on exact spanning reads that already show either a repeat-
    # sized length change or nearby local alignment stress.
    if synthetic_remap_mode == "off":
        return False, "disabled"
    if anchor_tier != "two_sided_exact":
        return False, "not_two_sided_exact"
    if synthetic_remap_mode == "all_two_sided_exact":
        return True, "all_two_sided_exact"
    if bam_signal != "reference_length":
        return True, "nonreference_length"
    if has_indel_near_orr:
        return True, "nearby_indel"
    if has_softclip_near_orr:
        return True, "nearby_softclip"
    return False, "reference_like_exact_read"


def build_read_row(
    sample_id: str,
    read: pysam.AlignedSegment,
    anchor_info: dict[str, object],
    cigar_info: dict[str, object],
    panel_info: dict[str, str],
    synthetic_info: dict[str, str],
    synthetic_reason: str,
) -> dict[str, str]:
    anchor_tier = str(anchor_info["anchor_tier"])
    middle_sequence = str(anchor_info["middle_sequence"])
    delta_bp, delta_repeats, bam_signal = classify_bam_signal(anchor_tier, middle_sequence)

    # Keep the per-read table intentionally wide: the raw BAM-derived signal,
    # panel score, and synthetic-remap score all stay side by side so manual
    # adjudication never has to reconstruct how a row was classified.
    row = {
        "sample_id": sample_id,
        "read_name": read.query_name,
        "is_read1": "yes" if read.is_read1 else "no",
        "is_read2": "yes" if read.is_read2 else "no",
        "strand": "-" if read.is_reverse else "+",
        "chrom": read.reference_name or "",
        "reference_start_1based": str(read.reference_start + 1),
        "reference_end_1based": str(read.reference_end),
        "mapq": str(read.mapping_quality),
        "query_length": str(read.query_length or 0),
        "anchor_tier": anchor_tier,
        "left_anchor_bases": str(anchor_info["left_anchor_bases"]),
        "right_anchor_bases": str(anchor_info["right_anchor_bases"]),
        "orr_sequence_length": str(len(middle_sequence)),
        "delta_bp_vs_reference": delta_bp,
        "delta_repeats_vs_reference": delta_repeats,
        "bam_signal": bam_signal,
        "has_indel_near_orr": "yes" if cigar_info["has_indel_near_orr"] else "no",
        "has_softclip_near_orr": "yes" if cigar_info["has_softclip_near_orr"] else "no",
        "indel_ops": str(cigar_info["indel_ops"]),
        "middle_sequence": middle_sequence,
        "panel_scored": panel_info.get("panel_scored", "no"),
        "best_allele_id": panel_info.get("best_allele_id", ""),
        "best_event_type": panel_info.get("best_event_type", ""),
        "best_architecture": panel_info.get("best_architecture", ""),
        "best_evidence_tier": panel_info.get("best_evidence_tier", ""),
        "best_source_label": panel_info.get("best_source_label", ""),
        "best_source_url": panel_info.get("best_source_url", ""),
        "best_hamming_distance": panel_info.get("best_hamming_distance", ""),
        "second_best_hamming_distance": panel_info.get("second_best_hamming_distance", ""),
        "best_minus_second_hamming": panel_info.get("best_minus_second_hamming", ""),
        "panel_confidence": panel_info.get("panel_confidence", ""),
        "synthetic_remap_attempted": synthetic_info.get("synthetic_remap_attempted", "no"),
        "synthetic_remap_reason": synthetic_reason,
        "synthetic_best_allele_id": synthetic_info.get("synthetic_best_allele_id", ""),
        "synthetic_best_event_type": synthetic_info.get("synthetic_best_event_type", ""),
        "synthetic_best_architecture": synthetic_info.get("synthetic_best_architecture", ""),
        "synthetic_best_evidence_tier": synthetic_info.get("synthetic_best_evidence_tier", ""),
        "synthetic_best_orientation": synthetic_info.get("synthetic_best_orientation", ""),
        "synthetic_best_edit_distance": synthetic_info.get("synthetic_best_edit_distance", ""),
        "synthetic_second_best_edit_distance": synthetic_info.get("synthetic_second_best_edit_distance", ""),
        "synthetic_edit_distance_gap": synthetic_info.get("synthetic_edit_distance_gap", ""),
        "synthetic_best_target_start_1based": synthetic_info.get("synthetic_best_target_start_1based", ""),
        "synthetic_best_target_end_1based": synthetic_info.get("synthetic_best_target_end_1based", ""),
        "synthetic_aligned_columns": synthetic_info.get("synthetic_aligned_columns", ""),
        "synthetic_remap_confidence": synthetic_info.get("synthetic_remap_confidence", ""),
        "synthetic_vs_panel_call": "",
    }

    if row["panel_scored"] != "yes" and bam_signal in {"expansion_like", "contraction_like", "complex_length"}:
        row["panel_confidence"] = "low"
    if row["panel_scored"] == "yes" and row["synthetic_remap_attempted"] == "yes":
        row["synthetic_vs_panel_call"] = (
            "match"
            if row["best_allele_id"] == row["synthetic_best_allele_id"]
            else "mismatch"
        )
    elif row["synthetic_remap_attempted"] == "yes":
        row["synthetic_vs_panel_call"] = "panel_not_scored"

    return row


def summarize_rows(
    sample_id: str,
    bam_path: Path,
    panel_path: Path,
    rows: list[dict[str, str]],
    sample_calls_row: dict[str, str] | None,
    gangstr_row: dict[str, str] | None,
) -> dict[str, str]:
    # Collapse the read table into one row that is easy to sort across a cohort.
    # The summary preserves the counts most useful for separating recurrent
    # background from a plausible mosaic allele.
    anchor_counts = Counter(row["anchor_tier"] for row in rows)
    bam_signal_counts = Counter(row["bam_signal"] for row in rows)
    delta_bp_counts = Counter(
        row["delta_bp_vs_reference"]
        for row in rows
        if row["delta_bp_vs_reference"] != ""
    )
    allele_counts = Counter(
        row["best_allele_id"]
        for row in rows
        if row["best_allele_id"]
    )
    synthetic_allele_counts = Counter(
        row["synthetic_best_allele_id"]
        for row in rows
        if row["synthetic_best_allele_id"]
    )
    synthetic_concordance_counts = Counter(
        row["synthetic_vs_panel_call"]
        for row in rows
        if row["synthetic_vs_panel_call"]
    )

    def encode_counter(counter: Counter[str]) -> str:
        if not counter:
            return ""
        return ",".join(f"{label}:{counter[label]}" for label in sorted(counter))

    high_conf_nonref = sum(
        1
        for row in rows
        if row["best_event_type"] not in {"", "reference"}
        and row["panel_confidence"] in {"high", "medium"}
    )
    synthetic_remapped_reads = sum(
        1 for row in rows if row["synthetic_remap_attempted"] == "yes"
    )
    synthetic_high_conf_nonref = sum(
        1
        for row in rows
        if row["synthetic_best_event_type"] not in {"", "reference"}
        and row["synthetic_remap_confidence"] in {"high", "medium"}
    )
    one_sided_nonref = sum(
        1
        for row in rows
        if row["anchor_tier"] in {"left_only", "right_only"}
        and (
            row["has_indel_near_orr"] == "yes"
            or row["has_softclip_near_orr"] == "yes"
        )
    )

    summary = {
        "sample_id": sample_id,
        "bam_path": str(bam_path),
        "panel_path": str(panel_path),
        "reference_architecture": REFERENCE_ARCHITECTURE,
        "total_overlapping_primary_reads": str(len(rows)),
        "two_sided_exact_reads": str(anchor_counts.get("two_sided_exact", 0)),
        "two_sided_window_reads": str(anchor_counts.get("two_sided_window", 0)),
        "left_only_reads": str(anchor_counts.get("left_only", 0)),
        "right_only_reads": str(anchor_counts.get("right_only", 0)),
        "unanchored_reads": str(anchor_counts.get("unanchored", 0)),
        "bam_signal_counts": encode_counter(bam_signal_counts),
        "delta_bp_counts": encode_counter(delta_bp_counts),
        "best_allele_counts": encode_counter(allele_counts),
        "synthetic_best_allele_counts": encode_counter(synthetic_allele_counts),
        "synthetic_vs_panel_counts": encode_counter(synthetic_concordance_counts),
        "high_or_medium_confidence_nonreference_reads": str(high_conf_nonref),
        "synthetic_remapped_reads": str(synthetic_remapped_reads),
        "synthetic_high_or_medium_confidence_nonreference_reads": str(synthetic_high_conf_nonref),
        "one_sided_indel_or_softclip_reads": str(one_sided_nonref),
        "eh_interpretation": "",
        "eh_total_repeat_counts": "",
        "gangstr_interpretation": "",
        "gangstr_total_repeat_counts": "",
    }

    # Propagate high-level caller context into the manual-review summary so the
    # cohort wrapper can compare evidence layers without another join step.
    if sample_calls_row is not None:
        summary["eh_interpretation"] = sample_calls_row.get("interpretation", "")
        summary["eh_total_repeat_counts"] = sample_calls_row.get("total_repeat_counts", "")
    if gangstr_row is not None:
        summary["gangstr_interpretation"] = gangstr_row.get("gangstr_interpretation", "")
        summary["gangstr_total_repeat_counts"] = gangstr_row.get("total_repeat_counts", "")

    return summary


def sort_read_rows(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    def sort_key(row: dict[str, str]) -> tuple[int, int, int, str]:
        # Surface the rows most worth eyeballing first: nonreference-length
        # exact reads before partial anchors, then better panel matches.
        interesting = 0 if row["bam_signal"] in {"expansion_like", "contraction_like", "complex_length"} else 1
        anchor_rank = {
            "two_sided_exact": 0,
            "two_sided_window": 1,
            "left_only": 2,
            "right_only": 2,
            "unanchored": 3,
        }.get(row["anchor_tier"], 3)
        best_distance = int(row["best_hamming_distance"]) if row["best_hamming_distance"] else 9999
        return (interesting, anchor_rank, best_distance, row["read_name"])

    return sorted(rows, key=sort_key)


def main() -> int:
    args = parse_args()
    sample_id = infer_sample_id(args)

    # Validate that the FASTA sequence at the configured ORR interval still
    # matches the expected canonical PRNP architecture before scoring any reads.
    chrom, orr_start, orr_end, _ = read_bed_interval(args.region_bed)
    panel_rows = load_panel(args.panel_tsv)
    panel_buckets = build_panel_buckets(panel_rows)

    with pysam.FastaFile(str(args.reference_fasta)) as fasta_handle:
        reference_orr = fasta_handle.fetch(chrom, orr_start, orr_end).upper()
        left_flank = fasta_handle.fetch(
            chrom, max(0, orr_start - args.synthetic_flank_size), orr_start
        ).upper()
        right_flank = fasta_handle.fetch(
            chrom, orr_end, orr_end + args.synthetic_flank_size
        ).upper()
    if reference_orr != REFERENCE_ORR_SEQUENCE:
        raise ValueError(
            "Reference FASTA PRNP ORR sequence does not match the expected "
            "R1-R2-R2-R3-R4 architecture."
        )
    synthetic_contigs = build_synthetic_contigs(panel_rows, left_flank, right_flank)

    # Load optional cohort-level context so the one-sample summary can carry
    # EH and GangSTR labels alongside the BAM-derived manual evidence.
    sample_calls = load_sample_context(args.sample_calls_tsv)
    gangstr_calls = load_sample_context(args.gangstr_calls_tsv)
    sample_calls_row = sample_calls.get(sample_id)
    gangstr_row = gangstr_calls.get(sample_id)

    read_rows: list[dict[str, str]] = []
    with pysam.AlignmentFile(str(args.bam), reference_filename=str(args.reference_fasta)) as bam_handle:
        fetch_start = max(0, orr_start - args.anchor_probe)
        fetch_end = orr_end + args.anchor_probe
        for read in bam_handle.fetch(chrom, fetch_start, fetch_end):
            # This helper is meant for interpretable manual review, so keep only
            # primary mapped reads with reasonable mapping quality.
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate or read.is_qcfail:
                continue
            if read.mapping_quality < args.min_mapq:
                continue
            if read.reference_name != chrom:
                continue

            anchor_info = extract_anchor_segment(
                read=read,
                orr_start=orr_start,
                orr_end=orr_end,
                anchor_probe=args.anchor_probe,
                min_anchor_bases=args.min_anchor_bases,
            )
            cigar_info = inspect_cigar_near_orr(
                read=read,
                orr_start=orr_start,
                orr_end=orr_end,
                nearby_window=args.nearby_window,
            )

            panel_info: dict[str, str] = {}
            if anchor_info["anchor_tier"] == "two_sided_exact":
                # Only exact two-sided reads are panel-scored. Partial anchors
                # stay in the table, but they are too ambiguous for direct
                # sequence-to-allele Hamming comparisons.
                delta_bp, _, _ = classify_bam_signal(
                    str(anchor_info["anchor_tier"]),
                    str(anchor_info["middle_sequence"]),
                )
                panel_info = score_against_panel(
                    middle_sequence=str(anchor_info["middle_sequence"]),
                    delta_bp=delta_bp,
                    panel_buckets=panel_buckets,
                )

            synthetic_info: dict[str, str] = {}
            # Synthetic rescoring is gated separately so we can prioritize only
            # the reads most likely to clarify a possible mosaic signal.
            synthetic_attempt, synthetic_reason = should_attempt_synthetic_remap(
                synthetic_remap_mode=args.synthetic_remap_mode,
                anchor_tier=str(anchor_info["anchor_tier"]),
                bam_signal=classify_bam_signal(
                    str(anchor_info["anchor_tier"]),
                    str(anchor_info["middle_sequence"]),
                )[2],
                has_indel_near_orr=bool(cigar_info["has_indel_near_orr"]),
                has_softclip_near_orr=bool(cigar_info["has_softclip_near_orr"]),
            )
            if synthetic_attempt:
                synthetic_info = synthetic_remap_read(
                    query_sequence=read.query_sequence or "",
                    contigs=synthetic_contigs,
                )

            read_rows.append(
                build_read_row(
                    sample_id=sample_id,
                    read=read,
                    anchor_info=anchor_info,
                    cigar_info=cigar_info,
                    panel_info=panel_info,
                    synthetic_info=synthetic_info,
                    synthetic_reason=synthetic_reason,
                )
            )

    read_rows = sort_read_rows(read_rows)
    summary_row = summarize_rows(
        sample_id=sample_id,
        bam_path=args.bam,
        panel_path=args.panel_tsv,
        rows=read_rows,
        sample_calls_row=sample_calls_row,
        gangstr_row=gangstr_row,
    )

    # Keep the output schemas explicit so cohort aggregation and ad hoc review
    # notebooks can rely on stable column names.
    read_fieldnames = [
        "sample_id",
        "read_name",
        "is_read1",
        "is_read2",
        "strand",
        "chrom",
        "reference_start_1based",
        "reference_end_1based",
        "mapq",
        "query_length",
        "anchor_tier",
        "left_anchor_bases",
        "right_anchor_bases",
        "orr_sequence_length",
        "delta_bp_vs_reference",
        "delta_repeats_vs_reference",
        "bam_signal",
        "has_indel_near_orr",
        "has_softclip_near_orr",
        "indel_ops",
        "middle_sequence",
        "panel_scored",
        "best_allele_id",
        "best_event_type",
        "best_architecture",
        "best_evidence_tier",
        "best_source_label",
        "best_source_url",
        "best_hamming_distance",
        "second_best_hamming_distance",
        "best_minus_second_hamming",
        "panel_confidence",
        "synthetic_remap_attempted",
        "synthetic_remap_reason",
        "synthetic_best_allele_id",
        "synthetic_best_event_type",
        "synthetic_best_architecture",
        "synthetic_best_evidence_tier",
        "synthetic_best_orientation",
        "synthetic_best_edit_distance",
        "synthetic_second_best_edit_distance",
        "synthetic_edit_distance_gap",
        "synthetic_best_target_start_1based",
        "synthetic_best_target_end_1based",
        "synthetic_aligned_columns",
        "synthetic_remap_confidence",
        "synthetic_vs_panel_call",
    ]
    summary_fieldnames = [
        "sample_id",
        "bam_path",
        "panel_path",
        "reference_architecture",
        "total_overlapping_primary_reads",
        "two_sided_exact_reads",
        "two_sided_window_reads",
        "left_only_reads",
        "right_only_reads",
        "unanchored_reads",
        "bam_signal_counts",
        "delta_bp_counts",
        "best_allele_counts",
        "synthetic_best_allele_counts",
        "synthetic_vs_panel_counts",
        "high_or_medium_confidence_nonreference_reads",
        "synthetic_remapped_reads",
        "synthetic_high_or_medium_confidence_nonreference_reads",
        "one_sided_indel_or_softclip_reads",
        "eh_interpretation",
        "eh_total_repeat_counts",
        "gangstr_interpretation",
        "gangstr_total_repeat_counts",
    ]

    # Write the detailed read ledger plus the one-row sample summary as separate
    # TSVs so cohort wrappers can aggregate summaries without losing the raw
    # evidence needed for manual follow-up.
    write_tsv(args.output_prefix.with_suffix(".reads.tsv"), read_rows, read_fieldnames)
    write_tsv(
        args.output_prefix.with_suffix(".summary.tsv"),
        [summary_row],
        summary_fieldnames,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
