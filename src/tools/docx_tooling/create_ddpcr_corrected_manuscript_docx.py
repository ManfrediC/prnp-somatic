#!/usr/bin/env python3
"""Create a clean revised DOCX for the ddPCR manuscript update.

The tracked-edit DOCX is produced afterwards with Microsoft Word's Compare
Documents feature, using this clean DOCX as the revised document.
"""

from __future__ import annotations

import copy
import csv
import math
import re
import shutil
import tempfile
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parents[3]
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript" / "v2"
ORIGINAL_DOCX = MANUSCRIPT_DIR / "2026_05_11_Carta_PRNP_v2.docx"
CLEAN_DOCX = MANUSCRIPT_DIR / "2026_05_11_Carta_PRNP_v2.ddpcr_corrected_clean.docx"
SNV_DATA_XLSX = PROJECT_ROOT / "results" / "ddPCR" / "SNV_data_final.xlsx"
PARTICIPANT_POOLED_CSV = (
    PROJECT_ROOT
    / "manuscript"
    / "figures"
    / "ddpcr_fractional_abundance_pooled"
    / "participant_pooled_rows.csv"
)
SAMPLE_REGION_HAPLOID_CSV = (
    PROJECT_ROOT / "results" / "ddPCR" / "haploid_genomes_surveyed" / "sample_region_haploid_genomes.csv"
)

W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
XML_NS = "http://www.w3.org/XML/1998/namespace"
SS_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"

ET.register_namespace("w", W_NS)


def qn(ns: str, tag: str) -> str:
    return f"{{{ns}}}{tag}"


def paragraph_text(p: ET.Element) -> str:
    parts: list[str] = []
    for node in p.iter():
        if node.tag in {qn(W_NS, "t"), qn(W_NS, "delText")}:
            parts.append(node.text or "")
    return "".join(parts).strip()


def text_run(text: str) -> ET.Element:
    run = ET.Element(qn(W_NS, "r"))
    t = ET.SubElement(run, qn(W_NS, "t"))
    if text.startswith(" ") or text.endswith(" "):
        t.set(qn(XML_NS, "space"), "preserve")
    t.text = text
    return run


def formula_paragraph(text: str) -> ET.Element:
    p = ET.Element(qn(W_NS, "p"))
    ppr = ET.SubElement(p, qn(W_NS, "pPr"))
    jc = ET.SubElement(ppr, qn(W_NS, "jc"))
    jc.set(qn(W_NS, "val"), "center")
    p.append(text_run(text))
    return p


FORMULAS = {
    "lambda": "λ_t = -ln(n_t- / N)",
    "fa": "FA (%) = 100 × λ_mut / (λ_mut + λ_ref)",
    "lob": "LoB_count = Q_0.95[Binomial(N, p_0)]",
    "gtotal": "G_total = -N ln(N_0,0 / N)",
}


REGION_LABELS = {
    "bg": "basal ganglia",
    "cb": "cerebellum",
    "fr": "frontal cortex",
    "hc": "hippocampus",
    "ps": "pons",
    "sn": "substantia nigra",
    "th": "thalamus",
    "pooled": "pooled",
}


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def column_index(cell_ref: str) -> int:
    letters = re.match(r"[A-Z]+", cell_ref)
    if letters is None:
        raise ValueError(f"Could not parse Excel cell reference: {cell_ref}")
    idx = 0
    for char in letters.group(0):
        idx = idx * 26 + ord(char) - ord("A") + 1
    return idx - 1


def read_xlsx_rows(path: Path) -> list[dict[str, object]]:
    with zipfile.ZipFile(path, "r") as zf:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            shared_root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for item in shared_root.findall(f".//{{{SS_NS}}}si"):
                shared_strings.append("".join(t.text or "" for t in item.findall(f".//{{{SS_NS}}}t")))

        workbook_root = ET.fromstring(zf.read("xl/workbook.xml"))
        first_sheet = workbook_root.find(f".//{{{SS_NS}}}sheet")
        if first_sheet is None:
            raise RuntimeError(f"No sheets found in {path}")
        rel_id = first_sheet.get(qn(REL_NS, "id"))
        if not rel_id:
            raise RuntimeError(f"First sheet in {path} has no relationship id")

        rels_root = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        target = None
        for rel in rels_root.findall(f".//{{{PKG_REL_NS}}}Relationship"):
            if rel.get("Id") == rel_id:
                target = rel.get("Target")
                break
        if target is None:
            raise RuntimeError(f"Could not resolve first worksheet in {path}")
        target = target.lstrip("/")
        sheet_path = target if target.startswith("xl/") else "xl/" + target
        sheet_root = ET.fromstring(zf.read(sheet_path))

    table: list[list[object]] = []
    for row in sheet_root.findall(f".//{{{SS_NS}}}row"):
        values: list[object] = []
        for cell in row.findall(f"{{{SS_NS}}}c"):
            ref = cell.get("r", "")
            idx = column_index(ref)
            while len(values) <= idx:
                values.append("")

            cell_type = cell.get("t")
            value_node = cell.find(f"{{{SS_NS}}}v")
            if cell_type == "inlineStr":
                value = "".join(t.text or "" for t in cell.findall(f".//{{{SS_NS}}}t"))
            elif value_node is None or value_node.text is None:
                value = ""
            elif cell_type == "s":
                value = shared_strings[int(value_node.text)]
            elif cell_type == "b":
                value = value_node.text == "1"
            else:
                raw = value_node.text
                try:
                    value = float(raw)
                    if value.is_integer():
                        value = int(value)
                except ValueError:
                    value = raw
            values[idx] = value
        table.append(values)

    if not table:
        return []
    headers = [str(value) for value in table[0]]
    rows: list[dict[str, object]] = []
    for raw_row in table[1:]:
        row = {header: raw_row[i] if i < len(raw_row) else "" for i, header in enumerate(headers)}
        rows.append(row)
    return rows


def as_true(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def as_float(value: object) -> float:
    if value in {None, "", "NA"}:
        return math.nan
    return float(value)


def format_decimal(value: object, digits: int = 3) -> str:
    number = as_float(value)
    return f"{number:.{digits}f}"


def format_integer(value: float) -> str:
    return f"{round(value):,}"


def participant_key(value: object) -> tuple[int, int, str]:
    text = str(value)
    match = re.search(r"\d+", text)
    number = int(match.group(0)) if match else 10_000
    if text.startswith("CJD"):
        group = 0
    elif text.startswith("Control"):
        group = 1
    else:
        group = 2
    return group, number, text


def join_english(items: list[str]) -> str:
    if not items:
        return "none"
    if len(items) == 1:
        return items[0]
    if len(items) == 2:
        return f"{items[0]} and {items[1]}"
    return f"{', '.join(items[:-1])}, and {items[-1]}"


def format_sample_call(row: dict[str, object]) -> str:
    region = REGION_LABELS.get(str(row["brain_region"]), str(row["brain_region"]))
    return (
        f"{row['participant']} {region} ({row['mutation']}; FA {format_decimal(row['fractional_abundance'])}%, "
        f"95% CI {format_decimal(row['ci_low'])}-{format_decimal(row['ci_high'])}%)"
    )


def format_participant_call(row: dict[str, object]) -> str:
    return (
        f"{row['participant']} {row['mutation']} (FA {format_decimal(row['fractional_abundance'])}%, "
        f"95% CI {format_decimal(row['ci_low'])}-{format_decimal(row['ci_high'])}%)"
    )


def build_ddpcr_summary() -> dict[str, object]:
    snv_rows = read_xlsx_rows(SNV_DATA_XLSX)
    participant_rows = read_csv_rows(PARTICIPANT_POOLED_CSV)
    haploid_rows = read_csv_rows(SAMPLE_REGION_HAPLOID_CSV)

    def is_cjd30_e200k(row: dict[str, object]) -> bool:
        return str(row.get("participant")) == "CJD30" and str(row.get("mutation")) == "E200K"

    non_germline_rows = [row for row in snv_rows if not is_cjd30_e200k(row)]
    sample_pass = [
        row
        for row in non_germline_rows
        if as_true(row.get("detected_above_LoB")) and as_true(row.get("detected_above_LoD"))
    ]
    non_pooled_pass = [row for row in sample_pass if not as_true(row.get("is_pooled"))]
    technical_pooled_pass = [row for row in sample_pass if as_true(row.get("is_pooled"))]
    e200k_technical_pooled_pass = [
        row
        for row in technical_pooled_pass
        if str(row.get("mutation")) == "E200K"
    ]
    e200k_technical_pooled_pass.sort(key=lambda row: participant_key(row.get("participant")))

    participant_pass = [
        row
        for row in participant_rows
        if not is_cjd30_e200k(row)
        and as_true(row.get("detected_above_lob"))
        and as_true(row.get("detected_above_lod"))
    ]
    participant_pass.sort(key=lambda row: participant_key(row.get("participant")))

    participant_lookup = {
        (str(row.get("participant")), str(row.get("mutation"))): row
        for row in participant_rows
    }
    pooled_follow_up_rows = [
        participant_lookup[(str(row["participant"]), str(row["mutation"]))]
        for row in e200k_technical_pooled_pass
        if (str(row["participant"]), str(row["mutation"])) in participant_lookup
        and not (
            as_true(participant_lookup[(str(row["participant"]), str(row["mutation"]))].get("detected_above_lob"))
            and as_true(participant_lookup[(str(row["participant"]), str(row["mutation"]))].get("detected_above_lod"))
        )
    ]

    lob_only_rows = [
        row
        for row in non_germline_rows
        if as_true(row.get("detected_above_LoB")) and not as_true(row.get("detected_above_LoD"))
    ]
    cjd30_e200k_rows = [row for row in snv_rows if is_cjd30_e200k(row)]
    mm2t_participants = sorted(
        {str(row.get("participant")) for row in snv_rows if str(row.get("histotype")) == "MM2T"},
        key=participant_key,
    )

    accepted_counts = [as_float(row["n_accepted_droplets"]) for row in haploid_rows]
    haploid_counts = [as_float(row["n_haploid_genomes_poisson"]) for row in haploid_rows]
    cjd30_fa_values = [as_float(row["fractional_abundance"]) for row in cjd30_e200k_rows]

    return {
        "sample_pass_count": len(sample_pass),
        "non_pooled_pass_count": len(non_pooled_pass),
        "technical_pooled_pass_count": len(technical_pooled_pass),
        "technical_pooled_pass_description": join_english(
            [format_sample_call(row) for row in e200k_technical_pooled_pass]
        ),
        "participant_pass_count": len(participant_pass),
        "participant_pass_description": join_english(
            [format_participant_call(row) for row in participant_pass]
        ),
        "pooled_follow_up_description": join_english(
            [
                f"{row['participant']} pooled FA {format_decimal(row['fractional_abundance'])}%"
                for row in pooled_follow_up_rows
            ]
        ),
        "lob_only_count": len(lob_only_rows),
        "mm2t_count": len(mm2t_participants),
        "accepted_min": format_integer(min(accepted_counts)),
        "accepted_max": format_integer(max(accepted_counts)),
        "haploid_min": format_integer(min(haploid_counts)),
        "haploid_max": format_integer(max(haploid_counts)),
        "cjd30_fa_min": f"{min(cjd30_fa_values):.1f}",
        "cjd30_fa_max": f"{max(cjd30_fa_values):.1f}",
    }


DDPCR_SUMMARY = build_ddpcr_summary()


REPLACEMENTS = {
    "Somatic mutation of the prion protein gene PRNP has been proposed": (
        "Somatic mutation of the prion protein gene PRNP has been proposed as a cause of sporadic "
        "Creutzfeldt-Jakob disease (sCJD), but existing evidence is based on limited anatomical "
        "sampling and incomplete coverage of PRNP. To test this model, we developed a "
        "decontamination procedure for DNA extraction from prion-infected brain tissue and combined "
        "high-sensitivity droplet digital PCR (ddPCR) across multiple brain regions with "
        "probe-based capture sequencing of the entire PRNP locus. Using validated ddPCR assays with "
        "empirical limits of detection of 0.056-0.13% variant allele fraction (VAF), we screened "
        "bulk DNA from six brain regions in 31 putative CJD cases and 8 non-prion controls for "
        "D178N, E200K, and P102L. After excluding an incidentally identified E200K germline "
        "carrier, no non-pooled sample-region measurement exceeded both the limit of blank and "
        "limit of detection. Two technically pooled sample-region E200K summaries exceeded both "
        f"criteria ({DDPCR_SUMMARY['technical_pooled_pass_description']}), and "
        f"{DDPCR_SUMMARY['participant_pass_description']} remained positive after participant-level "
        "pooling. No D178N or P102L measurement met both criteria, and no D178N signal was detected "
        "in sporadic fatal insomnia (thalamic CJD). Targeted sequencing identified no somatic variants in the PRNP "
        "coding region. Two rare intronic candidate variants were observed at low VAFs, one "
        "narrowly exceeding the sequencing limit of detection. Together, these findings argue "
        "against common, recurrent tissue-level mosaicism for canonical pathogenic PRNP missense "
        "variants as a trigger for sCJD in the brain regions studied."
    ),
    "Sample preparation: The ddPCR reactions were prepared": (
        "Sample preparation: The ddPCR reactions were prepared in a dedicated pre-PCR workspace "
        "using template gDNA, fluorescent hydrolysis probes, 2x ddPCR Supermix for Probes "
        "(Bio-Rad) and nuclease-free water. Template gDNA extracted from bulk tissue (130 ng) was "
        "digested with HindIII-HF (New England Biolabs, Ipswich, MA, USA) at 37 °C for 60 min. "
        "Each reaction plate contained PRNP mutant or single-gene controls (SGCs), where "
        "appropriate, together with wild-type (WT) controls and no-template controls (NTCs). WT "
        "controls were prepared with commercially available human gDNA containing a reference PRNP "
        "allele (Human Genomic DNA, mixed; Cat. Nr. G3041; Promega AG, Dübendorf, Switzerland). "
        "NTCs were prepared after completion of all other reaction mixes to control for "
        "cross-contamination. Droplets were produced with a QX200 Droplet Generator (Bio-Rad, "
        "Hercules, CA, USA) and PCRs were performed as published [29]: 10 min at 95 °C; 40 cycles "
        "of 30 s at 94 °C and 60 s at 60 °C; followed by 10 min at 98 °C and cooling to 12 °C. "
        "Droplet fluorescence data were acquired on a QX200 Droplet Reader and initially reviewed "
        "in QuantaSoft Analysis Pro (Bio-Rad). QuantaSoft was used for gating and well-level QC; "
        "downstream fractional-abundance and detection analyses were performed from exported "
        "project data in R. Wells showing inefficient amplification, poor droplet quality, or fewer "
        "than 10,000 accepted droplets were excluded from further analysis."
    ),
    "Single nucleotide variant analysis: To search for pathogenic PRNP": (
        "Single nucleotide variant analysis: To search for pathogenic PRNP single nucleotide "
        "variants (SNVs), we used commercially available 40x TaqMan SNP Assays (Thermo Fisher "
        "Scientific) targeting E200K (Chr20:4699818, Assay ID C__27531205_10) and D178N "
        "(Chr20:4699752, Assay ID C____574563_20). A customised hydrolysis probe assay was "
        "designed for P102L (Custom TaqMan SNP Genotyping Assay, human; Chr20:4699525; forward "
        "primer: 5′-GTGGCACCCACAGTCAGT-3′; reverse primer: 5′-GCAGCACCAGCCATGTG-3′; wild type "
        "probe sequence: VIC-5′-TTGGCTTACTCGGCTTGT-3′-NFQ, mutant probe sequence: "
        "FAM-5′-TGGCTTACTCAGCTTGT-3′-NFQ). DNA was quantified spectrophotometrically. Mutant "
        "controls were prepared with gDNA containing the appropriate heterozygous PRNP variant and "
        "amplified and analysed as described above. All wells containing mutant gDNA were prepared "
        "in a separate workspace to avoid contamination. Droplets were classified as REF-only, "
        "MUT-only, double-positive or double-negative using fluorescence amplitude thresholds set "
        "from NTC, WT and PRNP-mutant controls. These partition counts, rather than QuantaSoft "
        "copy-per-microlitre estimates, were used for the main fractional-abundance analysis."
    ),
    "Data processing: Raw droplet fluorescence data were exported": (
        "Data processing: QuantaSoft .ddpcr project files were batch-imported in R (v.4.3.2) using "
        "a pipeline adapted from digital MIQE guidelines [30] and the tidyverse and binom packages. "
        "For each active well, the pipeline extracted embedded QuantaSoft JSON metadata, "
        "reconstructed accepted-droplet channel classes, annotated sample identifiers and targeted "
        "variant (E200K, D178N or P102L), and cross-checked the reconstructed counts against the "
        "corresponding CSV exports. Fractional abundance (FA) was calculated on the "
        "Poisson-corrected molecule-occupancy scale, not as the raw fraction of mutant-positive "
        "droplets. For target t, occupancy was estimated from the target-negative accepted-droplet "
        "fraction, and FA was expressed as a percentage:"
    ),
    "Limit of blank: To rule out false positives": (
        "Limit of blank: To rule out false positives from background noise, the limit of blank "
        "(LoB) was computed on the mutant-droplet count scale. Droplets were gated manually using "
        "NTC, WT and mutant controls. WT control wells with at least 10,000 accepted droplets were "
        "used as the biological blanks for each assay and plate, because they contain genomic DNA "
        "background without the targeted mutation. NTCs were retained as contamination controls but "
        "were not pooled into the LoB background model. For each assay and plate, the blank "
        "false-positive rate p0 was estimated as the upper 95% exact binomial confidence bound for "
        "mutant-positive droplets among accepted WT-control droplets. If a sample was assembled "
        "from multiple plates, the highest p0 among contributing plates was used; if a plate lacked "
        "WT blanks, an assay-wide WT-control p0 was used as fallback. The LoB count was defined as "
        "the 95th percentile of the binomial blank distribution:"
    ),
    "Limit of detection: The analytical sensitivity": (
        "Limit of detection: The analytical sensitivity of the SNV assays was assessed "
        "experimentally using artificial mosaic samples for E200K, D178N, and P102L, generated by "
        "mixing gDNA from heterozygous mutation carriers with commercially available wild-type "
        "gDNA (Promega Human Genomic DNA, see above). Mixtures were prepared at predefined final "
        "variant allele frequencies of 1%, 0.33%, 0.1%, and 0.05% and analysed in parallel with "
        "controls. The LoD was defined as the lowest spike-in level whose 95% CI did not overlap "
        "the WT distribution, and LoD thresholds were applied on the percentage FA scale."
    ),
    "Participant-level pooling and quantification": (
        "Participant-level pooling and quantification: For summary analyses at the participant "
        "level, mutant-positive, reference-positive and accepted droplets were aggregated across "
        "all available tissue regions per participant and assay. Participant-level FA and 95% CIs "
        "were then recalculated from the pooled droplet counts using the same Poisson-corrected "
        "target-occupancy method. The incidentally discovered E200K-heterozygous germline carrier "
        "(CJD30) was excluded from E200K mosaicism analyses. For each participant-assay "
        "combination, the maximum p0 across contributing plates was used, and the participant-level "
        "LoB was defined as the 95th binomial quantile with N equal to that participant's total "
        "accepted droplets."
    ),
    "To express assay depth in genome-equivalent terms": (
        "To express assay depth in genome-equivalent terms, target occupancy was estimated from "
        "negative-droplet fractions using a Poisson droplet-occupancy model. The estimated total "
        "number of haploid genome equivalents surveyed was:"
    ),
    "where N denotes the accepted-droplet count": (
        "where N denotes the accepted-droplet count and N0,0 denotes the number of accepted "
        "droplets with no REF or MUT signal. Mutant and reference haploid genome equivalents were "
        "estimated analogously from the corresponding target-negative droplet counts. For pooled "
        "participant-level summaries, droplet counts were aggregated by participant before applying "
        "the same Poisson transformation to the aggregate negative-droplet fraction. Exact binomial "
        "CIs were calculated for each negative-droplet fraction and converted to genome-equivalent "
        "limits using the same transformation."
    ),
    "Definition of mutation-positive samples": (
        "Definition of mutation-positive samples: A sample or participant was considered "
        "mutation-positive only if it passed both decision layers: the observed mutant-positive "
        "droplet count had to exceed the count-scale LoB, and the lower bound of the 95% FA "
        "confidence interval had to exceed the assay-specific LoD."
    ),
    "ddPCR assays targeting D178N": (
        "ddPCR assays targeting D178N, E200K, and P102L showed clear separation between "
        "heterozygous positive controls and negative-control clusters. Low-level mutant-positive "
        "background droplets were present, as expected for rare-mutation ddPCR, and were handled "
        "through the empirical LoB procedure. A custom A117V assay showed poor amplification "
        "efficiency and inadequate cluster separation, likely due to local sequence complexity or "
        "secondary DNA structures, and was not used for analysis; such issues are a common "
        "limitation of ddPCR [47]. This limitation is specific to ddPCR and does not affect the "
        "sequencing-based analyses."
    ),
    "No confirmed tissue-level PRNP mosaicism by multi-region ddPCR": (
        "Technically pooled low-level E200K signals by multi-region ddPCR"
    ),
    "We analysed bulk DNA from six brain regions": (
        "We analysed bulk DNA from six brain regions in 31 CJD cases, comprising 30 sCJD cases and "
        "one case subsequently classified as E200K genetic CJD, together with 8 non-prion controls. "
        "Allele-specific ddPCR assays targeted the pathogenic PRNP variants D178N, E200K and P102L. "
        "A summary of results is provided in Table 3 and Figure 6, with detailed fractional "
        "abundance and technical data in Tables S1 and S2. Based on Poisson occupancy calculations, "
        f"the retained ddPCR data spanned {DDPCR_SUMMARY['accepted_min']} to "
        f"{DDPCR_SUMMARY['accepted_max']} accepted droplets and {DDPCR_SUMMARY['haploid_min']} to "
        f"{DDPCR_SUMMARY['haploid_max']} total haploid genome equivalents per sample-region assay. "
        "The incidentally identified "
        "E200K germline carrier CJD30 showed the expected approximately heterozygous E200K signal "
        f"across all six regions (FA approximately {DDPCR_SUMMARY['cjd30_fa_min']}-"
        f"{DDPCR_SUMMARY['cjd30_fa_max']}%) and was excluded from E200K somatic "
        "mosaicism analyses."
    ),
    "No sample met the prespecified positivity criterion": (
        "After excluding CJD30 from E200K somatic analyses, no non-pooled non-CJD30 sample-region "
        "measurement passed both ddPCR decision criteria. Two technically pooled sample-region "
        f"E200K summaries passed both criteria: {DDPCR_SUMMARY['technical_pooled_pass_description']}. "
        "No D178N or P102L sample-region summary passed both LoB and LoD, and no D178N result met "
        f"criteria for a positive call in the {DDPCR_SUMMARY['mm2t_count']} MM2T/thalamic CJD cases. "
        f"A further {DDPCR_SUMMARY['lob_only_count']} sample-region measurements exceeded the "
        "count-scale LoB but not the FA-scale LoD."
    ),
    "A post hoc subgroup analysis in which droplet counts": (
        f"In participant-level analyses, {DDPCR_SUMMARY['participant_pass_description']} remained "
        "positive after pooling all available regions. The CJD21 technically pooled thalamic E200K "
        "signal did not remain above both thresholds after participant-level pooling "
        f"({DDPCR_SUMMARY['pooled_follow_up_description']}). No participant was positive for D178N "
        "or P102L (Figure S2). Because regional pooling can dilute a focal mosaic signal, these "
        "participant-level analyses were used as a summary view rather than as the sole criterion "
        "for excluding regional mosaicism."
    ),
    "In this study, we tested whether pathogenic somatic PRNP mutations": (
        "In this study, we tested whether pathogenic somatic PRNP mutations contribute to sporadic "
        "Creutzfeldt-Jakob disease (sCJD), whose cause is currently unknown. We combined "
        "decontamination-compatible DNA extraction from prion-infected brain tissue, "
        "high-sensitivity ddPCR across multiple brain regions, and targeted sequencing of the "
        "entire PRNP locus. The corrected ddPCR analysis identified two technically pooled "
        "low-level E200K sample-region summaries, including one participant-level signal in a "
        "sCJD case, but no non-pooled non-CJD30 sample-region measurement passed both thresholds. "
        "The ddPCR data therefore did not identify recurrent or multi-variant evidence of "
        "tissue-level mosaicism for canonical pathogenic PRNP missense variants. Targeted "
        "sequencing identified no somatic coding PRNP variants, although rare low-VAF intronic "
        "candidate variants were observed."
    ),
    "To search for pathogenic PRNP variants throught": (
        "To search for pathogenic PRNP variants across multiple brain regions, we applied highly "
        "sensitive ddPCR assays with empirical limits of detection of 0.056-0.13% VAF. Digital "
        "partitioning allows tens to hundreds of thousands of haploid genome equivalents to be "
        "interrogated per sample, with up to single-allele resolution for the detection of mutant "
        "alleles - one to two orders of magnitude more alleles than are typically interrogated by "
        "current single-cell DNA sequencing experiments. A focal mosaic clone should, therefore, "
        "produce countable mutant alleles against the wild-type background, even in bulk tissue. "
        "The new analysis did detect countable low-level E200K signal in two technically pooled "
        "sample-region summaries, but the pattern was sparse, restricted to E200K, and absent from "
        "non-pooled non-CJD30 sample-region calls."
    ),
    "In some samples, mean fractional-abundance estimates": (
        "The low-level E200K threshold crossings require cautious interpretation. Both arose from "
        "technically pooled sample-region summaries; no corresponding non-pooled non-CJD30 "
        "sample-region measurement passed both thresholds, and only CJD4 persisted after "
        "participant-level pooling. We therefore report these as ddPCR-positive low-level E200K "
        "findings at the technically pooled summary level, not as confirmed somatic drivers of "
        "sCJD. At the low mutant copy numbers relevant to rare mosaicism, stochastic sampling and "
        "Poisson uncertainty can produce elevations that are technically real but biologically "
        "difficult to interpret [71]."
    ),
    "We aimed to profile multiple histological subtypes": (
        "We aimed to profile multiple histological subtypes of sCJD, including the exceedingly "
        "rare thalamic variant, MM2-T, also known as sporadic fatal insomnia (sFI), the most "
        "uncommon form of sporadic prion disease. In the largest published series of this subtype, "
        "only 31 proven cases were identified at major prion surveillance centres worldwide over "
        "several decades [73]. Given its clinicopathological similarity to the genetic prion "
        "disease fatal familial insomnia (FFI), which is caused by the D178N mutation in the "
        "appropriate codon 129 context, it is tempting to hypothesise that thalamic CJD could "
        "arise from somatic D178N mutations. However, even after corrected ddPCR reanalysis, we "
        "found no evidence for D178N mosaicism in any of the three MM2T cases whose samples we "
        "could obtain. This argues against detectable D178N mosaicism in the sampled regions as a "
        "common explanation for thalamic CJD, but does not preclude the possibility that a "
        "mutational focus was missed or that somatic mutations play a role in other patients."
    ),
    "Overall, this study provides no support": (
        "Overall, this study provides no support for common, recurrent tissue-level mosaicism of "
        "canonical pathogenic PRNP missense variants as a driver of sCJD in the samples studied. "
        "Our data are in line with, and extend, prior studies of somatic PRNP variation in sCJD, "
        "which sequenced the protein coding region with limited anatomical sampling. Murley et al. "
        "used high-depth amplicon sequencing of the PRNP coding region in frontal cortex and "
        "reported two candidate somatic missense variants at low VAF, while finding no overall "
        "excess of somatic PRNP variation in sCJD compared with controls [18]. Our ddPCR screen "
        "addressed a different question: whether canonical high-penetrance PRNP mutations recur as "
        "tissue-level mosaics across multiple brain regions. The technically pooled E200K findings "
        "reported here keep rare regional mosaicism on the table, but they do not support a common "
        "or disease-specific pattern. McDonough et al. sequenced the PRNP coding region, principally "
        "in cerebellar cortex, with additional regions being sampled in cases of the visual-onset "
        "Heidenhain variant, and did not identify bona fide somatic PRNP variants [19]. Taken "
        "together, the available data argue against common, reproducible tissue-level PRNP coding "
        "mosaicism in sCJD, while leaving unresolved the significance of technically pooled "
        "low-level E200K ddPCR signals, the rare intronic variants identified in our study, and the low-level ORR "
        "candidate calls reported by Murley et al."
    ),
    "Several limitations should be considered": (
        "Several limitations should be considered. First, analyses were performed on bulk tissue, "
        "so mosaicism restricted to rare cell clones may have been diluted below detection. "
        "Second, anatomical sampling was necessarily incomplete; even six regions represent only a "
        "small fraction of the brain, and the initiating site may vary between cases. Third, ddPCR "
        "profiling was limited to D178N, E200K and P102L because a robust A117V ddPCR assay could "
        "not be established. However, A117V-associated Gerstmann-Straussler-Scheinker (GSS) "
        "disease is clinically and histologically distinct from typical sCJD phenotypes [74,75], "
        "making A117V an unlikely explanation for typical sCJD. Fourth, the technically pooled "
        "low-level E200K ddPCR signals were not independently confirmed in the same anatomical "
        "regions by an orthogonal assay. Larger structural variants and PRNP "
        "copy-number variation were also not assessed. Finally, cohort size was constrained by the "
        "rarity of autopsy tissue from sCJD and by the resource intensity of multi-region analysis "
        "of infectious material."
    ),
    "Taken together, our findings argue against widespread": (
        "Taken together, our findings argue against widespread or recurrent tissue-level mosaicism "
        "for canonical pathogenic PRNP missense variants as a common initiating event in sCJD. "
        "They do not, however, exclude the possibility that somatic mutations in PRNP, or in other "
        "genes relevant to prion biology, act as causative factors in rare tissue foci and "
        "isolated cases. The technically pooled low-level E200K ddPCR findings sharpen this point: they are not "
        "sufficient to establish a common somatic PRNP mechanism, but they show why rare regional "
        "events remain difficult to rule out conclusively. Testing these possibilities will "
        "require larger cohorts, denser anatomical sampling, broader genomic approaches, such as "
        "targeted multi-gene panels or whole-genome sequencing, and eventually single-cell or "
        "ultra-high-depth analysis of large cell populations. The DNA extraction and "
        "decontamination protocol described here establishes a practical basis for further "
        "genetic analysis of prion-infected material. Overall, this study helps delineate the "
        "limits of the PRNP mosaicism model while highlighting unresolved questions that may guide "
        "future studies of sCJD pathogenesis."
    ),
    "Sequencing FASTQ files are available freely": (
        "Sequencing FASTQ files are available freely at xxxx (probably NCBI SRA). Raw ddPCR "
        "project-derived partition counts and processed ddPCR result tables are available freely "
        "at xxx (probably Zenodo). All other non-identifying data are available for academic "
        "purposes upon reasonable request to the corresponding author."
    ),
}


def replace_paragraph_by_prefix(paragraphs: list[ET.Element], prefix: str, new_text: str) -> ET.Element:
    matches = [p for p in paragraphs if paragraph_text(p).startswith(prefix)]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one paragraph starting {prefix!r}, found {len(matches)}")
    p = matches[0]
    ppr = copy.deepcopy(p.find(qn(W_NS, "pPr")))
    p.clear()
    if ppr is not None:
        p.append(ppr)
    p.append(text_run(new_text))
    return p


def replace_paragraph_with_formula(p: ET.Element, text: str) -> None:
    ppr = copy.deepcopy(p.find(qn(W_NS, "pPr")))
    p.clear()
    if ppr is not None:
        p.append(ppr)
    p.append(text_run(text))


def insert_after(body: ET.Element, anchor: ET.Element, new_p: ET.Element) -> ET.Element:
    children = list(body)
    idx = children.index(anchor)
    body.insert(idx + 1, new_p)
    return new_p


def count_words(text: str) -> int:
    return len(re.findall(r"\b\S+\b", text))


def all_paragraphs(root: ET.Element) -> list[ET.Element]:
    return root.findall(f".//{qn(W_NS, 'p')}")


def write_docx_from_tree(source: Path, dest: Path, root: ET.Element) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        with zipfile.ZipFile(source, "r") as zin:
            zin.extractall(tmpdir)
        document_path = tmpdir / "word" / "document.xml"
        ET.ElementTree(root).write(document_path, encoding="utf-8", xml_declaration=True)
        if dest.exists():
            dest.unlink()
        with zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED) as zout:
            for path in tmpdir.rglob("*"):
                if path.is_file():
                    zout.write(path, path.relative_to(tmpdir).as_posix())


def main() -> None:
    shutil.copy2(ORIGINAL_DOCX, CLEAN_DOCX)
    with zipfile.ZipFile(CLEAN_DOCX, "r") as z:
        root = ET.fromstring(z.read("word/document.xml"))

    body = root.find(qn(W_NS, "body"))
    if body is None:
        raise RuntimeError("Could not find Word document body")

    paragraphs = all_paragraphs(root)
    anchors: dict[str, ET.Element] = {}
    for prefix, new_text in REPLACEMENTS.items():
        anchors[prefix] = replace_paragraph_by_prefix(paragraphs, prefix, new_text)
        paragraphs = all_paragraphs(root)

    data_processing = anchors["Data processing: Raw droplet fluorescence data were exported"]
    lambda_p = insert_after(body, data_processing, formula_paragraph(FORMULAS["lambda"]))
    insert_after(body, lambda_p, formula_paragraph(FORMULAS["fa"]))

    lob = anchors["Limit of blank: To rule out false positives"]
    insert_after(body, lob, formula_paragraph(FORMULAS["lob"]))

    formula_intro = anchors["To express assay depth in genome-equivalent terms"]
    gtotal_old = [p for p in all_paragraphs(root) if paragraph_text(p).startswith("Gtotal =")]
    if len(gtotal_old) != 1:
        raise RuntimeError(f"Expected one Gtotal paragraph, found {len(gtotal_old)}")
    replace_paragraph_with_formula(gtotal_old[0], FORMULAS["gtotal"])

    texts = [paragraph_text(p) for p in all_paragraphs(root) if paragraph_text(p)]
    try:
        abstract_start = texts.index("Abstract")
        intro_start = texts.index("Introduction")
        acknowledgements_start = texts.index("Acknowledgements")
    except ValueError as exc:
        raise RuntimeError("Could not locate section boundaries for word counts") from exc

    abstract_count = count_words(" ".join(texts[abstract_start + 1:intro_start]))
    main_count = count_words(" ".join(texts[intro_start + 1:acknowledgements_start]))
    replace_paragraph_by_prefix(all_paragraphs(root), "Abstract word count:", f"Abstract word count: {abstract_count}")
    replace_paragraph_by_prefix(all_paragraphs(root), "Main text word count:", f"Main text word count: {main_count}")

    write_docx_from_tree(CLEAN_DOCX, CLEAN_DOCX, root)
    print(CLEAN_DOCX.as_posix())
    print(f"Abstract word count: {abstract_count}")
    print(f"Main text word count: {main_count}")
    print("Formula anchor:", paragraph_text(formula_intro))


if __name__ == "__main__":
    main()
