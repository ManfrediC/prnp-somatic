# Current goal

## Outcome

Create a new copy of `manuscript/v2/Carta_PRNP_v2_final.docx` that preserves
all existing content and formatting, inserts the 21 supplied replacement
passages immediately below their corresponding passages as tracked additions
authored by `Manfredi`, and highlights the duplicated LoD passage in orange.

## Verification surface

- Confirm all 21 replacement passages are present exactly once as `w:ins`
  revisions authored by `Manfredi`.
- Confirm no existing text is deleted or replaced and the source DOCX remains
  byte-for-byte unchanged.
- Confirm each tracked addition directly follows its corresponding source
  passage in document order.
- Confirm the duplicated LoD wording is highlighted in orange.
- Render the output to page PNGs and inspect every page for layout defects.

## Constraints and boundaries

- Do not alter the source DOCX.
- Do not accept or reject any existing revisions or remove existing comments.
- Preserve the manuscript's styles, relationships, figures, tables, fields,
  section structure, and metadata except where revision tracking requires an
  update.
- Use the exact supplied replacement wording, converting Markdown emphasis to
  equivalent Word bold/italic formatting rather than leaving Markdown marks.

## Completion evidence

- Source/output SHA-256 hashes recorded and source hash unchanged.
- Structural audit reports 21 new tracked paragraphs by `Manfredi`, zero new
  deletions, correct adjacency, exact replacement text, and orange duplicate
  highlighting.
- Final DOCX ZIP integrity check passes.
- Latest render contains a complete page set and every page is visually clean.

## Stop conditions

- A supplied passage cannot be matched unambiguously to manuscript content.
- Structural insertion would require removing or overwriting manuscript text.
- The document cannot be saved as a valid DOCX package.

## Verification result - 2026-08-03

- Created `manuscript/v2/Carta_PRNP_v2_final_wording_revisions.docx` without
  modifying the source manuscript.
- Source SHA-256 remained
  `c4e63a1f4e06e20552e0042506ac55df5107a578d2d5fe9d326cd6a6cb98438c`.
- Output SHA-256 is
  `46aec93047e8f8d6123bbd090a0ce80a3adecee9c4e1f8b6ae4704933776b919`.
- Structural audit found 21 tracked content insertions, 21 tracked paragraph
  marks, zero tracked deletions, exact supplied wording, correct adjacency,
  one orange-highlighted duplicated LoD sentence, and changes limited to
  `word/document.xml` and `word/settings.xml`.
- Microsoft Word opened the output read-only and reported revision tracking
  enabled, exactly 21 insertion revisions by `Manfredi`, all 19 existing
  comments preserved, and 33 pages.
- LibreOffice was unavailable. Microsoft Word exported a read-only PDF, all 33
  pages were rendered to PNG and visually inspected, and no layout defects
  were found. The blank final page is byte-identical to the source manuscript's
  blank final page.

## Current second-round goal - 2026-08-03

Create the second-round DOCX at
`manuscript/v2/Carta_PRNP_v2_final_wording_revisions2.docx` from its current
contents, adding the 17 supplied snippets immediately below the 17 current
yellow-highlighted source paragraphs as tracked insertions authored by
`Manfredi`.

Verification surface:

- Existing document text, comments and formatting remain intact.
- Exactly 17 new content insertions and 17 new paragraph-mark insertions are
  authored by `Manfredi`; no tracked deletions are added.
- The 17 inserted texts match the supplied snippets exactly and are adjacent
  to their corresponding yellow paragraphs.
- The output is a valid DOCX and renders without layout defects.

## Completion evidence - second round - 2026-08-03

- Created `manuscript/v2/Carta_PRNP_v2_final_wording_revisions2.docx` from the
  current second-round baseline without deleting or replacing existing text.
- Baseline SHA-256:
  `53d647cd74b093a61d89ff978aa357a2235d5738da1a3269332d4c40a417d787`.
- Output SHA-256:
  `408d4cebae2a5f697b12a2ae20d43f0d2c94516fddcb6add3999ba7b9d0e995f`.
- Structural audit passed: 17 new tracked content insertions and 17 paragraph
  marks authored by `Manfredi`, zero tracked deletions, exact text and correct
  adjacency; 19 existing comments were preserved.
- Microsoft Word reopened the output read-only with tracking enabled and 32
  pages. DOCX ZIP integrity passed.
- LibreOffice was unavailable; Word/PDF fallback rendered all 32 pages, which
  were visually inspected with no layout defects.

## Numbering-only harmonisation - 2026-08-03

Outcome:

- Make `manuscript/v2/Carta_PRNP_v2.docx` and the manuscript figure/table
  wrappers use one numbering scheme: Figures S1-S4 and Tables S1-S5.
- Change only the well-level DOCX references from Figure S3/Table S7 to
  Figure S4/Table S5; preserve the participant-pooled Figure S3 reference.
- Do not alter captions, scientific wording, numerical values, source data, or
  figure/table content.

Verification:

- DOCX context-specific reference audit and ZIP/comment-preservation check. The
  canonical DOCX contained no separate participant-pooled Figure S3 text
  reference; the only Figure S3 reference was the well-level sentence changed
  here, so no other manuscript Figure S3 reference was altered.
- Wrapper/crosswalk order and contiguous-label audit.
- PDF compilation and rendered-page inspection. DOCX rendering was attempted
  with LibreOffice and the Word fallback, but neither headless session produced
  a PDF in this environment.

## Regenerate Table 3 from current ddPCR outputs - 2026-08-03

Outcome:

- Make `src/ddpcr/ddpcr_sample_number.R` the sole Table 3 generator.
- Preserve the existing booktabs layout while replacing stale hard-coded rows
  in the combined tables bundle with generated TeX.
- Leave `manuscript/v2/Carta_PRNP_v2.docx`, ddPCR calculations, numbering and
  all other tables unchanged.

Verification:

- Generated CSV and TeX contain the same 41 rows and the expected CJD4,
  CJD21, CJD30 and aggregate values.
- The combined source contains one `tab:ddPCR-results` label and no hard-coded
  stale Table 3 block.
- The rebuilt Table 3 retains booktabs rules and spacing, fits on one page,
  and passes rendered visual inspection.

Stop conditions:

- Regeneration changes the source workbook or ddPCR calculation logic.
- The generated table cannot reproduce the current layout without modifying
  unrelated tables.

## Participant-pooled Figure S3 tracked reference - 2026-08-03

Outcome:

- Add one suitable participant-pooled Figure S3 reference to the canonical DOCX
  as a tracked insertion, without changing captions, values, or existing
  numbering references.

Verification:

- Confirmed one `<w:ins>` insertion with the new Figure S3 sentence,
  `w:trackRevisions`, preserved comments, no stale Table S7 reference, and
  valid DOCX ZIP integrity.

## Figure S3 participant-pooled legend correction - 2026-08-03

Outcome:

- Replace the obsolete all-pooled-below-LoD interpretation in both Figure S3
  legend sources with the verified CJD4-positive/CJD21-negative participant-
  pooling interpretation.

Verification:

- Source and participant-pooled CSV checks confirmed CJD4 FA 0.103% is above
  both thresholds and CJD21 FA 0.041% is above neither.
- Rebuilt the 12-page figure bundle; stale wording was absent, corrected text
  was present, and all rendered pages were visually inspected.

## Figure 7 sample-region legend correction - 2026-08-03

Outcome:

- Remove the obsolete Figure 7 no-positive-call claim and describe the verified
  pooled CJD4 pons and CJD21 thalamus E200K sample-region results.

Verification:

- Both Figure 7 source captions contain the corrected LoB/LoD interpretation
  and values; the stale claim is absent.
- The rebuilt 12-page figure bundle contains the corrected Figure 7 caption;
  all pages were rendered and visually inspected without layout defects.

## Canonical supplementary Table S5 - 2026-08-03

Outcome:

- Make the current E200K preparation-level R output the sole generator of the
  complete, booktabs-formatted supplementary Table S5.
- Include exactly five individual-preparation rows and two pooled rows for
  CJD4 pons and CJD21 thalamus, excluding germline carrier CJD30.
- Replace the hard-coded bundle copy with the generated TeX while leaving the
  ddPCR calculations, Figure S4, other tables and canonical DOCX unchanged.

Verification:

- Generated CSV/TeX row and value audit against `selected_wells.csv`, the raw
  manifest and `SNV_data_final.xlsx`, including pooled count reconciliation
  and LoB/LoD flag checks.
- One-label/one-input source audit and two-pass PDF compilation.
- Rendered visual inspection of Table S5 for booktabs styling and clipping.
- Before/after SHA-256 equality for `manuscript/v2/Carta_PRNP_v2.docx`.

Stop conditions:

- Any change to ddPCR calculation logic, Figure S4, another table or the DOCX.
- Failure to reproduce the seven specified rows or their pooled totals.

## Manuscript figure and table consistency audit - 2026-08-03

Outcome:

- Audit the current canonical DOCX, figures bundle and tables bundle without
  changing scientific content or manuscript artefacts.

Verification surface:

- Establish that both combined PDFs are freshly compiled from their current
  TeX wrappers and included assets.
- Inspect every rendered PDF page and reconcile captions, displayed values and
  figure/table labels against the bundle source, current source data and DOCX
  references.
- Check DOCX figure/table references for sequential primary and supplementary
  numbering, complete caption matching, and stale or missing references.
- Search figure/table generation code and produced assets for credible unused
  outputs, while distinguishing diagnostic/intermediate files from candidate
  manuscript figures and tables.

Constraints and boundaries:

- Manuscript checks only. Do not edit the DOCX, figures, tables, data,
  captions, or scientific interpretation during this audit.

Completion evidence:

- A page-level inspection record for both PDFs, source/build-log evidence,
  cross-reference audit results, value reconciliations, and a concise list of
  discrepancies or unresolved limitations.

Stop conditions:

- A required current source dataset cannot be identified, or a PDF has no
  traceable compilation/source relationship.

## Canonical manuscript figure/table remediation - 2026-08-03

Outcome:

- Make the canonical DOCX, reusable wrappers, combined bundles and crosswalk
  agree exactly on Figures 1--8/S1--S4 and Tables 1--7/S1--S5.
- Keep each promoted manuscript source and generated manuscript asset in one
  canonical location under `manuscript/`.
- Generate Table 2 and its audit CSV from one LoD dataset with ascending
  confidence-interval bounds.

Verification surface:

- Run every named figure/table generator and require each declared output to
  exist and be non-empty.
- Build both combined bundles twice from the repository root with portable
  MiKTeX; require deterministic second-pass outputs and 12 pages per bundle.
- Render and inspect all 24 PDF pages for numbering, captions, clipping,
  overlap and missing content.
- Run the permanent manuscript consistency checker and checksum verifier.
- Compare visible, non-deleted DOCX references with the wrappers, crosswalk and
  combined PDF labels.

Constraints and boundaries:

- Do not edit the canonical DOCX, manuscript body text, source data, biological
  methods, scientific logic or plotted values.
- Preserve existing user-owned VAF terminology changes in the dirty worktree.
- The only scientific-data presentation change is normalising reversed LoD CI
  endpoints into ascending order.
- Leave all changes uncommitted until the user explicitly authorises a commit.

Stop conditions:

- A named generator changes source data, calculations or plotted values.
- The canonical DOCX or an unrelated tracked file changes.
- Either bundle cannot be reproduced from the canonical wrappers and assets.

## Completion evidence - 2026-08-03

- Canonical wrappers, promoted sources, standalone legends, LoD pipeline,
  crosswalk, documentation, checker and manifest were updated. The DOCX was
  not modified.
- The lollipop pipeline now starts from the recovered R/trackViewer baseline
  SVG, retains its exact head path bodies, applies the measured reviewed edits
  in Python, and writes non-empty canonical SVG/PDF outputs. The output hashes
  for all four head paths match the recovered baseline hashes.
- All named generators completed successfully. Both bundles were rebuilt twice
  from the repository root with portable MiKTeX and explicit manuscript output
  directories; each bundle has 12 pages.
  All 24 final rendered pages were inspected. The consistency checker passed
  and checksum write/check passed.
- No commit was created; the worktree remains pending user review and explicit
  commit permission.

## Current goal - Figure 4 selectable vector text - 2026-08-05

Outcome:

- Rebuild Figure 4 from the five existing raster RT-QuIC panels so that all
  visible axis titles and tick labels are selectable vector text.
- Preserve the coloured traces as lossless raster crops at a measured effective
  resolution of at least 300 ppi in the final figures bundle.
- Preserve the current Figure 4 as an explicit PDF backup and leave every
  original source panel unchanged.

Verification surface:

- Five lossless trace-only crops match exact pixel slices of the embedded source
  images; no OCR, resampling, sharpening or scientific-data reconstruction.
- The generated hybrid panel PDF has five pages, embedded fonts, selectable
  axis text and no Type 3/unembedded fonts introduced by the change.
- `pdfimages -list` reports at least 300 x/y ppi for all five Figure 4 traces
  in the final 12-page figures bundle.
- The Figure 4 caption, subcaptions, labels, layout and scientific traces are
  unchanged; the combined bundle remains 12 pages and passes the manuscript
  consistency checker.
- The canonical DOCX and all original RT-QuIC PDFs retain their baseline
  SHA-256 hashes.

Constraints and boundaries:

- Do not alter the canonical DOCX, source data, plotted values or captions.
- Keep the existing user-owned changes in the reproducibility checksum and
  manifest files; update them only with targeted, reviewed edits.
- Keep the current raster Figure 4 in
  `manuscript/figures/rt-quic/backup/figure_rt_quic_decontamination_raster_original.pdf`.
- Leave all changes uncommitted until the user explicitly authorises a commit.

Stop conditions:

- A crop would remove or obscure a trace, or its data-coordinate mapping is
  uncertain.
- The final raster content cannot meet 300 ppi without unacceptable layout or
  legibility loss.
- Selectable text, embedded-font or full-bundle verification fails.

Completion evidence - 2026-08-05:

- The five crop outputs match the verified source SHA-256/dimensions and retain
  the traces without raster labels or axis spines.
- The five-page hybrid panel PDF and the canonical 12-page figures bundle were
  both compiled twice with portable MiKTeX. Page 4 exposes selectable axis
  labels and uses embedded Type 1 fonts; trace effective resolutions are
  303.1--318.8 ppi.
- The original RT-QuIC PDFs, unused CJD Trizol RNA PDF and canonical DOCX
  retain their baseline SHA-256 hashes. The raster Figure 4 backup remains at
  `manuscript/figures/rt-quic/backup/figure_rt_quic_decontamination_raster_original.pdf`.
- All 12 rendered bundle pages and Figure 4 at 300 dpi were visually checked;
  the manuscript consistency checker and global checksum write/check passed.
- The alignment correction uses one 29.4 x 26.3 mm plot rectangle and identical
  standalone page bounds for all five panels; the subplot upper edges now align
  in the final rendered Figure 4. The minimum final trace resolution is
  303.2 ppi.
- No files were staged or committed.

## Supplementary ddPCR panel labels and typography - 2026-08-05

Outcome:

- Make Supplementary Figure S2 explicitly map panels A--D to D178N,
  E--H to E200K and I--L to P102L, with the four gating stages identified.
- Add visible A--C labels to Supplementary Figure S3 for D178N, E200K and
  P102L.
- Increase plot-source typography so both figures are legible at their final
  manuscript size without changing plotted data, gates or thresholds.

Verification surface:

- Rebuild both figure assets and the complete supplementary-figure bundle from
  their canonical sources.
- Confirm the panel letters and caption mappings agree with the plotted order.
- Inspect final rendered pages for legibility, clipping and overlap; confirm
  fonts are embedded and figure content remains vector-based.
- Run the manuscript consistency checker and verify the ddPCR source workbooks
  and canonical DOCX retain their baseline hashes.

Constraints and boundaries:

- Preserve all current user-owned Figure 4 and reproducibility changes.
- Do not change assay ordering, selected wells, thresholds, VAF values,
  confidence intervals or LoB/LoD classifications.
- Leave all changes uncommitted.

Completion evidence:

- The dedicated S2 generator wrote 12 ordered source plots: four stages each
  for D178N, E200K and P102L. The assembled SVG contains A--L exactly once and
  the caption now states both assay and stage mappings.
- The S3 source writes visible A--C labels for D178N, E200K and P102L and uses
  Cairo PDF output so Arial fonts are embedded and subsetted.
- Both final pages were rendered at 200 dpi and visually inspected. All panel
  letters, titles, ticks, axes and legends are legible, with no clipping or
  overlap.
- Portable MiKTeX produced a 12-page bundle twice. A later canonical-PDF file
  lock prevented another in-place rebuild after generated-preview cleanup, so
  the latest sources were built twice in `C:\tmp`; raster hashes for pages 10
  and 11 matched the canonical bundle exactly.
- The consistency checker and diff checks passed. The ddPCR source workbooks,
  canonical DOCX and unrelated Figure 6 asset retained their baseline hashes.

## Current goal - Figure 8 plot and legend labels - 2026-08-05

Outcome:

- Add a visible `VAF (%)` y-axis label and label the dashed threshold line
  `LoD (0.81%)` in Figure 8.
- State that the horizontal coordinates use GRCh38 and replace `Node colour`
  with `Point colour` in the Figure 8 legend.

Verification surface:

- The canonical SVG contains each requested plot label exactly once as vector
  text, and the generated PDF displays both labels without clipping or overlap.
- The Figure 8 legend contains `GRCh38 genomic coordinates`, `Point colour`,
  and `Vertical axis: VAF (%)`, with the obsolete wording absent.
- The lollipop data CSV and all four plotted head-path hashes remain unchanged.
- The figures bundle remains 12 pages, Figure 8 renders cleanly at final page
  size, and the manuscript consistency and checksum checks pass.

Constraints and boundaries:

- Do not alter the canonical DOCX, source data, variant positions, VAF values,
  colours, point geometry, exon track or any other figure caption.
- Preserve existing user-owned worktree changes and leave this work
  uncommitted unless the user explicitly authorises a commit.

Stop conditions:

- Either requested plot label cannot be positioned legibly without changing
  the scientific plot geometry.
- Regeneration changes the lollipop data or any plotted head path.

Completion evidence:

- The canonical SVG and PDF contain selectable `VAF (%)` and `LoD (0.81%)`
  text exactly once. The y-axis title has dedicated left canvas padding, and
  the LoD label is bold and separated from the dashed line.
- The Figure 8 legend now uses `Point colour`, specifies GRCh38 genomic
  coordinates and describes the vertical axis as `VAF (%)`.
- The lollipop data CSV retained SHA-256
  `dea0cc99e0123cae0151723a4dc27bbda35dfc909f14f6a7b5ba52537d4c807f`,
  the reviewed template retained SHA-256
  `6bc3b7f6b1e7295e037de6105d102252cf1b87406fa57a0410ed2cd25bf80f8a`,
  and all four plotted head-path hashes remained unchanged.
- The 12-page figures bundle compiled twice with portable MiKTeX. Figure 8
  was rendered at final page size and inspected without clipping, overlap or
  caption-layout defects. The manuscript consistency checker passed.
- All five Figure 8 checksum entries match their current files. The global
  checksum verifier remains blocked only by an unrelated pre-existing stale
  checksum for `manuscript/config/figure_table_crosswalk.tsv`.
- No files were staged or committed.

Final refinement status:

- `VAF (%)` now uses explicit 12 px Arial, matching the rendered cap height of
  the `Exon1` key. The right edges of `LoD (0.81%)` and `CJD23` align within
  0.0005 SVG units in the canonical standalone SVG/PDF.
- A twice-compiled 12-page bundle candidate at
  `tmp/pdfs/figure8-alignment-bundle/figures_with_legends.pdf` passed text and
  visual inspection. Adobe Acrobat currently holds the canonical combined PDF
  open, so replacing and checksumming that one file remains pending.
