# Insert a tracked-change sentence into the canonical manuscript DOCX.
# Adds a reference to the new Supplementary Figure S3 and Table S7 in the
# Results paragraph that announces the E200K-positive ddPCR findings.
# Author of the tracked change: Manfredi.

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$projectRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..\..")
$canonicalPath = Join-Path $projectRoot "manuscript\v2\Carta_PRNP_v2_final.docx"
$logPath = Join-Path $projectRoot "deposition\review\word_revision_log.txt"

if (-not (Test-Path -LiteralPath $canonicalPath)) {
    throw "Canonical DOCX not found: $canonicalPath"
}

# The anchor sentence already present in the DOCX (Results paragraph 111).
# The new sentence is inserted immediately AFTER this sentence.
$anchorText = "and the overall distribution of sample-region fractional abundances is shown in Figure 7."
$insertText = " Well-level fractional-abundance estimates and confidence intervals for the two E200K-positive samples are shown in Figure S3, with the underlying per-well and pooled droplet counts in Table S7."

$author = "Manfredi"
$initials = "MC"

$wdCollapseEnd = 0
$wdDoNotSaveChanges = 0
$wdFindStop = 1

$word = $null
$documentsToClose = @()

try {
    Write-Host "Opening Microsoft Word"
    $word = New-Object -ComObject Word.Application
    $word.Visible = $false
    $word.DisplayAlerts = 0

    Write-Host "Opening canonical DOCX: $canonicalPath"
    $doc = $word.Documents.Open($canonicalPath)
    $documentsToClose += $doc

    # Enable track changes and set the author for this session.
    $doc.TrackRevisions = $true
    $word.UserName = $author
    $word.UserInitials = $initials

    # Find the anchor sentence in the document body.
    $range = $doc.Content.Duplicate
    $range.Find.ClearFormatting()
    $range.Find.Text = $anchorText
    $range.Find.Forward = $true
    $range.Find.Wrap = $wdFindStop
    if (-not $range.Find.Execute()) {
        throw "Anchor text not found in document: $anchorText"
    }

    # Collapse to end of the found sentence, then insert the new sentence.
    $range.Collapse($wdCollapseEnd)
    $range.InsertAfter($insertText)

    # Count revisions to confirm the insertion was tracked.
    $revisionCount = $doc.Revisions.Count
    Write-Host "Tracked revisions after insertion: $revisionCount"
    if ($revisionCount -lt 1) {
        throw "No tracked revision was created; the insertion was not recorded."
    }

    # Verify the inserted text is present in the accepted-content view.
    $acceptedText = $doc.Content.Text
    if (-not $acceptedText.Contains("Figure S3") -or -not $acceptedText.Contains("Table S7")) {
        throw "Inserted text (Figure S3 / Table S7) not found in document content after insertion."
    }

    Write-Host "Saving canonical DOCX with tracked change"
    $doc.Save()
    $doc.Close($wdDoNotSaveChanges)
    $documentsToClose = @()

    Write-Host "Done. Tracked change by '$author' inserted successfully."
}
finally {
    foreach ($document in $documentsToClose) {
        try { $document.Close($wdDoNotSaveChanges) } catch { }
    }
    if ($null -ne $word) { $word.Quit() }
}