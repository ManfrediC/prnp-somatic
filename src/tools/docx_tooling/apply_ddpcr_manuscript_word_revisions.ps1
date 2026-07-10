Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$projectRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..\..")
$manifestPath = Join-Path $projectRoot "manuscript\v2\ddpcr_manuscript_replacements.json"
$originalPath = Join-Path $projectRoot "manuscript\v2\Carta_PRNP_v2.user_reviewed_original_2026-05-11.docx"
$trackedPath = Join-Path $projectRoot "manuscript\v2\Carta_PRNP_v2.codex_tracked.docx"
$auditPath = Join-Path $projectRoot "manuscript\v2\Carta_PRNP_v2.codex_tracked.audit.json"
$tempStem = "Carta_PRNP_v2_codex_" + [System.Guid]::NewGuid().ToString("N")
$cleanPath = Join-Path ([System.IO.Path]::GetTempPath()) ($tempStem + "_clean.docx")
$trackedTempPath = Join-Path ([System.IO.Path]::GetTempPath()) ($tempStem + "_tracked.docx")
$logPath = Join-Path $projectRoot "deposition\review\word_revision_log.txt"

New-Item -ItemType Directory -Force -Path (Split-Path -Parent $logPath) | Out-Null

if (-not (Test-Path -LiteralPath $manifestPath)) {
    throw "Manifest not found: $manifestPath"
}

$manifest = Get-Content -LiteralPath $manifestPath -Raw -Encoding UTF8 | ConvertFrom-Json

function Write-Log {
    param([string]$Message)
    $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    Add-Content -LiteralPath $logPath -Value "$stamp $Message" -Encoding UTF8
}

if (Test-Path -LiteralPath $logPath) {
    Remove-Item -LiteralPath $logPath -Force
}
Write-Log "Starting ddPCR manuscript Word revision script"

if (Test-Path -LiteralPath $cleanPath) {
    Remove-Item -LiteralPath $cleanPath -Force
}
Copy-Item -LiteralPath $originalPath -Destination $cleanPath

$wdCollapseEnd = 0
$wdAlignParagraphCenter = 1
$wdCompareTargetNew = 2
$wdGranularityWordLevel = 1
$wdDoNotSaveChanges = 0
$commentThreadsToRestore = @(
    [ordered]@{
        comment_prefix = "How? Illumina? Nanopore?"
        paragraph_prefix = "Prior studies have begun to examine somatic PRNP variation in sCJD."
        anchor_text = "ultra-deep short-read approach based on multiple independent PCR reactions"
    },
    [ordered]@{
        comment_prefix = "I think that this kind of statement is too sweeping"
        paragraph_prefix = "DNA purified from CJD brain tissue was assessed for prion seeding activity"
        anchor_text = "no detectable seeding activity (Figure 4)"
    }
)

function Get-ParagraphText {
    param($Paragraph)
    return (($Paragraph.Range.Text -replace "[`r`a]+$", "").Trim())
}

function Find-ParagraphByPrefix {
    param(
        $Document,
        [string]$Prefix
    )

    $matches = @()
    foreach ($paragraph in @($Document.Paragraphs)) {
        $text = Get-ParagraphText $paragraph
        if ($text.StartsWith($Prefix, [System.StringComparison]::Ordinal)) {
            $matches += $paragraph
        }
    }

    if ($matches.Count -ne 1) {
        throw "Expected one paragraph starting '$Prefix', found $($matches.Count)"
    }
    return $matches[0]
}

function Replace-ParagraphByPrefix {
    param(
        $Document,
        [string]$Prefix,
        [string]$Text
    )

    $paragraph = Find-ParagraphByPrefix -Document $Document -Prefix $Prefix
    $range = $paragraph.Range
    $range.End = $range.End - 1
    $range.Text = $Text
}

function Delete-ParagraphContentsByPrefix {
    param(
        $Document,
        [string]$Prefix
    )

    $paragraph = Find-ParagraphByPrefix -Document $Document -Prefix $Prefix
    $range = $paragraph.Range
    $range.End = $range.End - 1
    $range.Text = ""
}

function Insert-FormulaAfterPrefix {
    param(
        $Document,
        [string]$Prefix,
        [string]$Text
    )

    $paragraph = Find-ParagraphByPrefix -Document $Document -Prefix $Prefix
    $range = $paragraph.Range.Duplicate
    $range.Collapse($wdCollapseEnd)
    $range.InsertAfter("$Text`r")
}

function Format-FormulaParagraph {
    param(
        $Document,
        [string]$Prefix
    )

    $paragraph = Find-ParagraphByPrefix -Document $Document -Prefix $Prefix
    $paragraph.Range.ParagraphFormat.Alignment = $wdAlignParagraphCenter
    $paragraph.Range.Font.Italic = $true
    return $true
}

function Get-CommentText {
    param($Comment)
    return [string]$Comment.Range.Text
}

function Find-CommentByTextPrefix {
    param(
        $Document,
        [string]$Prefix
    )

    $matches = @()
    foreach ($comment in @($Document.Comments)) {
        if ((Get-CommentText $comment).StartsWith($Prefix, [System.StringComparison]::Ordinal)) {
            $matches += $comment
        }
    }
    if ($matches.Count -ne 1) {
        throw "Expected one comment starting '$Prefix', found $($matches.Count)"
    }
    return $matches[0]
}

function Test-CommentTextPresent {
    param(
        $Document,
        [string]$Text
    )

    foreach ($comment in @($Document.Comments)) {
        if ((Get-CommentText $comment) -eq $Text) {
            return $true
        }
    }
    return $false
}

function Find-TextRangeInParagraph {
    param(
        $Document,
        [string]$ParagraphPrefix,
        [string]$Text
    )

    $paragraph = Find-ParagraphByPrefix -Document $Document -Prefix $ParagraphPrefix
    $range = $paragraph.Range.Duplicate
    $range.Find.ClearFormatting()
    $range.Find.Text = $Text
    if (-not $range.Find.Execute()) {
        throw "Could not find comment anchor '$Text' in paragraph '$ParagraphPrefix'"
    }
    return $range
}

function Restore-CommentThread {
    param(
        $OriginalDocument,
        $ComparisonDocument,
        [string]$CommentPrefix,
        [string]$ParagraphPrefix,
        [string]$AnchorText
    )

    $sourceComment = Find-CommentByTextPrefix -Document $OriginalDocument -Prefix $CommentPrefix
    $sourceText = Get-CommentText $sourceComment
    if (Test-CommentTextPresent -Document $ComparisonDocument -Text $sourceText) {
        return 0
    }

    $anchorRange = Find-TextRangeInParagraph -Document $ComparisonDocument -ParagraphPrefix $ParagraphPrefix -Text $AnchorText
    $restoredComment = $ComparisonDocument.Comments.Add($anchorRange, $sourceText)
    $restoredComment.Author = $sourceComment.Author
    $restoredComment.Initial = $sourceComment.Initial
    $restoredCount = 1

    foreach ($sourceReply in @($sourceComment.Replies)) {
        $replyText = Get-CommentText $sourceReply
        if (-not (Test-CommentTextPresent -Document $ComparisonDocument -Text $replyText)) {
            $restoredReply = $restoredComment.Replies.Add($anchorRange, $replyText)
            $restoredReply.Author = $sourceReply.Author
            $restoredReply.Initial = $sourceReply.Initial
            $restoredCount += 1
        }
    }
    return $restoredCount
}

$word = $null
$documentsToClose = @()
$formattedFormulaCount = 0
$restoredCommentCount = 0
$baseCommentTexts = @()

try {
    Write-Log "Opening Microsoft Word"
    $word = New-Object -ComObject Word.Application
    $word.Visible = $false
    $word.DisplayAlerts = 0

    Write-Log "Opening clean revised copy"
    $doc = $word.Documents.Open($cleanPath)
    $documentsToClose += $doc

    foreach ($replacement in $manifest.replacements) {
        Write-Log "Replacing paragraph prefix: $($replacement.prefix.Substring(0, [Math]::Min(60, $replacement.prefix.Length)))"
        Replace-ParagraphByPrefix -Document $doc -Prefix $replacement.prefix -Text $replacement.text
    }

    foreach ($formula in $manifest.formula_insertions) {
        Write-Log "Inserting formula after prefix: $($formula.after_prefix.Substring(0, [Math]::Min(60, $formula.after_prefix.Length)))"
        Insert-FormulaAfterPrefix -Document $doc -Prefix $formula.after_prefix -Text $formula.text
    }

    foreach ($formula in $manifest.formula_replacements) {
        Write-Log "Replacing formula prefix: $($formula.prefix)"
        Replace-ParagraphByPrefix -Document $doc -Prefix $formula.prefix -Text $formula.text
    }

    foreach ($replacement in $manifest.reviewer_replacements) {
        Write-Log "Applying reviewer replacement prefix: $($replacement.prefix.Substring(0, [Math]::Min(60, $replacement.prefix.Length)))"
        Replace-ParagraphByPrefix -Document $doc -Prefix $replacement.prefix -Text $replacement.text
    }

    foreach ($deletion in $manifest.reviewer_deletions) {
        Write-Log "Applying reviewer deletion prefix: $($deletion.prefix.Substring(0, [Math]::Min(60, $deletion.prefix.Length)))"
        Delete-ParagraphContentsByPrefix -Document $doc -Prefix $deletion.prefix
    }

    Write-Log "Updating manuscript word counts"
    Replace-ParagraphByPrefix -Document $doc -Prefix "Abstract word count:" -Text "Abstract word count: $($manifest.word_counts.abstract)"
    Replace-ParagraphByPrefix -Document $doc -Prefix "Main text word count:" -Text "Main text word count: $($manifest.word_counts.main)"

    foreach ($formula in $manifest.formula_insertions) {
        if (Format-FormulaParagraph -Document $doc -Prefix $formula.text) {
            $formattedFormulaCount += 1
        }
    }
    foreach ($formula in $manifest.formula_replacements) {
        if (Format-FormulaParagraph -Document $doc -Prefix $formula.text) {
            $formattedFormulaCount += 1
        }
    }

    Write-Log "Saving clean revised copy"
    $doc.Save()
    $doc.Close($wdDoNotSaveChanges)
    $documentsToClose = @()

    Write-Log "Opening original and clean revised copy for comparison"
    $originalDoc = $word.Documents.Open($originalPath, $false, $true)
    $revisedDoc = $word.Documents.Open($cleanPath, $false, $true)
    $documentsToClose += $originalDoc
    $documentsToClose += $revisedDoc
    $baseCommentTexts = @($originalDoc.Comments | ForEach-Object { Get-CommentText $_ })

    Write-Log "Creating tracked comparison"
    $comparisonDoc = $word.CompareDocuments(
        $originalDoc,
        $revisedDoc,
        $wdCompareTargetNew,
        $wdGranularityWordLevel,
        $false,
        $true,
        $true,
        $true,
        $true,
        $true,
        $true,
        $true,
        $true,
        $true,
        "Codex ddPCR",
        $true
    )
    $documentsToClose += $comparisonDoc
    Write-Log "Restoring comment threads whose anchors were replaced"
    foreach ($commentThread in $commentThreadsToRestore) {
        $restoredCommentCount += Restore-CommentThread `
            -OriginalDocument $originalDoc `
            -ComparisonDocument $comparisonDoc `
            -CommentPrefix $commentThread.comment_prefix `
            -ParagraphPrefix $commentThread.paragraph_prefix `
            -AnchorText $commentThread.anchor_text
    }
    $comparisonCommentTexts = @($comparisonDoc.Comments | ForEach-Object { Get-CommentText $_ })
    $missingCommentTexts = @($baseCommentTexts | Where-Object { $_ -notin $comparisonCommentTexts })
    if ($missingCommentTexts.Count -ne 0 -or $comparisonDoc.Comments.Count -ne $baseCommentTexts.Count) {
        throw "Comment preservation failed: expected $($baseCommentTexts.Count), found $($comparisonDoc.Comments.Count), missing $($missingCommentTexts.Count)."
    }
    Write-Log "Saving tracked comparison to temporary path"
    $comparisonDoc.SaveAs2([string]$trackedTempPath)

    $comparisonDoc.Close($wdDoNotSaveChanges)
    $revisedDoc.Close($wdDoNotSaveChanges)
    $originalDoc.Close($wdDoNotSaveChanges)
    $documentsToClose = @()

    Write-Log "Reopening tracked comparison for Word validation"
    $validationDoc = $word.Documents.Open($trackedTempPath, $false, $true)
    $documentsToClose += $validationDoc
    $validatedRevisionCount = $validationDoc.Revisions.Count
    $validatedCommentCount = $validationDoc.Comments.Count
    $validatedCommentTexts = @($validationDoc.Comments | ForEach-Object { Get-CommentText $_ })
    $validatedText = $validationDoc.Content.Text
    if ($validatedRevisionCount -lt 1) {
        throw "Word validation found no tracked revisions."
    }
    if (-not $validatedText.Contains("SRP714079") -or -not $validatedText.Contains("10.5281/zenodo.21227118")) {
        throw "Word validation could not find the final SRA study and Zenodo DOI text."
    }
    $missingValidatedCommentTexts = @($baseCommentTexts | Where-Object { $_ -notin $validatedCommentTexts })
    if ($validatedCommentCount -ne $baseCommentTexts.Count -or $missingValidatedCommentTexts.Count -ne 0) {
        throw "Word validation did not preserve all original comments."
    }
    $validationDoc.Close($wdDoNotSaveChanges)
    $documentsToClose = @()

    Move-Item -LiteralPath $trackedTempPath -Destination $trackedPath -Force

    $audit = [ordered]@{
        base = (Resolve-Path $originalPath).Path
        base_sha256 = (Get-FileHash $originalPath -Algorithm SHA256).Hash.ToLowerInvariant()
        output = (Resolve-Path $trackedPath).Path
        output_sha256 = (Get-FileHash $trackedPath -Algorithm SHA256).Hash.ToLowerInvariant()
        author = "Codex"
        replacement_count = @($manifest.replacements).Count
        formula_insertion_count = @($manifest.formula_insertions).Count
        formula_replacement_count = @($manifest.formula_replacements).Count
        reviewer_replacement_count = @($manifest.reviewer_replacements).Count
        reviewer_deletion_count = @($manifest.reviewer_deletions).Count
        word_revision_count = $validatedRevisionCount
        original_comment_count = $baseCommentTexts.Count
        restored_comment_count = $restoredCommentCount
        output_comment_count = $validatedCommentCount
        original_comment_texts_preserved = $true
        word_reopen_validation = "pass"
        final_sra_study_present = $true
        final_zenodo_doi_present = $true
        generated_utc = (Get-Date).ToUniversalTime().ToString("yyyy-MM-ddTHH:mm:ssZ")
    }
    $audit | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $auditPath -Encoding UTF8

    Write-Host "Tracked-edit DOCX: $trackedPath"
    Write-Host "Word revision count: $validatedRevisionCount"
    Write-Host "Reviewer comments preserved: $validatedCommentCount/$($baseCommentTexts.Count)"
    Write-Host "Formula paragraphs formatted as centred display lines: $formattedFormulaCount"
    Write-Log "Completed successfully"
}
finally {
    foreach ($document in $documentsToClose) {
        try {
            $document.Close($wdDoNotSaveChanges)
        }
        catch {
        }
    }
    if ($null -ne $word) {
        $word.Quit()
    }
    foreach ($temporaryPath in @($cleanPath, $trackedTempPath)) {
        if (Test-Path -LiteralPath $temporaryPath) {
            Remove-Item -LiteralPath $temporaryPath -Force
        }
    }
}
