Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$projectRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..\..")
$manifestPath = Join-Path $projectRoot "manuscript\v2\ddpcr_manuscript_replacements.json"
$originalPath = Join-Path $projectRoot "manuscript\v2\2026_05_11_Carta_PRNP_v2.docx"
$cleanPath = Join-Path $projectRoot "manuscript\v2\2026_05_11_Carta_PRNP_v2.ddpcr_corrected_clean_word.docx"
$trackedPath = Join-Path $projectRoot "manuscript\v2\2026_05_11_Carta_PRNP_v2.ddpcr_corrected_tracked.docx"
$logPath = Join-Path $projectRoot "results\ddPCR\validation\word_revision_log.txt"

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
if (Test-Path -LiteralPath $trackedPath) {
    Remove-Item -LiteralPath $trackedPath -Force
}
Copy-Item -LiteralPath $originalPath -Destination $cleanPath

$wdCollapseEnd = 0
$wdAlignParagraphCenter = 1
$wdCompareTargetNew = 2
$wdGranularityWordLevel = 1
$wdDoNotSaveChanges = 0

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

$word = $null
$documentsToClose = @()
$formattedFormulaCount = 0

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
        $false,
        $true,
        "Codex ddPCR",
        $true
    )
    $documentsToClose += $comparisonDoc
    Write-Log "Saving tracked comparison"
    $comparisonDoc.SaveAs2([string]$trackedPath)

    $comparisonDoc.Close($wdDoNotSaveChanges)
    $revisedDoc.Close($wdDoNotSaveChanges)
    $originalDoc.Close($wdDoNotSaveChanges)
    $documentsToClose = @()

    Write-Host "Clean revised DOCX: $cleanPath"
    Write-Host "Tracked-edit DOCX: $trackedPath"
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
}
