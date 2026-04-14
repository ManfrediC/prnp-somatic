[CmdletBinding()]
param(
    [string]$SureSelectRoot = 'C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Experiments',
    [string]$DdPcrRoot = 'C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\ddPCR',
    [string]$SamplesRoot = 'C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Samples',
    [string]$OutputRun = '',
    [string]$RepoRoot = '',
    [switch]$SkipPdfExtraction
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.IO.Compression.FileSystem

if (-not $OutputRun) {
    $OutputRun = Get-Date -Format 'yyyy-MM-dd_HHmmss'
}

if (-not $RepoRoot) {
    $RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
}

$OutputDir = Join-Path $RepoRoot "results\dna_quality\$OutputRun"
New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null

function Normalize-Whitespace {
    param([string]$Value)
    if ($null -eq $Value) { return '' }
    return (($Value -replace "`u00A0", ' ') -replace '\s+', ' ').Trim()
}

function Normalize-Key {
    param([string]$Value)
    return ((Normalize-Whitespace $Value).ToLowerInvariant() -replace '[^a-z0-9]+', '')
}

function Get-AliasKey {
    param([string]$Value)
    return ((Normalize-Whitespace $Value).ToUpperInvariant() -replace '[^A-Z0-9]+', '')
}

function Get-SampleIdHint {
    param([string]$Text)
    $text = Normalize-Whitespace $Text
    if (-not $text) { return '' }
    if ($text -match '(?i)(CJD\d+)(?=$|[^0-9A-Z]|_)') { return $Matches[1].ToUpperInvariant() }
    if ($text -match '(?i)(CTRL\d+)(?=$|[^0-9A-Z]|_)') { return $Matches[1].ToUpperInvariant() }
    if ($text -match '(?i)(NA[0-9A-Z_]+)(?=$|[^0-9A-Z]|_)') { return $Matches[1].ToUpperInvariant() }
    return ''
}

function Get-FirstPropertyValue {
    param(
        $Row,
        [string[]]$Candidates
    )

    foreach ($candidate in $Candidates) {
        $target = Normalize-Key $candidate
        foreach ($property in $Row.PSObject.Properties) {
            if ((Normalize-Key $property.Name) -eq $target) {
                $value = Normalize-Whitespace ([string]$property.Value)
                if ($value) { return $value }
            }
        }
    }
    return ''
}

function Get-ObjectPropertyValue {
    param(
        $Row,
        [string]$Field
    )

    if ($null -eq $Row) { return '' }
    $property = $Row.PSObject.Properties | Where-Object { $_.Name -eq $Field } | Select-Object -First 1
    if (-not $property) { return '' }
    return $property.Value
}

function Get-RunTokens {
    param([string]$RunId)
    if (-not $RunId) { return @() }
    $stopWords = @('SS', 'QC', 'SEQUENCE', 'SAMPLE', 'SAMPLES', 'PART', 'RUN', 'FINAL', 'FIRST', 'SECOND', 'THIRD', 'DAY', 'BATCH', 'OF', 'THE', 'MIRKA')
    return @($RunId.ToUpperInvariant().Split('_') | Where-Object { $_ -and $_ -notin $stopWords } | Select-Object -Unique)
}

function Get-RunOverlapScore {
    param(
        [string]$Left,
        [string]$Right
    )

    if (-not $Left -or -not $Right) { return 0 }
    $leftTokens = Get-RunTokens -RunId $Left
    $rightTokens = Get-RunTokens -RunId $Right
    $rightLookup = @{}
    foreach ($token in $rightTokens) { $rightLookup[$token] = $true }

    $score = 0
    foreach ($token in $leftTokens) {
        if ($rightLookup.ContainsKey($token)) {
            if ($token -match '^\d{4}$') { $score += 3 }
            elseif ($token -match '^\d+$') { $score += 2 }
            else { $score += 1 }
        }
    }
    return $score
}

function Get-StageForPath {
    param([string]$Path)
    $text = $Path.ToLowerInvariant()
    if ($text.Contains('1st qc') -or $text.Contains('first day')) { return 'pre_capture' }
    if ($text.Contains('post-capture')) { return 'post_capture' }
    if ($text.Contains('qc before sample submission') -or $text.Contains('sample submission')) { return 'submission_qc' }
    if ($text.Contains(' second sequencing run ') -or $text.EndsWith('second sequencing run')) { return 'submission_qc' }
    return 'unknown'
}

function Get-RunId {
    param(
        [string]$Path,
        [string]$Root
    )

    $parent = Split-Path -Parent $Path
    if ($parent.StartsWith($Root, [System.StringComparison]::OrdinalIgnoreCase)) {
        $relative = $parent.Substring($Root.Length).TrimStart('\')
        if ($relative) {
            return (($relative -replace '[^A-Za-z0-9]+', '_').Trim('_'))
        }
    }
    return (([System.IO.Path]::GetFileNameWithoutExtension($Path) -replace '[^A-Za-z0-9]+', '_').Trim('_'))
}

function Get-ZipEntryText {
    param(
        $Archive,
        [string]$EntryName
    )

    $entry = $Archive.GetEntry($EntryName)
    if (-not $entry) { return '' }
    $stream = $entry.Open()
    $reader = $null
    try {
        $reader = New-Object System.IO.StreamReader($stream)
        return $reader.ReadToEnd()
    }
    finally {
        if ($reader) { $reader.Dispose() } elseif ($stream) { $stream.Dispose() }
    }
}

function Get-SpreadsheetSharedStrings {
    param($Archive)

    $text = Get-ZipEntryText -Archive $Archive -EntryName 'xl/sharedStrings.xml'
    if (-not $text) { return @() }
    $xml = [xml]$text
    $values = New-Object System.Collections.Generic.List[string]
    foreach ($node in $xml.SelectNodes("//*[local-name()='si']")) {
        $parts = @($node.SelectNodes(".//*[local-name()='t']") | ForEach-Object { $_.InnerText })
        $values.Add((Normalize-Whitespace ($parts -join ''))) | Out-Null
    }
    return $values.ToArray()
}

function Get-SpreadsheetSheets {
    param($Archive)

    $workbookText = Get-ZipEntryText -Archive $Archive -EntryName 'xl/workbook.xml'
    $relsText = Get-ZipEntryText -Archive $Archive -EntryName 'xl/_rels/workbook.xml.rels'
    if (-not $workbookText -or -not $relsText) { return @() }

    $workbook = [xml]$workbookText
    $rels = [xml]$relsText
    $relMap = @{}
    foreach ($relationship in $rels.SelectNodes("//*[local-name()='Relationship']")) {
        $relMap[$relationship.Id] = $relationship.Target
    }

    $sheets = New-Object System.Collections.Generic.List[object]
    foreach ($sheet in $workbook.SelectNodes("//*[local-name()='sheet']")) {
        $rid = $sheet.GetAttribute('id', 'http://schemas.openxmlformats.org/officeDocument/2006/relationships')
        if (-not $rid) {
            $idAttr = $sheet.Attributes | Where-Object { $_.LocalName -eq 'id' } | Select-Object -First 1
            if ($idAttr) { $rid = $idAttr.Value }
        }
        $target = $relMap[$rid]
        if (-not $target) { continue }
        $entryPath = $target.TrimStart('/')
        if (-not $entryPath.StartsWith('xl/')) { $entryPath = "xl/$entryPath" }
        $sheets.Add([pscustomobject]@{
                Name = Normalize-Whitespace ([string]$sheet.name)
                EntryPath = $entryPath
            }) | Out-Null
    }
    return $sheets
}

function Get-CellColumnIndex {
    param([string]$Reference)
    if ($Reference -notmatch '^([A-Z]+)') { return 0 }
    $index = 0
    foreach ($letter in $Matches[1].ToCharArray()) {
        $index = ($index * 26) + ([int][char]$letter - [int][char]'A' + 1)
    }
    return $index
}

function Get-CellText {
    param(
        $Cell,
        [string[]]$SharedStrings
    )

    $cellType = ''
    if ($Cell.Attributes['t']) {
        $cellType = [string]$Cell.Attributes['t'].Value
    }

    if ($cellType -eq 'inlineStr') {
        $parts = @($Cell.SelectNodes(".//*[local-name()='t']") | ForEach-Object { $_.InnerText })
        return Normalize-Whitespace ($parts -join '')
    }

    $valueNode = $Cell.SelectSingleNode("*[local-name()='v']")
    if (-not $valueNode) { return '' }
    $raw = Normalize-Whitespace $valueNode.InnerText
    if ($cellType -eq 's' -and $raw -match '^\d+$') {
        $index = [int]$raw
        if ($index -lt $SharedStrings.Length) { return $SharedStrings[$index] }
    }
    return $raw
}

function Get-XlsxRows {
    param([string]$Path)

    $archive = [System.IO.Compression.ZipFile]::OpenRead($Path)
    try {
        $sharedStrings = Get-SpreadsheetSharedStrings -Archive $archive
        $sheets = Get-SpreadsheetSheets -Archive $archive
        $rows = New-Object System.Collections.Generic.List[object]

        foreach ($sheet in $sheets) {
            $sheetText = Get-ZipEntryText -Archive $archive -EntryName $sheet.EntryPath
            if (-not $sheetText) { continue }
            $sheetXml = [xml]$sheetText
            $headerNames = $null
            foreach ($row in $sheetXml.SelectNodes("//*[local-name()='sheetData']/*[local-name()='row']")) {
                $cellMap = @{}
                foreach ($cell in $row.SelectNodes("*[local-name()='c']")) {
                    $columnIndex = Get-CellColumnIndex -Reference ([string]$cell.r)
                    if ($columnIndex -le 0) { continue }
                    $cellMap[$columnIndex] = Get-CellText -Cell $cell -SharedStrings $sharedStrings
                }
                if ($cellMap.Count -eq 0) { continue }

                $maxColumn = ($cellMap.Keys | Measure-Object -Maximum).Maximum
                $values = for ($column = 1; $column -le $maxColumn; $column++) {
                    if ($cellMap.ContainsKey($column)) { $cellMap[$column] } else { '' }
                }

                if (-not $headerNames) {
                    $headerNames = @()
                    $seen = @{}
                    for ($index = 0; $index -lt (@($values)).Count; $index++) {
                        $header = Normalize-Whitespace $values[$index]
                        if (-not $header) { $header = "column_$($index + 1)" }
                        $candidate = $header
                        $suffix = 2
                        while ($seen.ContainsKey($candidate.ToLowerInvariant())) {
                            $candidate = "{0}_{1}" -f $header, $suffix
                            $suffix++
                        }
                        $seen[$candidate.ToLowerInvariant()] = $true
                        $headerNames += $candidate
                    }
                    continue
                }

                if ((@($values | Where-Object { $_ })).Count -eq 0) { continue }
                $record = [ordered]@{
                    source_path = $Path
                    sheet_name = $sheet.Name
                    row_index = [int]$row.r
                }
                for ($index = 0; $index -lt (@($headerNames)).Count; $index++) {
                    $record[$headerNames[$index]] = if ($index -lt (@($values)).Count) { $values[$index] } else { '' }
                }
                $rows.Add([pscustomobject]$record) | Out-Null
            }
        }
        return $rows
    }
    finally {
        $archive.Dispose()
    }
}

function Write-Tsv {
    param(
        [string]$Path,
        $Rows,
        [string[]]$Columns
    )

    $materializedRows = @($Rows | ForEach-Object { $_ })
    if ((@($materializedRows)).Count -eq 0) {
        Set-Content -Path $Path -Value ($Columns -join "`t")
        return
    }

    $orderedRows = foreach ($row in $materializedRows) {
        $record = [ordered]@{}
        foreach ($column in $Columns) {
            $value = Get-ObjectPropertyValue -Row $row -Field $column
            $record[$column] = if ($null -ne $value) { $value } else { '' }
        }
        [pscustomobject]$record
    }
    $orderedRows | Export-Csv -Path $Path -Delimiter "`t" -Encoding UTF8 -NoTypeInformation
}

function Convert-ToNumberOrBlank {
    param([string]$Value)
    $text = Normalize-Whitespace $Value
    if (-not $text) { return '' }
    $text = $text -replace ',', '.'
    [double]$number = 0
    if ([double]::TryParse($text, [System.Globalization.NumberStyles]::Float, [System.Globalization.CultureInfo]::InvariantCulture, [ref]$number)) {
        return $number
    }
    return ''
}

function Get-PdfLines {
    param([string]$Path)
    if ($SkipPdfExtraction) { return @() }
    $raw = & gswin64c -q -dNOPAUSE -dBATCH -sDEVICE=txtwrite "-sOutputFile=-" $Path
    return @($raw -split "`r?`n")
}

function Get-RunNativePair {
    param([System.IO.FileInfo[]]$Files)
    $native = $Files | Where-Object { $_.Extension -in @('.D1000', '.HSD1000') } | Select-Object -First 1
    if ($native) { return $native.FullName }
    return ''
}

function Parse-TapeStationCsvBundle {
    param(
        [System.IO.DirectoryInfo]$Directory,
        [string]$RunId,
        [string]$Stage
    )

    $sampleCsv = Get-ChildItem -Path $Directory.FullName -Filter '*_sampleTable.csv' -File | Select-Object -First 1
    $peakCsv = Get-ChildItem -Path $Directory.FullName -Filter '*_compactPeakTable.csv' -File | Select-Object -First 1
    $regionCsv = Get-ChildItem -Path $Directory.FullName -Filter '*_compactRegionTable.csv' -File | Select-Object -First 1
    $pdfPair = Get-ChildItem -Path $Directory.FullName -Filter '*.pdf' -File | Select-Object -First 1
    $nativePair = Get-RunNativePair -Files (Get-ChildItem -Path $Directory.FullName -File)

    if (-not $sampleCsv) { return @() }

    $samples = Import-Csv -Path $sampleCsv.FullName
    $peaks = if ($peakCsv) { Import-Csv -Path $peakCsv.FullName } else { @() }
    $regions = if ($regionCsv) { Import-Csv -Path $regionCsv.FullName } else { @() }

    $peakMap = @{}
    foreach ($peak in $peaks | Where-Object { $_.Well -ne 'EL1' -and $_.'Peak Comment' -notmatch 'Marker' }) {
        $key = $peak.Well
        if (-not $peakMap.ContainsKey($key)) { $peakMap[$key] = @() }
        $peakMap[$key] += $peak
    }

    $regionMap = @{}
    foreach ($region in $regions) {
        $regionMap[$region.WellId] = $region
    }

    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($sample in $samples | Where-Object { $_.Well -ne 'EL1' }) {
        $peakRows = if ($peakMap.ContainsKey($sample.Well)) { $peakMap[$sample.Well] } else { @() }
        $dominantPeak = $peakRows | Sort-Object { [double]($_.'% Integrated Area') } -Descending | Select-Object -First 1
        $region = if ($regionMap.ContainsKey($sample.Well)) { $regionMap[$sample.Well] } else { $null }
        $rows.Add([pscustomobject]@{
                run_id = $RunId
                stage = $Stage
                instrument = 'D1000'
                parser_source = 'csv'
                source_path = $sampleCsv.FullName
                native_source_path = $nativePair
                pdf_source_path = if ($pdfPair) { $pdfPair.FullName } else { '' }
                well = Normalize-Whitespace $sample.Well
                sample_description = Normalize-Whitespace $sample.'Sample Description'
                reported_concentration = Convert-ToNumberOrBlank (Get-FirstPropertyValue -Row $sample -Candidates @('Conc. [ng/µl]', 'Conc. [ng/ul]'))
                reported_concentration_unit = 'ng/ul'
                peak_row_count = (@($peakRows)).Count
                dominant_peak_size_bp = if ($dominantPeak) { Convert-ToNumberOrBlank $dominantPeak.'Size [bp]' } else { '' }
                dominant_peak_calibrated_concentration = if ($dominantPeak) { Convert-ToNumberOrBlank (Get-FirstPropertyValue -Row $dominantPeak -Candidates @('Calibrated Conc. [ng/µl]', 'Calibrated Conc. [ng/ul]')) } else { '' }
                dominant_peak_molarity = if ($dominantPeak) { Convert-ToNumberOrBlank $dominantPeak.'Peak Molarity [nmol/l]' } else { '' }
                dominant_peak_area_pct = if ($dominantPeak) { Convert-ToNumberOrBlank $dominantPeak.'% Integrated Area' } else { '' }
                region_from_bp = if ($region) { Convert-ToNumberOrBlank $region.'From [bp]' } else { '' }
                region_to_bp = if ($region) { Convert-ToNumberOrBlank $region.'To [bp]' } else { '' }
                region_average_size_bp = if ($region) { Convert-ToNumberOrBlank $region.'Average Size [bp]' } else { '' }
                region_concentration = if ($region) { Convert-ToNumberOrBlank (Get-FirstPropertyValue -Row $region -Candidates @('Conc. [ng/µl]', 'Conc. [ng/ul]')) } else { '' }
                region_concentration_unit = 'ng/ul'
                region_molarity = if ($region) { Convert-ToNumberOrBlank $region.'Region Molarity [nmol/l]' } else { '' }
                region_molarity_unit = 'nmol/l'
                region_percent_total = if ($region) { Convert-ToNumberOrBlank $region.'% of Total' } else { '' }
                alert = Normalize-Whitespace $sample.Alert
                observations = Normalize-Whitespace $sample.Observations
            }) | Out-Null
    }
    return $rows
}

function Parse-TapeStationPdf {
    param(
        [string]$Path,
        [string]$RunId,
        [string]$Stage
    )

    $lines = Get-PdfLines -Path $Path
    if (-not $lines -or (@($lines)).Count -eq 0) { return @() }
    $instrument = if (($lines | Where-Object { $_ -match 'High Sensitivity D1000' } | Select-Object -First 1)) { 'HSD1000' } else { 'D1000' }
    $concUnit = if ($instrument -eq 'HSD1000') { 'pg/ul' } else { 'ng/ul' }
    $molarityUnit = if ($instrument -eq 'HSD1000') { 'pmol/l' } else { 'nmol/l' }
    $nativePair = Get-RunNativePair -Files (Get-ChildItem -Path (Split-Path -Parent $Path) -File)
    $rows = New-Object System.Collections.Generic.List[object]
    $current = $null
    $section = ''

    foreach ($line in $lines) {
        $trimmed = Normalize-Whitespace $line
        if (-not $trimmed) { continue }

        if ($trimmed -match '^(EL1|[A-H]\d{1,2}):\s*(.+)$') {
            if ($current -and $current.well -ne 'EL1') {
                $rows.Add([pscustomobject]$current) | Out-Null
            }
            $current = [ordered]@{
                run_id = $RunId
                stage = $Stage
                instrument = $instrument
                parser_source = 'pdf'
                source_path = $Path
                native_source_path = $nativePair
                pdf_source_path = $Path
                well = $Matches[1]
                sample_description = Normalize-Whitespace $Matches[2]
                reported_concentration = ''
                reported_concentration_unit = $concUnit
                peak_row_count = 0
                dominant_peak_size_bp = ''
                dominant_peak_calibrated_concentration = ''
                dominant_peak_molarity = ''
                dominant_peak_area_pct = ''
                region_from_bp = ''
                region_to_bp = ''
                region_average_size_bp = ''
                region_concentration = ''
                region_concentration_unit = $concUnit
                region_molarity = ''
                region_molarity_unit = $molarityUnit
                region_percent_total = ''
                alert = ''
                observations = ''
            }
            $section = ''
            continue
        }

        if (-not $current) { continue }
        if ($trimmed -eq 'Sample Table') { $section = 'sample'; continue }
        if ($trimmed -eq 'Peak Table') { $section = 'peak'; continue }
        if ($trimmed -eq 'Region Table') { $section = 'region'; continue }
        if ($trimmed -like 'TapeStation Analysis Software*' -or $trimmed -like '*ScreenTape*' -or $trimmed -like 'Filename:*' -or $trimmed -eq 'Sample Info' -or $trimmed -eq 'Default image (Contrast 100%)' -or $trimmed -eq 'Default image (Contrast 50%), Image is Scaled to Sample') { continue }

        switch ($section) {
            'sample' {
                if ($trimmed -match "^(EL1|[A-H]\d{1,2})\s+(.*)$" -and $Matches[1] -eq $current.well) {
                    $tail = Normalize-Whitespace $Matches[2]
                    if ($tail -match '^([0-9]+(?:\.[0-9]+)?)\s+') {
                        $current.reported_concentration = Convert-ToNumberOrBlank $Matches[1]
                    }
                    if ($tail -match '(Caution!.*)$') {
                        $current.observations = Normalize-Whitespace $Matches[1]
                    }
                }
            }
            'peak' {
                if ($trimmed -match 'Marker') { continue }
                $tokens = @($trimmed -split '\s+')
                if ($instrument -eq 'D1000' -and (@($tokens)).Count -ge 5 -and $tokens[0] -match '^\d+(\.\d+)?$' -and $tokens[4] -match '^\d+(\.\d+)?$') {
                    $area = [double](Convert-ToNumberOrBlank $tokens[4])
                    $current.peak_row_count = [int]$current.peak_row_count + 1
                    if (-not $current.dominant_peak_area_pct -or $area -gt [double]$current.dominant_peak_area_pct) {
                        $current.dominant_peak_size_bp = Convert-ToNumberOrBlank $tokens[0]
                        $current.dominant_peak_calibrated_concentration = Convert-ToNumberOrBlank $tokens[1]
                        $current.dominant_peak_molarity = Convert-ToNumberOrBlank $tokens[3]
                        $current.dominant_peak_area_pct = $area
                    }
                }
                elseif ($instrument -eq 'HSD1000' -and (@($tokens)).Count -ge 4 -and $tokens[0] -match '^\d+(\.\d+)?$' -and $tokens[3] -match '^\d+(\.\d+)?$') {
                    $area = [double](Convert-ToNumberOrBlank $tokens[3])
                    $current.peak_row_count = [int]$current.peak_row_count + 1
                    if (-not $current.dominant_peak_area_pct -or $area -gt [double]$current.dominant_peak_area_pct) {
                        $current.dominant_peak_calibrated_concentration = Convert-ToNumberOrBlank $tokens[0]
                        $current.dominant_peak_molarity = Convert-ToNumberOrBlank $tokens[2]
                        $current.dominant_peak_area_pct = $area
                    }
                }
            }
            'region' {
                if ($trimmed -match '^(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)') {
                    $current.region_from_bp = Convert-ToNumberOrBlank $Matches[1]
                    $current.region_to_bp = Convert-ToNumberOrBlank $Matches[2]
                    $current.region_average_size_bp = Convert-ToNumberOrBlank $Matches[3]
                    $current.region_concentration = Convert-ToNumberOrBlank $Matches[4]
                    $current.region_molarity = Convert-ToNumberOrBlank $Matches[5]
                    $current.region_percent_total = Convert-ToNumberOrBlank $Matches[6]
                }
            }
        }
    }

    if ($current -and $current.well -ne 'EL1') {
        $rows.Add([pscustomobject]$current) | Out-Null
    }
    return $rows
}

function Get-MetadataRows {
    param(
        [System.IO.FileInfo[]]$Files,
        [string]$SourceFamily
    )

    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($file in $Files) {
        $sourceRunId = switch ($SourceFamily) {
            'sureselect_prep' { Get-RunId -Path $file.FullName -Root $SureSelectRoot }
            'ddpcr_quantity' { Get-RunId -Path $file.FullName -Root $DdPcrRoot }
            'sample_manifest' { Get-RunId -Path $file.FullName -Root $SamplesRoot }
            default { '' }
        }
        $workbookRows = Get-XlsxRows -Path $file.FullName
        foreach ($row in $workbookRows) {
            $samplePieces = @(
                (Get-FirstPropertyValue -Row $row -Candidates @('new_name', 'new name', 'sample name', 'sample', '# sample / code', 'sample code', 'name', 'patient')),
                (Get-FirstPropertyValue -Row $row -Candidates @('code')),
                (Get-FirstPropertyValue -Row $row -Candidates @('sample description'))
            ) | Where-Object { $_ }

            $sampleId = ''
            foreach ($piece in $samplePieces) {
                $candidate = Get-SampleIdHint -Text $piece
                if ($candidate) { $sampleId = $candidate; break }
            }
            if (-not $sampleId) {
                foreach ($property in $row.PSObject.Properties) {
                    $candidate = Get-SampleIdHint -Text ([string]$property.Value)
                    if ($candidate) { $sampleId = $candidate; break }
                }
            }

            $record = [ordered]@{
                source_family = $SourceFamily
                source_run_id = $sourceRunId
                source_path = $file.FullName
                sheet_name = $row.sheet_name
                row_index = $row.row_index
                raw_sample_label = (Get-FirstPropertyValue -Row $row -Candidates @('new_name', 'new name', 'sample name', 'sample', '# sample / code', 'sample code', 'name', 'patient'))
                sample_id = $sampleId
                code = (Get-FirstPropertyValue -Row $row -Candidates @('code'))
                group = (Get-FirstPropertyValue -Row $row -Candidates @('group'))
                region = (Get-FirstPropertyValue -Row $row -Candidates @('brain region', 'region'))
                well_or_barcode = (Get-FirstPropertyValue -Row $row -Candidates @('sureselect barcode', 'barcode', 'index', 'well'))
                dna_input_ng_ul = Convert-ToNumberOrBlank (Get-FirstPropertyValue -Row $row -Candidates @('DNA fr (ng/ul)', 'DNA (ng/ul)', 'DNA (ng/µl)', 'DNA ng/ul'))
                clean_dna_ng_ul = Convert-ToNumberOrBlank (Get-FirstPropertyValue -Row $row -Candidates @('clean DNA (ng/ul)', 'clean DNA'))
                qubit_ng_ul = Convert-ToNumberOrBlank (Get-FirstPropertyValue -Row $row -Candidates @('qubit', 'qubit conc. (ng/ul)', 'qubit conc', 'DNA conc qubit'))
                pre_capture_pcr_ng_ul = Convert-ToNumberOrBlank (Get-FirstPropertyValue -Row $row -Candidates @('pre-capt. PCR product conc. (ng/ul)', 'pre capt PCR product conc. (ng/ul)'))
                molarity_nm = Convert-ToNumberOrBlank (Get-FirstPropertyValue -Row $row -Candidates @('final sample conc. / molarity [nm]', 'molarity', 'molarity (nm)'))
                elution_volume_ul = Convert-ToNumberOrBlank (Get-FirstPropertyValue -Row $row -Candidates @('eluted in (µl)', 'eluted in (ul)'))
                notes = (Get-FirstPropertyValue -Row $row -Candidates @('comments', 'comment', 'notes', 'observations'))
            }

            if ($record.raw_sample_label -or $record.sample_id -or $record.well_or_barcode -or $record.dna_input_ng_ul -or $record.clean_dna_ng_ul -or $record.qubit_ng_ul -or $record.pre_capture_pcr_ng_ul -or $record.molarity_nm) {
                $rows.Add([pscustomobject]$record) | Out-Null
            }
        }
    }
    return $rows
}

function Get-MaxNumericValue {
    param(
        $Rows,
        [string]$Field
    )
    $values = @($Rows | ForEach-Object {
            $value = Get-ObjectPropertyValue -Row $_ -Field $Field
            if ($value -ne '') { [double]$value }
        })
    if ((@($values)).Count -eq 0) { return '' }
    return ($values | Measure-Object -Maximum).Maximum
}

function Get-MeanNumericValue {
    param(
        $Rows,
        [string]$Field
    )
    $values = @($Rows | ForEach-Object {
            $value = Get-ObjectPropertyValue -Row $_ -Field $Field
            if ($value -ne '') { [double]$value }
        })
    if ((@($values)).Count -eq 0) { return '' }
    return [Math]::Round((($values | Measure-Object -Average).Average), 3)
}

function Get-MedianNumericValue {
    param(
        $Rows,
        [string]$Field
    )
    $values = @($Rows | ForEach-Object {
            $value = Get-ObjectPropertyValue -Row $_ -Field $Field
            if ($value -ne '') { [double]$value }
        } | Sort-Object)
    if ((@($values)).Count -eq 0) { return '' }
    $count = $values.Count
    if ($count % 2 -eq 1) { return $values[[int][Math]::Floor($count / 2)] }
    return [Math]::Round((($values[($count / 2) - 1] + $values[$count / 2]) / 2), 3)
}

function Get-DistinctJoinedValue {
    param(
        $Rows,
        [string]$Field
    )
    $values = @($Rows | ForEach-Object { Normalize-Whitespace ([string](Get-ObjectPropertyValue -Row $_ -Field $Field)) } | Where-Object { $_ } | Sort-Object -Unique)
    return ($values -join '; ')
}

function Resolve-AliasMatch {
    param(
        [System.Collections.Generic.List[object]]$AliasEntries,
        [string[]]$CandidateValues,
        [string]$RunId
    )

    foreach ($candidateValue in $CandidateValues) {
        $key = Get-AliasKey $candidateValue
        if (-not $key) { continue }

        $matches = @($AliasEntries | Where-Object { $_.alias_key -eq $key })
        if ($matches.Count -eq 0) { continue }

        $uniqueSampleIds = @($matches | Select-Object -ExpandProperty sample_id -Unique)
        if ($uniqueSampleIds.Count -eq 1) {
            return [pscustomobject]@{
                sample_id = $uniqueSampleIds[0]
                source = $matches[0].source
            }
        }

        $bestMatch = $null
        $bestScore = -1
        $bestSampleId = ''
        $tie = $false
        foreach ($match in $matches) {
            $score = Get-RunOverlapScore -Left $RunId -Right $match.source_run_id
            if ($score -gt $bestScore) {
                $bestScore = $score
                $bestMatch = $match
                $bestSampleId = $match.sample_id
                $tie = $false
            }
            elseif ($score -eq $bestScore -and $match.sample_id -ne $bestSampleId) {
                $tie = $true
            }
        }

        if ($bestMatch -and -not $tie -and $bestScore -gt 0) {
            return [pscustomobject]@{
                sample_id = $bestMatch.sample_id
                source = "{0}_context" -f $bestMatch.source
            }
        }
    }

    return $null
}

$sureSelectFiles = @(Get-ChildItem -Path $SureSelectRoot -Recurse -File | Where-Object {
        $_.Extension -in @('.pdf', '.csv', '.xlsx', '.D1000', '.HSD1000')
    })
$sureSelectXlsx = @($sureSelectFiles | Where-Object { $_.Extension -eq '.xlsx' })
$csvBundleDirectories = @($sureSelectFiles | Where-Object { $_.Name -like '*_sampleTable.csv' } | Select-Object -ExpandProperty DirectoryName -Unique)
$pdfFilesForParsing = @($sureSelectFiles | Where-Object { $_.Extension -eq '.pdf' -and $_.DirectoryName -notin $csvBundleDirectories })

$ddpcrFiles = @(Get-ChildItem -Path $DdPcrRoot -Recurse -File | Where-Object {
        $_.Extension -in @('.xlsx', '.docx') -and $_.Name -match '(?i)Covaris and Tapestation|summary table patients ddPCR|qubit|clean DNA|qPCR - samples'
    })
$ddpcrXlsx = @($ddpcrFiles | Where-Object { $_.Extension -eq '.xlsx' })
$sampleFiles = @(Get-ChildItem -Path $SamplesRoot -Recurse -File | Where-Object { $_.Extension -eq '.xlsx' })

$fileInventory = New-Object System.Collections.Generic.List[object]
foreach ($file in $sureSelectFiles) {
    $fileInventory.Add([pscustomobject]@{
            source_family = 'sureselect'
            run_id = Get-RunId -Path $file.FullName -Root $SureSelectRoot
            stage = Get-StageForPath -Path $file.FullName
            instrument = if ($file.Extension -eq '.HSD1000') { 'HSD1000' } elseif ($file.Extension -eq '.D1000' -or $file.Name -match 'D1000') { 'D1000' } else { '' }
            file_type = $file.Extension.TrimStart('.')
            parser_role = if ($file.Name -like '*_sampleTable.csv') { 'sample_table' } elseif ($file.Name -like '*_compactPeakTable.csv') { 'peak_table' } elseif ($file.Name -like '*_compactRegionTable.csv') { 'region_table' } elseif ($file.Extension -eq '.pdf') { 'report' } elseif ($file.Extension -in @('.D1000', '.HSD1000')) { 'native' } else { 'metadata' }
            parse_strategy = if ($file.Extension -eq '.pdf' -and $file.DirectoryName -notin $csvBundleDirectories) { 'pdf_text' } elseif ($file.Extension -eq '.csv') { 'csv_bundle' } elseif ($file.Extension -in @('.D1000', '.HSD1000')) { 'inventory_only' } else { 'spreadsheet' }
            path = $file.FullName
        }) | Out-Null
}
foreach ($file in $ddpcrFiles) {
    $fileInventory.Add([pscustomobject]@{
            source_family = 'ddpcr'
            run_id = Get-RunId -Path $file.FullName -Root $DdPcrRoot
            stage = 'upstream_quantity'
            instrument = if ($file.Extension -eq '.docx') { 'documentation' } else { 'spreadsheet' }
            file_type = $file.Extension.TrimStart('.')
            parser_role = if ($file.Extension -eq '.docx') { 'protocol_note' } else { 'metadata' }
            parse_strategy = if ($file.Extension -eq '.xlsx') { 'spreadsheet' } else { 'inventory_only' }
            path = $file.FullName
        }) | Out-Null
}
foreach ($file in $sampleFiles) {
    $fileInventory.Add([pscustomobject]@{
            source_family = 'samples'
            run_id = Get-RunId -Path $file.FullName -Root $SamplesRoot
            stage = 'sample_manifest'
            instrument = 'spreadsheet'
            file_type = $file.Extension.TrimStart('.')
            parser_role = 'metadata'
            parse_strategy = 'spreadsheet'
            path = $file.FullName
        }) | Out-Null
}

$libraryRows = New-Object System.Collections.Generic.List[object]
foreach ($directoryName in $csvBundleDirectories) {
    $directory = Get-Item -Path $directoryName
    $runId = Get-RunId -Path $directory.FullName -Root $SureSelectRoot
    $stage = Get-StageForPath -Path $directory.FullName
    foreach ($row in Parse-TapeStationCsvBundle -Directory $directory -RunId $runId -Stage $stage) {
        $libraryRows.Add($row) | Out-Null
    }
}
foreach ($pdf in $pdfFilesForParsing) {
    $runId = Get-RunId -Path $pdf.FullName -Root $SureSelectRoot
    $stage = Get-StageForPath -Path $pdf.FullName
    foreach ($row in Parse-TapeStationPdf -Path $pdf.FullName -RunId $runId -Stage $stage) {
        $libraryRows.Add($row) | Out-Null
    }
}

$prepMetadata = New-Object System.Collections.Generic.List[object]
foreach ($row in Get-MetadataRows -Files $sureSelectXlsx -SourceFamily 'sureselect_prep') {
    $prepMetadata.Add($row) | Out-Null
}
foreach ($row in Get-MetadataRows -Files $ddpcrXlsx -SourceFamily 'ddpcr_quantity') {
    $prepMetadata.Add($row) | Out-Null
}
foreach ($row in Get-MetadataRows -Files $sampleFiles -SourceFamily 'sample_manifest') {
    $prepMetadata.Add($row) | Out-Null
}

$aliasEntries = New-Object System.Collections.Generic.List[object]
$aliasMap = @{}
function Add-AliasEntry {
    param(
        [string]$Alias,
        [string]$SampleId,
        [string]$Source,
        [string]$SourceRunId,
        [string]$Notes
    )
    $key = Get-AliasKey $Alias
    if (-not $key -or -not $SampleId) { return }
    if ($aliasMap.ContainsKey($key)) {
        if ($aliasMap[$key].sample_id -ne $SampleId) { $aliasMap[$key].ambiguous = $true }
    }
    else {
        $aliasMap[$key] = [ordered]@{
            sample_id = $SampleId
            source = $Source
            ambiguous = $false
        }
    }
    $aliasEntries.Add([pscustomobject]@{
            alias = Normalize-Whitespace $Alias
            alias_key = $key
            sample_id = $SampleId
            source = $Source
            source_run_id = $SourceRunId
            notes = Normalize-Whitespace $Notes
        }) | Out-Null
}

foreach ($row in $prepMetadata | Where-Object { $_.sample_id }) {
    Add-AliasEntry -Alias $row.sample_id -SampleId $row.sample_id -Source $row.source_family -SourceRunId $row.source_run_id -Notes $row.source_path
    Add-AliasEntry -Alias $row.raw_sample_label -SampleId $row.sample_id -Source $row.source_family -SourceRunId $row.source_run_id -Notes $row.source_path
    Add-AliasEntry -Alias $row.code -SampleId $row.sample_id -Source $row.source_family -SourceRunId $row.source_run_id -Notes $row.source_path
    Add-AliasEntry -Alias $row.well_or_barcode -SampleId $row.sample_id -Source $row.source_family -SourceRunId $row.source_run_id -Notes $row.source_path
}

$preprocessingManifestPath = Join-Path $RepoRoot 'config\preprocessing_samples.tsv'
$preprocessingRows = if (Test-Path $preprocessingManifestPath) { Import-Csv -Path $preprocessingManifestPath -Delimiter "`t" } else { @() }
$batchMap = @{}
foreach ($row in $preprocessingRows) {
    $batchMap[$row.sample_id] = $row.batch_id
    Add-AliasEntry -Alias $row.sample_id -SampleId $row.sample_id -Source 'config_preprocessing' -SourceRunId '' -Notes $row.batch_id
}

$overridePath = Join-Path $RepoRoot 'config\dna_quality_sample_alias_overrides.tsv'
if (Test-Path $overridePath) {
    $overrideRows = Import-Csv -Path $overridePath -Delimiter "`t"
    foreach ($row in $overrideRows) {
        Add-AliasEntry -Alias $row.alias -SampleId $row.sample_id -Source 'alias_override' -SourceRunId '' -Notes $row.notes
    }
}

$resolvedPrepMetadata = New-Object System.Collections.Generic.List[object]
foreach ($row in $prepMetadata) {
    $sampleId = $row.sample_id
    if (-not $sampleId) {
        $directCandidates = @($row.raw_sample_label, $row.code, $row.well_or_barcode)
        foreach ($candidate in $directCandidates) {
            $directSampleId = Get-SampleIdHint -Text $candidate
            if ($directSampleId) {
                $sampleId = $directSampleId
                break
            }
        }
    }
    if (-not $sampleId) {
        $aliasMatch = Resolve-AliasMatch -AliasEntries $aliasEntries -CandidateValues @($row.raw_sample_label, $row.code, $row.well_or_barcode) -RunId $row.source_run_id
        if ($aliasMatch) {
            $sampleId = $aliasMatch.sample_id
        }
    }

    $resolvedPrepMetadata.Add([pscustomobject]@{
            source_family = $row.source_family
            source_run_id = $row.source_run_id
            source_path = $row.source_path
            sheet_name = $row.sheet_name
            row_index = $row.row_index
            raw_sample_label = $row.raw_sample_label
            sample_id = $sampleId
            code = $row.code
            group = $row.group
            region = $row.region
            well_or_barcode = $row.well_or_barcode
            dna_input_ng_ul = $row.dna_input_ng_ul
            clean_dna_ng_ul = $row.clean_dna_ng_ul
            qubit_ng_ul = $row.qubit_ng_ul
            pre_capture_pcr_ng_ul = $row.pre_capture_pcr_ng_ul
            molarity_nm = $row.molarity_nm
            elution_volume_ul = $row.elution_volume_ul
            notes = $row.notes
        }) | Out-Null
}

$resolvedLibraryRows = New-Object System.Collections.Generic.List[object]
foreach ($row in $libraryRows) {
    $sampleId = Get-SampleIdHint -Text $row.sample_description
    $sampleIdSource = ''
    if ($sampleId) {
        $sampleIdSource = 'inline_label'
    }
    else {
        $aliasMatch = Resolve-AliasMatch -AliasEntries $aliasEntries -CandidateValues @($row.sample_description, $row.well) -RunId $row.run_id
        if ($aliasMatch) {
            $sampleId = $aliasMatch.sample_id
            $sampleIdSource = $aliasMatch.source
        }
    }

    $resolvedLibraryRows.Add([pscustomobject]@{
            run_id = $row.run_id
            stage = if ($row.stage -eq 'unknown' -and $row.instrument -eq 'HSD1000') { 'submission_qc' } else { $row.stage }
            instrument = $row.instrument
            parser_source = $row.parser_source
            source_path = $row.source_path
            native_source_path = $row.native_source_path
            pdf_source_path = $row.pdf_source_path
            well = $row.well
            sample_description = $row.sample_description
            sample_id = $sampleId
            sample_id_source = $sampleIdSource
            batch_id = if ($sampleId -and $batchMap.ContainsKey($sampleId)) { $batchMap[$sampleId] } else { '' }
            reported_concentration = $row.reported_concentration
            reported_concentration_unit = $row.reported_concentration_unit
            peak_row_count = $row.peak_row_count
            dominant_peak_size_bp = $row.dominant_peak_size_bp
            dominant_peak_calibrated_concentration = $row.dominant_peak_calibrated_concentration
            dominant_peak_molarity = $row.dominant_peak_molarity
            dominant_peak_area_pct = $row.dominant_peak_area_pct
            region_from_bp = $row.region_from_bp
            region_to_bp = $row.region_to_bp
            region_average_size_bp = $row.region_average_size_bp
            region_concentration = $row.region_concentration
            region_concentration_unit = $row.region_concentration_unit
            region_molarity = $row.region_molarity
            region_molarity_unit = $row.region_molarity_unit
            region_percent_total = $row.region_percent_total
            alert = $row.alert
            observations = $row.observations
        }) | Out-Null
}

$prepBySample = @{}
foreach ($row in $resolvedPrepMetadata | Where-Object { $_.sample_id }) {
    if (-not $prepBySample.ContainsKey($row.sample_id)) { $prepBySample[$row.sample_id] = @() }
    $prepBySample[$row.sample_id] += $row
}

$seqFiles = @(Get-ChildItem -Path (Join-Path $RepoRoot 'results') -Recurse -Filter 'sequencing_metrics_per_sample.tsv' -File | Sort-Object LastWriteTime -Descending)
$preprocessingSampleIds = @($preprocessingRows | Select-Object -ExpandProperty sample_id -Unique)
$sequencingBySample = @{}
foreach ($file in $seqFiles) {
    foreach ($row in Import-Csv -Path $file.FullName -Delimiter "`t") {
        if ($row.sample_id -notin $preprocessingSampleIds) { continue }
        if (-not $sequencingBySample.ContainsKey($row.sample_id)) {
            $row | Add-Member -NotePropertyName source_path -NotePropertyValue $file.FullName
            $sequencingBySample[$row.sample_id] = $row
        }
    }
}

$inputDnaQuantity = @($resolvedPrepMetadata | Where-Object {
        $_.dna_input_ng_ul -ne '' -or $_.clean_dna_ng_ul -ne '' -or $_.qubit_ng_ul -ne '' -or $_.pre_capture_pcr_ng_ul -ne '' -or $_.molarity_nm -ne ''
    })

$sampleQualityMaster = New-Object System.Collections.Generic.List[object]
foreach ($row in $resolvedLibraryRows) {
    $prepRows = if ($row.sample_id -and $prepBySample.ContainsKey($row.sample_id)) { $prepBySample[$row.sample_id] } else { @() }
    $seq = if ($row.sample_id -and $sequencingBySample.ContainsKey($row.sample_id)) { $sequencingBySample[$row.sample_id] } else { $null }

    $libraryScore = 0
    if ($row.reported_concentration -ne '' -or $row.region_concentration -ne '') { $libraryScore++ }
    if ($row.region_average_size_bp -ne '' -and [double]$row.region_average_size_bp -ge 300 -and [double]$row.region_average_size_bp -le 500) { $libraryScore++ }
    if ($row.region_percent_total -ne '' -and [double]$row.region_percent_total -ge 85) { $libraryScore++ }
    if ($row.peak_row_count -ne '' -and [int]$row.peak_row_count -le 2) { $libraryScore++ }

    $sequencingBand = ''
    if ($seq) {
        $meanDepth = Convert-ToNumberOrBlank $seq.target_mean_depth
        $onTarget = Convert-ToNumberOrBlank $seq.on_target_fraction
        $fold80 = Convert-ToNumberOrBlank $seq.target_fold80
        if ($meanDepth -ne '' -and $onTarget -ne '' -and $fold80 -ne '' -and $meanDepth -ge 3000 -and $onTarget -ge 0.70 -and $fold80 -le 1.50) {
            $sequencingBand = 'high'
        }
        elseif ($meanDepth -ne '' -and $onTarget -ne '' -and $meanDepth -ge 2000 -and $onTarget -ge 0.60) {
            $sequencingBand = 'moderate'
        }
        else {
            $sequencingBand = 'low'
        }
    }

    $sampleQualityMaster.Add([pscustomobject]@{
            sample_id = $row.sample_id
            sample_id_source = $row.sample_id_source
            batch_id = $row.batch_id
            run_id = $row.run_id
            stage = $row.stage
            instrument = $row.instrument
            parser_source = $row.parser_source
            source_path = $row.source_path
            well = $row.well
            sample_description = $row.sample_description
            reported_concentration = $row.reported_concentration
            reported_concentration_unit = $row.reported_concentration_unit
            peak_row_count = $row.peak_row_count
            dominant_peak_size_bp = $row.dominant_peak_size_bp
            dominant_peak_area_pct = $row.dominant_peak_area_pct
            region_average_size_bp = $row.region_average_size_bp
            region_concentration = $row.region_concentration
            region_concentration_unit = $row.region_concentration_unit
            region_molarity = $row.region_molarity
            region_percent_total = $row.region_percent_total
            library_qc_heuristic_score = $libraryScore
            dna_input_ng_ul = Get-MaxNumericValue -Rows $prepRows -Field 'dna_input_ng_ul'
            clean_dna_ng_ul = Get-MaxNumericValue -Rows $prepRows -Field 'clean_dna_ng_ul'
            qubit_ng_ul = Get-MaxNumericValue -Rows $prepRows -Field 'qubit_ng_ul'
            pre_capture_pcr_ng_ul = Get-MaxNumericValue -Rows $prepRows -Field 'pre_capture_pcr_ng_ul'
            prep_molarity_nm = Get-MaxNumericValue -Rows $prepRows -Field 'molarity_nm'
            prep_source_families = Get-DistinctJoinedValue -Rows $prepRows -Field 'source_family'
            prep_regions = Get-DistinctJoinedValue -Rows $prepRows -Field 'region'
            sequencing_metrics_path = if ($seq) { $seq.source_path } else { '' }
            sequencing_on_target_fraction = if ($seq) { Convert-ToNumberOrBlank $seq.on_target_fraction } else { '' }
            sequencing_target_mean_depth = if ($seq) { Convert-ToNumberOrBlank $seq.target_mean_depth } else { '' }
            sequencing_target_fold80 = if ($seq) { Convert-ToNumberOrBlank $seq.target_fold80 } else { '' }
            sequencing_pct_duplication = if ($seq) { Convert-ToNumberOrBlank $seq.pct_duplication } else { '' }
            sequencing_estimated_library_size = if ($seq) { Convert-ToNumberOrBlank $seq.estimated_library_size } else { '' }
            sequencing_outcome_band = $sequencingBand
            alert = $row.alert
            observations = $row.observations
        }) | Out-Null
}

$scorecardRows = New-Object System.Collections.Generic.List[object]
foreach ($sampleId in ($sampleQualityMaster | Where-Object { $_.sample_id } | Select-Object -ExpandProperty sample_id -Unique | Sort-Object)) {
    $rows = @($sampleQualityMaster | Where-Object { $_.sample_id -eq $sampleId })
    $scorecardRows.Add([pscustomobject]@{
            sample_id = $sampleId
            batch_id = Get-DistinctJoinedValue -Rows $rows -Field 'batch_id'
            library_run_count = (@($rows)).Count
            stages = Get-DistinctJoinedValue -Rows $rows -Field 'stage'
            mean_reported_concentration = Get-MeanNumericValue -Rows $rows -Field 'reported_concentration'
            mean_region_average_size_bp = Get-MeanNumericValue -Rows $rows -Field 'region_average_size_bp'
            mean_region_percent_total = Get-MeanNumericValue -Rows $rows -Field 'region_percent_total'
            max_clean_dna_ng_ul = Get-MaxNumericValue -Rows $rows -Field 'clean_dna_ng_ul'
            max_qubit_ng_ul = Get-MaxNumericValue -Rows $rows -Field 'qubit_ng_ul'
            max_pre_capture_pcr_ng_ul = Get-MaxNumericValue -Rows $rows -Field 'pre_capture_pcr_ng_ul'
            mean_library_qc_heuristic_score = Get-MeanNumericValue -Rows $rows -Field 'library_qc_heuristic_score'
            sequencing_on_target_fraction = Get-MaxNumericValue -Rows $rows -Field 'sequencing_on_target_fraction'
            sequencing_target_mean_depth = Get-MaxNumericValue -Rows $rows -Field 'sequencing_target_mean_depth'
            sequencing_target_fold80 = Get-MaxNumericValue -Rows $rows -Field 'sequencing_target_fold80'
            sequencing_pct_duplication = Get-MaxNumericValue -Rows $rows -Field 'sequencing_pct_duplication'
            sequencing_outcome_band = Get-DistinctJoinedValue -Rows $rows -Field 'sequencing_outcome_band'
        }) | Out-Null
}

$analysisSummaryRows = New-Object System.Collections.Generic.List[object]
foreach ($group in @($scorecardRows | Group-Object sequencing_outcome_band | Sort-Object Name)) {
    $rows = @($group.Group)
    foreach ($metric in @('mean_region_average_size_bp', 'mean_region_percent_total', 'max_clean_dna_ng_ul', 'mean_library_qc_heuristic_score', 'sequencing_on_target_fraction', 'sequencing_target_mean_depth')) {
        $numericValues = @($rows | ForEach-Object {
                $value = Get-ObjectPropertyValue -Row $_ -Field $metric
                if ($value -ne '') { [double]$value }
            } | Sort-Object)
        $analysisSummaryRows.Add([pscustomobject]@{
                summary_type = 'scorecard_by_outcome'
                group_name = $group.Name
                metric = $metric
                n = $numericValues.Count
                mean = Get-MeanNumericValue -Rows $rows -Field $metric
                median = Get-MedianNumericValue -Rows $rows -Field $metric
                minimum = if ($numericValues.Count -gt 0) { $numericValues[0] } else { '' }
                maximum = if ($numericValues.Count -gt 0) { $numericValues[-1] } else { '' }
            }) | Out-Null
    }
}

foreach ($group in @($sampleQualityMaster | Group-Object stage | Sort-Object Name)) {
    $rows = @($group.Group)
    foreach ($metric in @('reported_concentration', 'region_average_size_bp', 'region_percent_total', 'library_qc_heuristic_score')) {
        $numericValues = @($rows | ForEach-Object {
                $value = Get-ObjectPropertyValue -Row $_ -Field $metric
                if ($value -ne '') { [double]$value }
            } | Sort-Object)
        $analysisSummaryRows.Add([pscustomobject]@{
                summary_type = 'library_rows_by_stage'
                group_name = $group.Name
                metric = $metric
                n = $numericValues.Count
                mean = Get-MeanNumericValue -Rows $rows -Field $metric
                median = Get-MedianNumericValue -Rows $rows -Field $metric
                minimum = if ($numericValues.Count -gt 0) { $numericValues[0] } else { '' }
                maximum = if ($numericValues.Count -gt 0) { $numericValues[-1] } else { '' }
            }) | Out-Null
    }
}

Write-Tsv -Path (Join-Path $OutputDir 'file_inventory.tsv') -Rows $fileInventory -Columns @('source_family', 'run_id', 'stage', 'instrument', 'file_type', 'parser_role', 'parse_strategy', 'path')
Write-Tsv -Path (Join-Path $OutputDir 'library_qc.tsv') -Rows $resolvedLibraryRows -Columns @('run_id', 'stage', 'instrument', 'parser_source', 'source_path', 'native_source_path', 'pdf_source_path', 'well', 'sample_description', 'sample_id', 'sample_id_source', 'batch_id', 'reported_concentration', 'reported_concentration_unit', 'peak_row_count', 'dominant_peak_size_bp', 'dominant_peak_calibrated_concentration', 'dominant_peak_molarity', 'dominant_peak_area_pct', 'region_from_bp', 'region_to_bp', 'region_average_size_bp', 'region_concentration', 'region_concentration_unit', 'region_molarity', 'region_molarity_unit', 'region_percent_total', 'alert', 'observations')
Write-Tsv -Path (Join-Path $OutputDir 'prep_metadata.tsv') -Rows $resolvedPrepMetadata -Columns @('source_family', 'source_run_id', 'source_path', 'sheet_name', 'row_index', 'raw_sample_label', 'sample_id', 'code', 'group', 'region', 'well_or_barcode', 'dna_input_ng_ul', 'clean_dna_ng_ul', 'qubit_ng_ul', 'pre_capture_pcr_ng_ul', 'molarity_nm', 'elution_volume_ul', 'notes')
Write-Tsv -Path (Join-Path $OutputDir 'input_dna_quantity.tsv') -Rows $inputDnaQuantity -Columns @('source_family', 'source_run_id', 'source_path', 'sheet_name', 'row_index', 'raw_sample_label', 'sample_id', 'code', 'group', 'region', 'well_or_barcode', 'dna_input_ng_ul', 'clean_dna_ng_ul', 'qubit_ng_ul', 'pre_capture_pcr_ng_ul', 'molarity_nm', 'elution_volume_ul', 'notes')
Write-Tsv -Path (Join-Path $OutputDir 'sample_aliases.tsv') -Rows $aliasEntries -Columns @('alias', 'alias_key', 'sample_id', 'source', 'source_run_id', 'notes')
Write-Tsv -Path (Join-Path $OutputDir 'sample_quality_master.tsv') -Rows $sampleQualityMaster -Columns @('sample_id', 'sample_id_source', 'batch_id', 'run_id', 'stage', 'instrument', 'parser_source', 'source_path', 'well', 'sample_description', 'reported_concentration', 'reported_concentration_unit', 'peak_row_count', 'dominant_peak_size_bp', 'dominant_peak_area_pct', 'region_average_size_bp', 'region_concentration', 'region_concentration_unit', 'region_molarity', 'region_percent_total', 'library_qc_heuristic_score', 'dna_input_ng_ul', 'clean_dna_ng_ul', 'qubit_ng_ul', 'pre_capture_pcr_ng_ul', 'prep_molarity_nm', 'prep_source_families', 'prep_regions', 'sequencing_metrics_path', 'sequencing_on_target_fraction', 'sequencing_target_mean_depth', 'sequencing_target_fold80', 'sequencing_pct_duplication', 'sequencing_estimated_library_size', 'sequencing_outcome_band', 'alert', 'observations')
Write-Tsv -Path (Join-Path $OutputDir 'dna_quality_scorecard.tsv') -Rows $scorecardRows -Columns @('sample_id', 'batch_id', 'library_run_count', 'stages', 'mean_reported_concentration', 'mean_region_average_size_bp', 'mean_region_percent_total', 'max_clean_dna_ng_ul', 'max_qubit_ng_ul', 'max_pre_capture_pcr_ng_ul', 'mean_library_qc_heuristic_score', 'sequencing_on_target_fraction', 'sequencing_target_mean_depth', 'sequencing_target_fold80', 'sequencing_pct_duplication', 'sequencing_outcome_band')
Write-Tsv -Path (Join-Path $OutputDir 'analysis_summary.tsv') -Rows $analysisSummaryRows -Columns @('summary_type', 'group_name', 'metric', 'n', 'mean', 'median', 'minimum', 'maximum')

$unresolvedLibraryCount = @($resolvedLibraryRows | Where-Object { -not $_.sample_id }).Count
$reportLines = @(
    '# DNA Quality Analysis',
    '',
    "Output run: $OutputRun",
    '',
    '## Scope',
    "- SureSelect files inventoried: $($sureSelectFiles.Count)",
    "- ddPCR files inventoried: $($ddpcrFiles.Count)",
    "- Samples metadata files inventoried: $($sampleFiles.Count)",
    "- Library QC rows parsed: $($resolvedLibraryRows.Count)",
    "- Prep metadata rows extracted: $($resolvedPrepMetadata.Count)",
    "- Library QC rows still unmatched to `sample_id`: $unresolvedLibraryCount",
    "- Samples with scorecards: $($scorecardRows.Count)",
    '',
    '## Notes',
    '- CSV sidecars were preferred over PDFs when both existed in the same Tapestation run directory.',
    '- Native `.D1000` / `.HSD1000` files are inventoried as provenance in `file_inventory.tsv` but are not parsed directly in this v1 workflow.',
    '- `library_qc_heuristic_score` is a transparent helper score, not a validated assay threshold.',
    '',
    '## Outputs',
    '- `file_inventory.tsv`: all useful external files discovered for this analysis.',
    '- `library_qc.tsv`: per-run Tapestation-derived QC rows.',
    '- `prep_metadata.tsv`: extracted supporting metadata from SureSelect, ddPCR, and Samples spreadsheets.',
    '- `input_dna_quantity.tsv`: subset of metadata rows with upstream quantity-related values.',
    '- `sample_aliases.tsv`: alias map used for sample-resolution joins.',
    '- `sample_quality_master.tsv`: joined library QC + upstream quantity + sequencing outcome table.',
    '- `dna_quality_scorecard.tsv`: per-sample summary for manuscript/plots.',
    '- `analysis_summary.tsv`: grouped descriptive summaries by sequencing outcome and QC stage.'
)
Set-Content -Path (Join-Path $OutputDir 'report.md') -Value $reportLines

Write-Host "DNA quality outputs written to $OutputDir"
