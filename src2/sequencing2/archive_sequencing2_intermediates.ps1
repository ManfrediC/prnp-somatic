[CmdletBinding()]
param(
    [string]$RepositoryRoot = 'C:\Projects\prnp-somatic',
    [string]$ArchiveRoot = 'C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD intermediates',
    [string]$RunId = '',
    [switch]$RemoveLocal,
    [switch]$MakeOnlineOnly
)

$ErrorActionPreference = 'Stop'

function Get-ArchiveFiles {
    param([string]$Root)

    Get-ChildItem -LiteralPath $Root -File -Recurse | Where-Object {
        $_.FullName -notmatch '[\\/]bam_work[\\/]'
    }
}

$repository = (Resolve-Path -LiteralPath $RepositoryRoot).Path
$source = Join-Path $repository 'results2\sequencing2\runs'
$results = Join-Path $repository 'results2\sequencing2\results'
$archiveParent = (Resolve-Path -LiteralPath $ArchiveRoot).Path

if (-not (Test-Path -LiteralPath $source -PathType Container)) {
    throw "Intermediate directory does not exist: $source"
}

$requiredResults = @(
    'mutect2_controls_no_pon\variant_qc\filtered_variants.tsv',
    'mutect2_controls_pon\panel_of_normals\CJD_controls_PoN.vcf.gz',
    'mutect2_cjd_dilutions_with_pon\variant_qc\cjd\filtered_variants.tsv',
    'mutect2_cjd_dilutions_with_pon\variant_qc\dilutions\filtered_variants.tsv'
)
foreach ($relativePath in $requiredResults) {
    $resultPath = Join-Path $results $relativePath
    if (-not (Test-Path -LiteralPath $resultPath -PathType Leaf) -or
        (Get-Item -LiteralPath $resultPath).Length -eq 0) {
        throw "Required final output is missing or empty: $resultPath"
    }
}

$commit = (& git -C $repository rev-parse --short HEAD).Trim()
if ([string]::IsNullOrWhiteSpace($RunId)) {
    $RunId = 'sequencing2_{0}_{1}' -f (Get-Date -Format 'yyyy-MM-dd_HHmmss'), $commit
}
if ($RunId -notmatch '^[A-Za-z0-9._-]+$') {
    throw 'RunId may contain only letters, numbers, dots, underscores and hyphens.'
}

$destination = Join-Path $archiveParent $RunId
if (Test-Path -LiteralPath $destination) {
    throw "Archive destination already exists: $destination"
}
New-Item -ItemType Directory -Path $destination | Out-Null

$sourceFiles = @(Get-ArchiveFiles -Root $source)
if ($sourceFiles.Count -eq 0) {
    throw "No intermediate files found under: $source"
}

$manifestRows = foreach ($file in $sourceFiles) {
    [pscustomobject]@{
        Path = $file.FullName.Substring($source.Length).TrimStart('\', '/')
        Bytes = $file.Length
        SHA256 = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash
    }
}

$excludedDirectories = @(
    Get-ChildItem -LiteralPath $source -Directory -Recurse |
        Where-Object Name -eq 'bam_work' |
        ForEach-Object FullName
)
$robocopyArguments = @($source, $destination, '/E', '/COPY:DAT', '/DCOPY:DAT', '/R:2', '/W:2', '/XJ')
if ($excludedDirectories.Count -gt 0) {
    $robocopyArguments += '/XD'
    $robocopyArguments += $excludedDirectories
}

& robocopy @robocopyArguments
if ($LASTEXITCODE -gt 7) {
    throw "Robocopy failed with exit code $LASTEXITCODE"
}

$destinationFiles = @(Get-ArchiveFiles -Root $destination)
if ($destinationFiles.Count -ne $manifestRows.Count) {
    throw "Archive file count mismatch: source $($manifestRows.Count), destination $($destinationFiles.Count)"
}
foreach ($row in $manifestRows) {
    $copiedPath = Join-Path $destination $row.Path
    if (-not (Test-Path -LiteralPath $copiedPath -PathType Leaf)) {
        throw "Archived file is missing: $copiedPath"
    }
    $copiedFile = Get-Item -LiteralPath $copiedPath
    $copiedHash = (Get-FileHash -LiteralPath $copiedPath -Algorithm SHA256).Hash
    if ($copiedFile.Length -ne $row.Bytes -or $copiedHash -ne $row.SHA256) {
        throw "Archived file failed verification: $copiedPath"
    }
}

$manifestPath = Join-Path $destination 'MANIFEST.tsv'
$manifestRows | Export-Csv -LiteralPath $manifestPath -Delimiter "`t" -NoTypeInformation

$status = (& git -C $repository status --porcelain) -join '; '
$readme = @(
    'Sequencing2 intermediate archive',
    "Created: $((Get-Date).ToString('o'))",
    "Source: $source",
    "Git commit: $commit",
    "Git status: $status",
    'Excluded: readcount_qc/bam_work directories (regenerable links to original BAMs)',
    'Retained locally: results/final_bam and results2/sequencing2/results'
)
Set-Content -LiteralPath (Join-Path $destination 'README.txt') -Value $readme -Encoding utf8

if ($RemoveLocal) {
    Remove-Item -LiteralPath $source -Recurse -Force
}
if ($MakeOnlineOnly) {
    & attrib.exe +U -P $destination
    & attrib.exe +U -P (Join-Path $destination '*') /S /D
    if ($LASTEXITCODE -ne 0) {
        throw "Could not mark archive online-only: $destination"
    }
}

Write-Host "Verified archive: $destination"
Write-Host "Files archived: $($manifestRows.Count)"
Write-Host "Local intermediates removed: $($RemoveLocal.IsPresent)"
Write-Host "Online-only requested: $($MakeOnlineOnly.IsPresent)"
