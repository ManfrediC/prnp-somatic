#!/usr/bin/awk -f
# parse_flagstat_consistent.awk
#
# This script parses samtools flagstat files and produces a summary table with one row per file.
# It extracts only the first numeric value from each matching metric line.
# The output columns are:
# Sample, Total, Primary, Secondary, Supplementary, Duplicates, PrimaryDuplicates,
# Mapped, PrimaryMapped, Paired, Read1, Read2, ProperlyPaired, WithItselfAndMateMapped,
# Singletons, WithMateDifferentChr, WithMateDifferentChrMapQ5
#
# Usage:
#   awk -f parse_flagstat_consistent.awk /add/results/pipeline_comparison/2025-02-08_flagstat_samples/*.txt > summary_consistent.tsv

BEGIN {
    OFS = "\t"
    print "Sample", "Total", "Primary", "Secondary", "Supplementary", "Duplicates", "PrimaryDuplicates", "Mapped", "PrimaryMapped", "Paired", "Read1", "Read2", "ProperlyPaired", "WithItselfAndMateMapped", "Singletons", "WithMateDifferentChr", "WithMateDifferentChrMapQ5"
}

# Helper function to extract only the first number from a line.
function getNumber(line) {
    if (match(line, /^([0-9]+)/, arr))
        return arr[1]
    return "NA"
}

# Reset metrics at the beginning of each file.
FNR == 1 {
    total = primary = secondary = supplementary = duplicates = primaryDuplicates = mapped = primaryMapped = paired = read1 = read2 = properlyPaired = withMateMapped = singletons = withMateDifferentChr = withMateDifferentChrMapQ5 = "NA"
}

# Extract the first number for each metric line.
$0 ~ / in total/         { total = getNumber($0) }
$0 ~ / primary$/ && $0 !~ /primary duplicates/ { primary = getNumber($0) }
$0 ~ / secondary/         { secondary = getNumber($0) }
$0 ~ / supplementary/     { supplementary = getNumber($0) }
$0 ~ /duplicates/ && $0 ~ /primary duplicates/ { primaryDuplicates = getNumber($0) }
$0 ~ /duplicates/ && $0 !~ /primary duplicates/ { duplicates = getNumber($0) }
$0 ~ / mapped/ && $0 ~ /\(.*%/ && $0 !~ /primary mapped/ { mapped = getNumber($0) }
$0 ~ / primary mapped/   { primaryMapped = getNumber($0) }
$0 ~ / paired in sequencing/ { paired = getNumber($0) }
$0 ~ / read1/            { read1 = getNumber($0) }
$0 ~ / read2/            { read2 = getNumber($0) }
$0 ~ / properly paired/  { properlyPaired = getNumber($0) }
$0 ~ / with itself and mate mapped/ { withMateMapped = getNumber($0) }
$0 ~ / singletons/       { singletons = getNumber($0) }
$0 ~ / with mate mapped to a different chr/ && $0 !~ /mapQ>=5/ { withMateDifferentChr = getNumber($0) }
$0 ~ / with mate mapped to a different chr/ && $0 ~ /mapQ>=5/ { withMateDifferentChrMapQ5 = getNumber($0) }

ENDFILE {
    # Extract sample name from FILENAME:
    sample = FILENAME
    sub(".*/", "", sample)
    sub(/\.flagstat(\.txt)?$/, "", sample)
    print sample, total, primary, secondary, supplementary, duplicates, primaryDuplicates, mapped, primaryMapped, paired, read1, read2, properlyPaired, withMateMapped, singletons, withMateDifferentChr, withMateDifferentChrMapQ5
}
