#!/usr/bin/awk -f
# parse_flagstat.awk
#
# This script parses a samtools flagstat file and produces a summary table.
# For each file it extracts the following metrics:
#   Sample, Total, Primary, Secondary, Supplementary, Duplicates, PrimaryDuplicates,
#   Mapped, PrimaryMapped, Paired, Read1, Read2, ProperlyPaired,
#   WithItselfAndMateMapped, Singletons, WithMateDifferentChr, WithMateDifferentChrMapQ5
#
# Usage:
#   awk -f parse_flagstat.awk *.flagstat.txt
#
# The script extracts the sample name from the filename by removing the path and the suffix.
#
# NOTE: This script uses the ENDFILE block available in GNU awk 4.1+.

BEGIN {
    # Print header
    print "Sample\tTotal\tPrimary\tSecondary\tSupplementary\tDuplicates\tPrimaryDuplicates\tMapped\tPrimaryMapped\tPaired\tRead1\tRead2\tProperlyPaired\tWithItselfAndMateMapped\tSingletons\tWithMateDifferentChr\tWithMateDifferentChrMapQ5"
}

# Reset variables at the start of each file.
FNR==1 {
    total = primary = secondary = supplementary = duplicates = primaryDuplicates = mapped = primaryMapped = paired = read1 = read2 = properlyPaired = withMateMapped = singletons = withMateDifferentChr = withMateDifferentChrMapQ5 = ""
}

# Lines are assumed to be in the format:
#   <count> + <fail> <text...>
# We capture the first field ($1) from lines matching our keywords.

$0 ~ / in total/         { total = $1 }
$0 ~ / primary$/ && $0 !~ /primary duplicates/ { primary = $1 }
$0 ~ / secondary/         { secondary = $1 }
$0 ~ / supplementary/     { supplementary = $1 }
$0 ~ /duplicates/ && $0 ~ /primary duplicates/ { primaryDuplicates = $1 }
$0 ~ /duplicates/ && $0 !~ /primary duplicates/ { duplicates = $1 }
$0 ~ / mapped/ && $0 ~ /\(.*%/ && $0 !~ /primary mapped/ { mapped = $1 }
$0 ~ / primary mapped/   { primaryMapped = $1 }
$0 ~ / paired in sequencing/ { paired = $1 }
$0 ~ / read1/            { read1 = $1 }
$0 ~ / read2/            { read2 = $1 }
$0 ~ / properly paired/  { properlyPaired = $1 }
$0 ~ / with itself and mate mapped/ { withMateMapped = $1 }
$0 ~ / singletons/       { singletons = $1 }
$0 ~ / with mate mapped to a different chr/ && $0 !~ /mapQ>=5/ { withMateDifferentChr = $1 }
$0 ~ / with mate mapped to a different chr/ && $0 ~ /mapQ>=5/ { withMateDifferentChrMapQ5 = $1 }

ENDFILE {
    # Extract sample name from FILENAME: remove any path and remove the trailing .flagstat or .flagstat.txt suffix.
    sample = FILENAME
    sub(".*/", "", sample)
    sub(/\.flagstat(\.txt)?$/, "", sample)
    # Print a row with all columns
    print sample "\t" total "\t" primary "\t" secondary "\t" supplementary "\t" duplicates "\t" primaryDuplicates "\t" mapped "\t" primaryMapped "\t" paired "\t" read1 "\t" read2 "\t" properlyPaired "\t" withMateMapped "\t" singletons "\t" withMateDifferentChr "\t" withMateDifferentChrMapQ5
}
