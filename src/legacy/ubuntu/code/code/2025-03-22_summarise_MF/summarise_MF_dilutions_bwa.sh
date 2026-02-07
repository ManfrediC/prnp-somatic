#!/bin/bash

# Create the output directory if needed
mkdir -p /add/results/2025-03-22_MF_dilution_summary/

output="/add/results/2025-03-22_MF_dilution_summary/MF_dilution_summary_bwa.tsv"

# Write header to TSV file
echo -e "sample\tchr\tposition\tmutation\tREF\tALT\tdepth\tREF_count\tALT_count\tconflict_count\tvariant_type\tVAF(%)" > "$output"

# Loop through each sample directory in the BWA dilution directory
for dir in /add/results/2025-02-10_MF_loop_bwa_dilutions/*_bwa; do
    sample=$(basename "$dir")
    file="${dir}/${sample}_all_candidates"

    # Check if the candidate file exists
    if [[ -f "$file" ]]; then
        while read -r line; do
            # Extract fields from the candidate file
            sample_col=$(echo "$line" | awk '{print $1}')
            chr=$(echo "$line" | awk '{print $2}')
            pos=$(echo "$line" | awk '{print $3}')
            ref=$(echo "$line" | awk '{print $4}')
            alt=$(echo "$line" | awk '{print $5}')
            ref_count=$(echo "$line" | awk '{print $6}')
            alt_count=$(echo "$line" | awk '{print $7}')
            conflict_count=$(echo "$line" | awk '{print $8}')
            variant_type=$(echo "$line" | awk '{print $9}')

            # Calculate depth = REF_reads + ALT_reads + conflict_count
            depth=$(echo "$ref_count + $alt_count + $conflict_count" | bc)

            # Calculate VAF(%) with high precision then round to 2 decimals
            if [[ "$depth" -gt 0 ]]; then
                raw_vaf=$(echo "scale=4; ($alt_count / $depth) * 100" | bc)
                vaf=$(printf "%.2f" "$raw_vaf")
            else
                vaf="NA"
            fi

            # Determine mutation name based on the position
            case $pos in
                4699818) mutation="E200K" ;;
                4699570) mutation="A117V" ;;
                4699752) mutation="D178N" ;;
                4699525) mutation="P102L" ;;
                *) mutation="Unknown" ;;
            esac

            # Write the output row with the desired columns
            echo -e "${sample_col}\t${chr}\t${pos}\t${mutation}\t${ref}\t${alt}\t${depth}\t${ref_count}\t${alt_count}\t${conflict_count}\t${variant_type}\t${vaf}" >> "$output"
        done < "$file"
    fi
done
