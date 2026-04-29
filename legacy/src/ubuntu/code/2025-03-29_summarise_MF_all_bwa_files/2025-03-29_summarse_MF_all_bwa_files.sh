#!/bin/bash
# This script summarizes the other MosaicForecast output files 
# (all_2x2table, all.phasing, all.phasing_2by2, all.merged.inforSNPs.pos, and multiple_inforSNPs.log)
# from sample directories under /add/results/2025-02-10_MF_loop_bwa/
# and writes the summary files to /add/results/2025-03-29_MF_files_summary_bwa_samples/

# Create the output directory if it doesn't exist
mkdir -p /add/results/2025-03-29_MF_files_summary_bwa_samples/

#########################################
# 1. Summarise all_2x2table files
#########################################
output_2x2="/add/results/2025-03-29_MF_files_summary_bwa_samples/MF_files_summary_2x2table.tsv"
echo -e "sample\tchr\tpos\tREF\tALT\tinforSNP_pos\tinforSNP_ref\tinforSNP_alt\tconflict_count\tmajor_het1\tmajor_het2\tminor_het1\tminor_het2\tvariant_type" > "$output_2x2"

for dir in /add/results/2025-02-10_MF_loop_bwa/*_bwa; do
    sample=$(basename "$dir")
    file="${dir}/${sample}_all_2x2table"
    if [[ -f "$file" ]]; then
        # Append file contents as-is (assuming space-delimited)
        cat "$file" >> "$output_2x2"
    fi
done

#########################################
# 2. Summarise all.phasing files
#########################################
output_phasing="/add/results/2025-03-29_MF_files_summary_bwa_samples/MF_files_summary_phasing.tsv"
echo -e "sample\tchr\tpos\tREF\tALT\tphasing\tconflicting_reads\tmappability\tvariant_type" > "$output_phasing"

for dir in /add/results/2025-02-10_MF_loop_bwa/*_bwa; do
    sample=$(basename "$dir")
    file="${dir}/${sample}_all.phasing"
    if [[ -f "$file" ]]; then
        while read -r line; do
            # Skip header if present
            if [[ "$line" == sample* ]]; then
                continue
            fi
            echo -e "$line" >> "$output_phasing"
        done < "$file"
    fi
done

#########################################
# 3. Summarise all.phasing_2by2 files
#########################################
output_phasing2="/add/results/2025-03-29_MF_files_summary_bwa_samples/MF_files_summary_phasing2by2.tsv"
echo -e "sample\tchr\tpos\tREF\tALT\tphasing_info" > "$output_phasing2"

for dir in /add/results/2025-02-10_MF_loop_bwa/*_bwa; do
    sample=$(basename "$dir")
    file="${dir}/${sample}_all.phasing_2by2"
    if [[ -f "$file" ]]; then
        while read -r line; do
            # Assume the first token is a key in the format sample;chr;pos;REF;ALT, followed by phasing info
            key=$(echo "$line" | awk '{print $1}')
            phasing_info=$(echo "$line" | cut -d' ' -f2-)
            sample_field=$(echo "$key" | cut -d';' -f1)
            chr_field=$(echo "$key" | cut -d';' -f2)
            pos_field=$(echo "$key" | cut -d';' -f3)
            ref_field=$(echo "$key" | cut -d';' -f4)
            alt_field=$(echo "$key" | cut -d';' -f5)
            echo -e "${sample_field}\t${chr_field}\t${pos_field}\t${ref_field}\t${alt_field}\t${phasing_info}" >> "$output_phasing2"
        done < "$file"
    fi
done

#########################################
# 4. Summarise all.merged.inforSNPs.pos files
#########################################
output_merged="/add/results/2025-03-29_MF_files_summary_bwa_samples/MF_files_summary_merged_inforSNPs.tsv"
echo -e "sample\tchr\tpos\tREF\tALT\tinforSNP_pos\tinforSNP_ref\tinforSNP_alt\tconflict_count\tvariant_type" > "$output_merged"

for dir in /add/results/2025-02-10_MF_loop_bwa/*_bwa; do
    sample=$(basename "$dir")
    file="${dir}/${sample}_all.merged.inforSNPs.pos"
    if [[ -f "$file" ]]; then
        while read -r line; do
            # Expected fields:
            # sample, chr, pos, REF, ALT, duplicate_chr, inforSNP_pos, inforSNP_ref, inforSNP_alt, conflict_count, variant_type
            sample_field=$(echo "$line" | awk '{print $1}')
            chr_field=$(echo "$line" | awk '{print $2}')
            pos_field=$(echo "$line" | awk '{print $3}')
            ref_field=$(echo "$line" | awk '{print $4}')
            alt_field=$(echo "$line" | awk '{print $5}')
            inforSNP_pos=$(echo "$line" | awk '{print $7}')
            inforSNP_ref=$(echo "$line" | awk '{print $8}')
            inforSNP_alt=$(echo "$line" | awk '{print $9}')
            conflict_count=$(echo "$line" | awk '{print $10}')
            variant_type=$(echo "$line" | awk '{print $11}')
            echo -e "${sample_field}\t${chr_field}\t${pos_field}\t${ref_field}\t${alt_field}\t${inforSNP_pos}\t${inforSNP_ref}\t${inforSNP_alt}\t${conflict_count}\t${variant_type}" >> "$output_merged"
        done < "$file"
    fi
done

#########################################
# 5. Summarise multiple_inforSNPs.log files
#########################################
output_log="/add/results/2025-03-29_MF_files_summary_bwa_samples/MF_files_summary_inforSNPs_log.tsv"
echo -e "sample\tchr\tpos\tREF\tALT\tlog_details" > "$output_log"

for dir in /add/results/2025-02-10_MF_loop_bwa/*_bwa; do
    sample=$(basename "$dir")
    file="${dir}/${sample}_multiple_inforSNPs.log"
    if [[ -f "$file" ]]; then
        while read -r line; do
            # Expect a key in the format sample;chr;pos;REF;ALT followed by log details
            key=$(echo "$line" | awk '{print $1}')
            log_details=$(echo "$line" | cut -d' ' -f2-)
            sample_field=$(echo "$key" | cut -d';' -f1)
            chr_field=$(echo "$key" | cut -d';' -f2)
            pos_field=$(echo "$key" | cut -d';' -f3)
            ref_field=$(echo "$key" | cut -d';' -f4)
            alt_field=$(echo "$key" | cut -d';' -f5)
            echo -e "${sample_field}\t${chr_field}\t${pos_field}\t${ref_field}\t${alt_field}\t${log_details}" >> "$output_log"
        done < "$file"
    fi
done

echo "Summary files created in /add/results/2025-03-29_MF_files_summary_bwa_samples/"
