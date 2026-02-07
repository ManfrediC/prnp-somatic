#!/bin/bash

# Directory with filtered VCF files from dilution samples
VCF_DIR="/add/results/2025-04-01_filtered_prnp_vcf/dilutions"

# Output combined TSV table
OUTPUT_TSV="/add/results/2025-04-05_combined_variant_tables/2025-04-05_dilutions_variant_table.tsv"

# Temporary directory for intermediate files
TEMP_DIR="/tmp"

# Remove previous combined output if it exists
rm -f "$OUTPUT_TSV"

first=1

# Loop over each VCF file in the directory (assumed to end with _filtered_prnp.vcf)
for vcf in "$VCF_DIR"/*.vcf; do
  # Extract sample name from the filename (adjust the sed expression as needed)
  sample=$(basename "$vcf" | sed 's/_filtered_prnp.vcf//')
  echo "Processing sample: $sample from file: $vcf"
  
  # Run GATK VariantsToTable to convert the VCF to a temporary TSV file
  temp_tsv="${TEMP_DIR}/${sample}.tsv"
  gatk VariantsToTable \
  --variant "$vcf" \
  --output "$temp_tsv" \
  #--fields CHROM --fields POS --fields REF --fields ALT --fields FILTER \
  #--fields AS_FilterStatus --fields AS_SB_TABLE --fields DP --fields GERMQ --fields MBQ --fields MFRL --fields MMQ --fields MPOS --fields POPAF --fields TLOD \
  #--genotypeFields GT --genotypeFields AD --genotypeFields AF --genotypeFields DP --genotypeFields F1R2 --genotypeFields F2R1 --genotypeFields FAD

  
  # For the first file, write the header (with an added Sample column)
  if [ $first -eq 1 ]; then
      # Prepend "Sample" to the header row and write the complete file to the output
      sed '1s/^/Sample\t/' "$temp_tsv" > "$OUTPUT_TSV"
      # Append the data rows with the sample name inserted as the first column
      tail -n +2 "$temp_tsv" | awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUTPUT_TSV"
      first=0
  else
      # For subsequent files, skip the header and prepend the sample name to each row before appending
      tail -n +2 "$temp_tsv" | awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUTPUT_TSV"
  fi
done

echo "Combined variant table created at: $OUTPUT_TSV"