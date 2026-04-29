# Navigate to the directory where the VCF files are located
cd /add/results/2024-10-13_count_D178N/

# Compress and index the VCF files
for file in *vcf; do
  bgzip -c "$file" > "/add/results/2024-10-13_D178N_VCF_summary/${file}.gz"
  tabix -p vcf "/add/results/2024-10-13_D178N_VCF_summary/${file}.gz"
done

# Move to the output directory where compressed VCFs are stored
cd /add/results/2024-10-13_D178N_VCF_summary/

# Merge the compressed VCF files into a single merged VCF file
bcftools merge *.vcf.gz > merged_D178N.vcf

# Extract and format variant information
{
  printf 'Sample\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tFORMAT-DP\tAD-REF\tAD-ALT\n';
  bcftools query --format "[%SAMPLE ]\t%CHROM\t%POS0\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT ]\t[%DP ]\t[%AD ] " merged_D178N.vcf |
  awk -F'\t' -v OFS="\t" '
  { 
    split($1, samples, " "); 
    split($9, genotypes, " "); 
    split($10, dps, " ");
    split($11, ads, " "); 
    for(s in samples){ 
       split(ads[s], sampleAd, ",");  
       print samples[s], $2, $3, $4, $5, $6, $7, $8, genotypes[s], dps[s], sampleAd[1], sampleAd[2]
     }
  }';
} | column -s$'\t' -t > D178N_summary.tsv
