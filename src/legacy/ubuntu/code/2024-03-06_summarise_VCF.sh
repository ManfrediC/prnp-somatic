

cd /add/2024-03-04_force_mutect/VCFoutput/

# compress and index the VCF
for file in *vcf; do
  bgzip -c "$file" > "/add/2024-03-06_summarise_VCF/output/${file}.gz"
  tabix -p vcf "/add/2024-03-06_summarise_VCF/output/${file}.gz"
done

cd /add/2024-03-06_summarise_VCF/output/

# merge VCF -> file with single variant line
bcftools merge *.vcf.gz > merged_E200K.vcf


# format

{
  printf 'Sample\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tFORMAT-DP\tAD-REF\tAD-ALT\n';
  bcftools query --format "[%SAMPLE ]\t%CHROM\t%POS0\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT ]\t[%DP ]\t[%AD ] " merged_E200K.vcf |
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
} | column -s$'\t' -t > E200K_summary.tsv
