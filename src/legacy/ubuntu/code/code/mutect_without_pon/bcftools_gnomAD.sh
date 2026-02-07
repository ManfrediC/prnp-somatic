# 2025-04-19
# this script is a copy of the script designed for the PoN pipeline
# this step annotates the VCFs with the AF field from gnomAD, allowing to estimate if a variant is germline or not

for vcf in /add/results/no_PoN/annot/*.func.vcf.gz; do
    base=$(basename "$vcf" .func.vcf.gz)
    bcftools annotate \
        -a /home/mcarta/databases/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
        -c INFO/AF \
        -h <(echo '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="gnomAD_AF added post hoc">') \
        -Oz -o /add/results/no_PoN/annot_with_gnomAD/${base}.func.af.vcf.gz \
        "$vcf"
    
    # index the output  
    tabix -f -p vcf /add/results/no_PoN/annot_with_gnomAD/${base}.func.af.vcf.gz
done
