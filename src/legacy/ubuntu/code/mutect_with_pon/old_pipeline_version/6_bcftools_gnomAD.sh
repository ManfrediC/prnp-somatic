# 2025-04-19
# this step annotates the VCFs with the AF field from gnomAD, allowing to estimate if a variant is germline or not

for vcf in /add/results/annotation_output/funcotator_PoN/*.funcotated.vcf.gz; do
    base=$(basename "$vcf" .funcotated.vcf.gz)
    bcftools annotate \
        -a /home/mcarta/databases/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
        -c INFO/AF \
        -h <(echo '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="gnomAD_AF added post hoc">') \
        -Oz -o /add/results/annotation_output/funcotator_PoN/${base}.funcotated.af.vcf.gz \
        "$vcf"
    
    # index the output  
    tabix -f -p vcf /add/results/annotation_output/funcotator_PoN/${base}.funcotated.af.vcf.gz
done
