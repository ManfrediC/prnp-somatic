

#databases

dbDir="/home/mcarta/databases"

fastaFile="$dbDir/chr2_chr4_chr20.fasta"
dbsnpDB="$dbDir/dbsnp_146.hg38.vcf.gz"
indelsDB="$dbDir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
regionFile="$dbDir/chr2_chr4_chr20.bed"
germlineResource="$dbDir/somatic-hg38_af-only-gnomad.hg38.vcf"
panelofNormals="$dbDir/somatic-hg38_1000g_pon.hg38.vcf.vcf"


#input files for FilterMutectCalls

inputVariants="312092_01-CJD5_A01.Mutect.vcf.gz"

inputDir="/add/mutect_output/2023-06-18_CJD_16_samples/"
outputDir="/add/mutect_filtered/"

myname=$(basename $inputVariants .vcf.gz)
echo "myname risulta $myname"


#filter

gatk FilterMutectCalls \
   -R $fastaFile \
   -V $inputDir/$inputVariants \
   --f-score-beta 5 \
   -O $outputDir/$myname.filteredF5.vcf.gz
