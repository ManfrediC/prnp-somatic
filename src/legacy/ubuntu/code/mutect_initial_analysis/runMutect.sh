
#export PATH="/usr/lib/jvm/java-17-openjdk-amd64/bin:/home/mcarta/gatk-4.5.0.0:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/home/mcarta/miniconda3/bin"



dbDir="/home/mcarta/databases"

fastaFile="$dbDir/chr2_chr4_chr20.fasta"

dbsnpDB="$dbDir/dbsnp_146.hg38.vcf.gz"
indelsDB="$dbDir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
regionFile="$dbDir/chr2_chr4_chr20.bed"

germlineResource="$dbDir/somatic-hg38_af-only-gnomad.hg38.vcf"
panelofNormals="$dbDir/somatic-hg38_1000g_pon.hg38.vcf.vcf"



#now: 1 file at a time. Edit inputFile

inputFile="/add/2021-10-19_first_CJD_seq/20211027.0-o26424_1_3-CJD6_fr.picard.markedDup.recal.bam"



resultDir="/add/mutect_output/2021-10-19_first_CJD_seq/"

gatk Mutect2 \
   -R $fastaFile \
   -I $inputFile \
   --germline-resource $germlineResource \
   -O $resultDir/testoutput.vcf.gz
   
exit
