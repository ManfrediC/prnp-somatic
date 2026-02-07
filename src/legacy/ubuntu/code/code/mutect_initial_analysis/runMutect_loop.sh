
# this command can be used to define the versions of java and gatk
#export PATH="/usr/lib/jvm/java-17-openjdk-amd64/bin:/home/mcarta/gatk-4.5.0.0:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/home/mcarta/miniconda3/bin"

#resources
dbDir="/home/mcarta/databases"

fastaFile="$dbDir/chr2_chr4_chr20.fasta"

dbsnpDB="$dbDir/dbsnp_146.hg38.vcf.gz"
indelsDB="$dbDir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
regionFile="$dbDir/chr2_chr4_chr20.bed"

germlineResource="$dbDir/somatic-hg38_af-only-gnomad.hg38.vcf"
panelofNormals="$dbDir/somatic-hg38_1000g_pon.hg38.vcf.vcf"


#loop
inputDir=$1
cd $inputDir

# to run, e.g. :
# bash runMutect_loop.sh "/add/2021-02-03_sequencing_of_dilutions/"
# bash runMutect_loop.sh "/add/2021-10-19_first_CJD_seq/"
# bash runMutect_loop.sh "/add/2023-06-18_CJD_16_samples/"
# bash runMutect_loop.sh "/add/2023-06-18_CJD_8_samples/"



#count number of rows that contain "recal" = number of samples to analyse

nC=`cat mutect_dataset.tsv | grep recal | wc -l`
echo "nC risulta $nC"

k=1 


#start loop; number of loops = number of samples (nC)
while [ $k -le $nC ]

do


sampleNames=(dum $(cut -f 1 mutect_dataset.tsv | grep -v Name))
echo "sampleNames risulta $sampleNames"

case1=(dum $(cut -f 2 mutect_dataset.tsv | grep -v InputFile))
echo "case1 risulta $case1"

fileName=$(basename ${case1[$k]} .picard.markedDup.recal.bam)
echo "fileName risulta $fileName"

sampleName=${sampleNames[$k]}
echo "sampleName risulta $sampleName"


#Mutect2

resultDir="/add/mutect_output/2021-10-19_first_CJD_seq/"


inputFile="$fileName.picard.markedDup.recal.bam"
echo "inputFile risulta $inputFile"

gatk Mutect2 \
   -R $fastaFile \
   -I $inputFile\
   --germline-resource $germlineResource \
   -O $resultDir/$sampleName.Mutect.vcf.gz
   

#move to next line
k=$(( $k + 1 ))

done

exit
