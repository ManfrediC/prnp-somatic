

#bash /add/runFilterMutectCalls_loop.sh /add/mutect_output/2021-10-19_first_CJD_seq/

#bash /add/runFilterMutectCalls_loop.sh "/add/mutect_output/2021-02-03_sequencing_of_dilutions/"
#bash /add/runFilterMutectCalls_loop.sh "/add/mutect_output/2023-06-18_CJD_16_samples/"
#bash /add/runFilterMutectCalls_loop.sh "/add/mutect_output/2023-06-18_CJD_8_samples/"


#databases

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


#dataset.tsv with files to be processed
datasetFile="filter_dataset.tsv"

#count number of rows that contain "recal" = number of samples to analyse
nC=`cat $datasetFile | grep Mutect | wc -l`
echo "nC risulta $nC"

k=1 


#start loop; number of loops = number of samples (nC)
while [ $k -le $nC ]

do


sampleNames=(dum $(cut -f 1 $datasetFile | grep -v Name))
echo "sampleNames risulta $sampleNames"

case1=(dum $(cut -f 2 $datasetFile | grep -v InputFile))
echo "case1 risulta $case1"

fileName=$(basename ${case1[$k]} .vcf.gz)
echo "fileName risulta $fileName"

sampleName=${sampleNames[$k]}
echo "sampleName risulta $sampleName"


#input files for FilterMutectCalls

inputFile="$fileName.vcf.gz"
echo "inputFile risulta $inputFile"

outputDir="/add/mutect_filtered/"

myname=$(basename $inputVariants .vcf.gz)
echo "myname risulta $myname"


#filter F1

gatk FilterMutectCalls \
   -R $fastaFile \
   -V $inputDir/$inputFile \
   --f-score-beta 1 \
   -O $outputDir/$sampleName.filteredF1.vcf.gz

#filter F2

gatk FilterMutectCalls \
   -R $fastaFile \
   -V $inputDir/$inputFile \
   --f-score-beta 2 \
   -O $outputDir/$sampleName.filteredF2.vcf.gz

#filter F5

gatk FilterMutectCalls \
   -R $fastaFile \
   -V $inputDir/$inputFile \
   --f-score-beta 5 \
   -O $outputDir/$sampleName.filteredF5.vcf.gz

#move to next line
k=$(( $k + 1 ))

done

exit