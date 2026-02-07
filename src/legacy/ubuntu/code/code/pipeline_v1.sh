#!/bin/bash
### Input directory parsed as variable
### $1 is the first positional parameter
inputDir=$1
cd $inputDir


#count number of rows that contain "fastq" = number of samples
nC=`cat dataset.tsv | grep fastq | wc -l`
echo "nC risulta $nC"

k=1 


#start loop; number of loops = number of samples (nC)
while [ $k -le $nC ]

do

# fetch names of samples to be processed from list in dataset.tsv
sampleNames=(dum $(cut -f 1 dataset.tsv | grep -v Name))
case1=(dum $(cut -f 2 dataset.tsv | grep -v Read1))
case2=(dum $(cut -f 3 dataset.tsv | grep -v Read2))
fileName=$(basename ${case1[$k]} _R1.fastq.gz)
sampleName=${sampleNames[$k]}


echo "the sample is $sampleName"


cores=4

### DATABASES DECLARATION

### Reference Fasta File, dict  and indexing 

###############regionFile - .bed is missing, but I don't think we need it

dbDir="/home/mcarta/databases"

fastaFile="$dbDir/chr2_chr4_chr20.fasta"

dbsnpDB="$dbDir/dbsnp_146.hg38.vcf.gz"
indelsDB="$dbDir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
regionFile="$dbDir/chr2_chr4_chr20.bed"

germlineResource="$dbDir/somatic-hg38_af-only-gnomad.hg38.vcf"
panelofNormals="$dbDir/somatic-hg38_1000g_pon.hg38.vcf.gz"


### PIPELINE STARTS


### Trimmomatic string declaration 

adapterFile="/home/mcarta/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"

trimmString="ILLUMINACLIP:$adapterFile:1:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:5:20 AVGQUAL:30  HEADCROP:0 MINLEN:80"
read1=${case1[$k]}
read2=${case2[$k]}

### Trimming/preprocessing: Trimmomatic
trimmomatic   PE -threads $cores -phred33 $read1 $read2  \
$fileName.R1.trimmed.fastq  $fileName.R1.unpaired \
$fileName.R2.trimmed.fastq $fileName.R2.unpaired $trimmString



### Mapping: Minimap2

##### -t indicates the number of threads. 8 is preset, but I've only got 6 CPUs -> change?

minimap2 -t 6 -a -x sr $fastaFile $fileName.R1.trimmed.fastq \
$fileName.R2.trimmed.fastq \
> $fileName.sam


### Fixing bam files: samtools
samtools view -hbS $fileName.sam |samtools sort - >  $fileName.bam
samtools index $fileName.bam


### Picard: add groups and mark duplicates  
picard AddOrReplaceReadGroups I=$fileName.bam \
O=$fileName.picard.bam SORT_ORDER=coordinate RGID=$fileName  \
RGLB=Paired_end RGPL=illumina RGSM=$fileName RGPU=illumina \
USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
     
picard MarkDuplicates I=$fileName.picard.bam \
O=$fileName.picard.markedDup.bam M=$fileName.picard.markedDup.metrics \
USE_JDK_DEFLATER=true USE_JDK_INFLATER=true


### GATK base recalibrator step 1
gatk --java-options -Xmx8g BaseRecalibrator \
   -I $fileName.picard.markedDup.bam \
   -R  $fastaFile \
   --known-sites $dbsnpDB \
   --known-sites $indelsDB \
   -O $fileName.recal_data.table
   
### GATK base recalibrator step 2
gatk --java-options -Xmx8g ApplyBQSR \
   -R $fastaFile \
   -I $fileName.picard.markedDup.bam \
   --bqsr-recal-file $fileName.recal_data.table \
   -O $fileName.picard.markedDup.recal.bam 

#move to next line
k=$(( $k + 1 ))

done

exit