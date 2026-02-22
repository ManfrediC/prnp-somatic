#!/bin/bash
### Input directory parsed as variable
### $1 is the first positional parameter
inputDir=$1
cd $inputDir

echo $1
echo $inputDir

k=1 

sampleNames=(dum $(cut -f 1 dataset.tsv | grep -v Name))
case1=(dum $(cut -f 2 dataset.tsv | grep -v Read1))
case2=(dum $(cut -f 3 dataset.tsv | grep -v Read2))
fileName=$(basename ${case1[$k]} _R1.fastq.gz)
sampleName=${sampleNames[$k]}


cores=4

### Reference Fasta File, dict  and indexing 
#fastaName="chr2_chr4_chr20"
#fastaFile="$fastaName.fasta"

fastaFile="/home/mcarta/databases/chr2_chr4_chr20.fasta"

### DATABASES DECLARATION

###############regionFile - .bed is missing!!!!!

dbDir="/home/mcarta/databases/"
dbsnpDB="$dbDir/dbsnp_146.hg38.vcf.gz"
indelsDB="$dbDir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
regionFile="$dbDir/chr2_chr4_chr20.bed"


### PIPELINE STARTS

### Trimmomatic string declaration 
trimmString="ILLUMINACLIP:adapters.fa:1:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:5:20 AVGQUAL:30  HEADCROP:0 MINLEN:80"
read1=${case1[$k]}
read2=${case2[$k]}

### Trimming/preprocessing: Trimmomatic
trimmomatic   PE -threads $cores -phred33 $read1 $read2  \
$fileName.R1.trimmed.fastq  $fileName.R1.unpaired \
$fileName.R2.trimmed.fastq $fileName.R2.unpaired $trimmString

##no errors up to here


### Mapping: Minimap2

##### -t indicates the number of threads. 8 is preset, but I've only got 6 CPUs -> change?

minimap2 -t 6 -a -x sr $fastaFile $fileName.R1.trimmed.fastq \
$fileName.R2.trimmed.fastq \
> $fileName.sam


### Fixing bam files: samtools
conda activate samtools
samtools view -hbS $fileName.sam |samtools sort - >  $fileName.bam
samtools index $fileName.bam
conda deactivate