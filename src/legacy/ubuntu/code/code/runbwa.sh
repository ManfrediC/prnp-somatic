#the FASTA has already been indexed: /mcarta/databases/chr2_chr4_chr20.fasta.fai

#same input as minimap2





#alignment with bwa mem
#https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/alignment_tools/bwa/running_bwa_commands/
bwa mem index_prefix [input_reads.fastq|input_reads_pair_1.fastq input_reads_pair_2.fastq] [options]


minimap2 -t 6 -a -x sr $fastaFile $fileName.R1.trimmed.fastq \
$fileName.R2.trimmed.fastq \
> $fileName.sam
