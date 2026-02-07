
inputdir="/add/seq_data/2025-02-07_dilutions_bwa/"
outputdir="/add/results/pipeline_comparison/2025-02-08_flagstat_dilutions/"

samplename="NA100_undil"
sample="20210203.0-o23928_1_1-NA100_undil_D01.picard.markedDup.recal.bam"

samtools flagstat "$inputdir$sample" > "${outputdir}${samplename}.flagstat.txt"
