inputdir="/add/seq_data/2025-02-07_dilutions_bwa/"
outputdir="/add/results/pipeline_comparison/2025-02-08_flagstat_dilutions/"

samplename="NA100_undil"
sample="NA100_undil.bwa.picard.markedDup.bam"

mkdir -p "$outputdir"
samtools flagstat "${inputdir}${sample}" > "${outputdir}${samplename}.flagstat.txt"
