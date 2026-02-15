cd /add/results/pipeline_comparison/2025-02-08_flagstat_samples/
awk -f /add/code/pipeline_comparison/2025-02-08_flagstat_samples/parse_flagstat_all.awk *.flagstat.txt > summary_flagstat_all.tsv