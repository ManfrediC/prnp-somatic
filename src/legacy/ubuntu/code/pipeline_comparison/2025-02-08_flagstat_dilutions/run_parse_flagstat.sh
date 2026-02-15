cd /add/results/pipeline_comparison/2025-02-08_flagstat_dilutions/
awk -f "/add/code/pipeline_comparison/2025-02-08_flagstat_dilutions/parse_flagstat.awk" *.flagstat.txt > summary_flagstat.tsv
