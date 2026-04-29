These are the scripts to force Mutect2 to call specific variants, even if they aren't present. 

samples.tsv is a list of all the CJD/Ctrl samples with the correct .bam to process.

Current status:
- E200K: works
- D178N: doesn't work, output VCFs are empty (perhaps because there are no variants?? But it works on samples without variants with E200K)
- M129V: 