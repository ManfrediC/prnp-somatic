## 03.11.2023

I’ve checked the pipeline (pipeline.sh). The latest version of this file will always be saved on Windows. Following things still need to be done:
- Varscan: I don’t think the parameters are correct. We need parameters that call somatic mutations without a matched sample with a reasonable threshold (probably VAF 0.05% - 0.5% ; that’s what I wrote in the thesis)
- GATK germline caller: check parameters. Maybe we don’t need this.
- Mutect doesn’t work yet.
- When these things have been fixed: create loop to analyse multiple files.

## 20.11.2023

- I’m expanding the pipeline by adding code to test.sh step by step.

- Mutect
- I need tumor-only mode, as I haven’t got normal reference samples.
- This mode seems to require a panel of normals (PoN) and a germline resource.
- It might be good to use the multi-sample mode in a second step (https://gatk.broadinstitute.org/hc/en-us/articles/360050722212-FAQ-for-Mutect2)

Sample from tutorial (https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
  gatk Mutect2 \
  -R reference.fa \
  -I sample.bam \
  --germline-resource af-only-gnomad.vcf.gz \
  --panel-of-normals pon.vcf.gz \
  -O single_sample.vcf.gz

- Germline-resource: I downloaded somatic-hg38_af-only-gnomad.hg38.vcf.vcf (Google Cloud Best Practices folder.)
- Panel of normals:  somatic-hg38_1000g_pon.hg38.vcf
- Why are they called “.vcf.vcf” ??
- Added shortcuts to these two files in the list at the top of the pipeline ($germlineResource, $panelofNormals)
- To save time I’m creating test-mutect.sh without the prior pipeline steps.

Mutect2 gives me an error: tumor samples are required (“—tumor-sample”” under required arguments). According to the tutorial I don’t need this  unsolved issue
How does one detect somatic mutations without matched tissues? See papers in “useful papers”
- Shi et al: unmatched analyses produce many technical articfacts (34-80%)
- Huang and Lee: list of tools for somatic mutations calling.
- Optimised for bulk DNA, non-cancer, without matched control: MosaicHunter 2017, EM-Mosaic 2020, MosaicForecast 2020, RePlow 2019
- MuTect should also work, but not optimised for this purpose.
- Ha 2023 (benchmark for mosaic detection): 
- SNV: MosaicForecast and Mutect2 are best for <5% VAF in single samples (fig 3f and Fig 5).
- Indel: MosaicForecast ist best, but no algorithms works well for <4%
- Finding overlaps from different callers to achieve high confidence is not recommended. use of multiple callers should be composed in a way of assigning one to the best-performing VAF area.

Created loop for GATK basic steps: pipeline.sh.
Ran basic steps on all raw files except folder CJD_16 samples (not enough disc space).

## 21.11.23

- IGV didn’t work (wrong version of Java): Installed Java 17 and changed PATH (as it was directed towards the old Java version in the Miniconda folder)  changing PATH solved problem.
- I viewed a few results: PRNP has high coverage (build has to be GRCh38/hg38)
- Install MosaicForecast (https://github.com/parklab/MosaicForecast/blob/master/FAQ.md)

### Generate variant call set as input of Mosaic Forecast with Mutect2

https://github.com/parklab/MosaicForecast/blob/master/FAQ.md
- Step1: Mutect2 in tumor-only mode for each normal sample
- It would make sense to provide an interval list that restricts the analysis to the regions targeted by sureselect (first I’ll try without, make this later)
- https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists

```
gatk Mutect2 \
   -R ref_fasta.fa \
   -I normal1.bam \
   -tumor normal1_sample_name \
   --germline-resource af-only-gnomad.vcf.gz \
   -L intervals.list \
   --interval-padding 100 \
   -O normal1_for_pon.vcf.gz

As of v4.1, there is no longer a need to specify the tumor sample name with -tumor. You need only specify the normal sample name with -normal, if you include a normal.

I seem to have version 4.0.5.1, so maybe this is why it wants -tumor


```

To do:
PATH to gatk might be wrong. It goes to miniconda (gatk4 4.0.5.1; using which gatk, whereas it should go to /home/mcarta/gatk
This might explain why mutect doesn’t work in tumor only mode (old version of mutect)

## 22.11.23

- Fixed path to use Java 17 and Gatk 4.4 (PATH directed it towards old versions in the miniconda3 directory
- This didn’t resolve the Java problem, so I changed “java” in Miniconda to “java-old”  solved
- Created /add disc with 200GB for extra space

## 23.11.23

- Trimmomatic gives an error: “adapters.fa not found”
- Installed fastqc

![Embedded image](analysis_journal_media/image1.png)

- In fact, the “trimmed” file still contains the same amount of Illumina adapters as the raw data  trimming didn’t work, probably because ILLUMINACLIP isn’t defined correctly.
- https://bioinformatics-core-shared-training.github.io/Bulk_RNAseq_Course_Nov22/Bulk_RNAseq_Course_Base/Markdowns/S3_Trimming_Reads.html

Old:

### Trimmomatic string declaration 
trimmString="ILLUMINACLIP:adapters.fa:1:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:5:20 AVGQUAL:30  HEADCROP:0 MINLEN:80"
read1=${case1[$k]}
read2=${case2[$k]}

### Trimming/preprocessing: Trimmomatic
trimmomatic   PE -threads $cores -phred33 $read1 $read2  \
$fileName.R1.trimmed.fastq  $fileName.R1.unpaired \
$fileName.R2.trimmed.fastq $fileName.R2.unpaired $trimmString

On the website they use the following to remove the Illumina universal adapter:

```
ILLUMINACLIP:trimmomatic/adapters/TruSeq3-SE.fa:2:30:7
```

New, with directory of adapter file for Paired Ends (PE), saved as testTrim.sh:
### Trimmomatic string declaration 
adapterFile="/home/mcarta/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"
trimmString="ILLUMINACLIP:$adapterFile:1:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:5:20 AVGQUAL:30  HEADCROP:0 MINLEN:80"
read1=${case1[$k]}
read2=${case2[$k]}

### Trimming/preprocessing: Trimmomatic
trimmomatic   PE -threads $cores -phred33 $read1 $read2  \
$fileName.R1.trimmed.fastq  $fileName.R1.unpaired \
$fileName.R2.trimmed.fastq $fileName.R2.unpaired $trimmString
- I just tried to run the entire pipeline on the correctly trimmed data. Picard didn’t work, as it seems to use the version of java in miniconda, whose name I changed to java-old  changed back.
- Now both Picard and IGV work

Mutect2
- Mutect with basic settings works now.
- To do: look up advanced settings with PoN and References, create PoN for MosaicF

## 07.01.2024

Mutect2
- Created loop: /add/runMutect_loop.sh
- The files need to be specified in mutect_dataset.tsv, which is inside the folder with the sequencing files
- ERROR: Unable to trim uncertain bases without flow order information
- https://gatk.broadinstitute.org/hc/en-us/community/posts/21837856660251-Mutect2-Unable-to-trim-uncertain-bases-without-flow-order-information
- It might be possible to resolve this by using bwa instead of minimap, but wait for an answer first (https://gatk.broadinstitute.org/hc/en-us/community/posts/11440622639387-Unable-to-trim-uncertain-bases-without-flow-order-information?page=1#community_comment_12925222020763)
- To do: once this is resolved, run with basic settings, then look up advanced settings with PoN and References, create PoN for MosaicF

BWA
- /add/runbwa.sh
- I don’t need this after all, as gatk 4.5 does not contain the issue with Minimap2

## 09.01.2024

- The error above (07.01.2024) is due to a bug in gatk 4.4, which is resolved in gatk 4.5
- I downloaded gatk 4.5 to /home/mcarta
- Modified /mcarta/bashrc and /etc/environment to correct path (/mcarta/gatk-4.5.0.0)
- However, Mutect2 still invokes v4.4
- We resolved this by renaming /gatk-4.5.0.0 to /gatk-4.4.0.0
- I also renamed the old version to gatk-4.4.0.0_old
- Created soft link
- ln -s gatk-4.4.0.0 gatk-4.5.0.0

### Mutect2

- I’ve analysed the dilution / spike-in set of samples
- In NA99A1_undil_D02, the A117V spike-in is detected (red arrow):

![Embedded image](analysis_journal_media/image2.png)

- Spike-in also detected in NA995A05_undil: allele fraction is also about 1% (rather than 0.5%)
- Todo: it might be good to make a list of the Mutect2-detected variants in the spike in and compare them with the “true” data (sequences of NA DNA and spike-in) to see how reliable Mutect is. If reliable, the unclear variants below might be less likely to be real.

![Embedded image](analysis_journal_media/image3.png)

Image: two non-confirmed variants in CJD3
I’ve inspected the Mutect results. CJD30 was heterozygous for E200K. Occasionally there were 1% variants (not flagged by Mutect) that were seen in several samples, probably artifacts. Mutect flagged some variants with extremely low prevalence in CJD4, but they probably aren’t pathogenic (most of them non-coding, plus a c129 variant that could be contamination or would be irrelevant if genuine).
I haven’t inspected CJD 6, 23, 25, 27 yet, because the vcf gave me a Java error. Cause: the pipeline worked with empty files (FASTQ were missing, probably because the disc was full – fixed now)

## 10.01.24

- Ran CJD 6, 23, 25, 27 in standard settings: no mutations

## 11.01.24

- FilterMutectCalls
- I’ve raised -f-score-beta (less false negatives, more false positives) to see whether the output is different

![Embedded image](analysis_journal_media/image4.png)

    - CJD27 with Fbeta = 1, 2, 5. In all inspected samples (first_CJD_seq) the result looks exactly the same

## To do Mutect

- Review the vcfs with higher Fbeta (I’ve run through all samples)
- Run again with lower threshold.
- In FilterMutectCalls you can increase the -f-score-beta parameter from its default of 1 to increase sensitivity at the expense of precision.
- More information is available in the Mutect2 Tool Index
- https://bioinformatics.stackexchange.com/questions/22033/how-to-adjust-cutoffs-in-mutect2
- Could be useful: --intervals / -L : One or more genomic intervals over which to operate (e.g. only do ORF)
- Run again in multi-sample mode and/or PoN.
- https://gatk.broadinstitute.org/hc/en-us/articles/360050722212-FAQ-for-Mutect2
- Can I specifically look for E200K and correlate the result with the ddPCR data?
- My question: https://www.biostars.org/p/9584544/
- Also look for the FFI mutation in the thalamic CJD samples
- Or look for all mutations that are pathogenic according to the Minikel paper

## Funcotator

- Funcotator (FUNCtional annOTATOR) analyzes given variants for their function (as retrieved from a set of data sources) and produces the analysis in a specified output file.
- This tool is a functional annotation tool that allows a user to add annotations to called variants based on a set of data sources, each with its own matching criteria.

## To do MosaicForecast

- Install MosaicForecast (https://github.com/parklab/MosaicForecast/blob/master/FAQ.md)

### Generate variant call set as input of Mosaic Forecast with Mutect2

https://github.com/parklab/MosaicForecast/blob/master/FAQ.md
- Step1: Mutect2 in tumor-only mode for each normal sample
- It would make sense to provide an interval list that restricts the analysis to the regions targeted by sureselect (first I’ll try without, make this later)
- https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists

## 03.03.2024 (mutation_finder.sh)

Find reads containing specific mutations (position and VCF notation given).
M129V: chr20:4699605 A>G

E220K
- Chromosome 20:4699818
- VCF:20  4699818  rs28933385  G  A

I’ve created a script find_E200K.sh which takes the BAM file for a specific sample and extracts the bases at a position of interest (here: E200K, denoted by the coordinates above). It then outputs the allele depth (AD) of both the REF (here: G) and ALT (here: A) bases into the file defined as $variantVCF. The AD is within the long string at the bottom of the VCF file. 
Different approach. This might be easier to loop:
https://gatk.broadinstitute.org/hc/en-us/community/posts/23439787101595-Is-there-a-GATK-tool-for-quantifying-a-particular-base
https://bioinformatics.stackexchange.com/questions/22237/how-do-i-quantify-a-specific-somatic-variant/22238#22238

Next steps:
- Create loop that processes all samples for a given mutation.
- The variable definitions should be derived from a file that contains the source BAM files (format: path/to/bam) and sample name.
- Then enclose it in a second loop. 
- Derive the definitions from a file with mutation name, position, REF and ALT.
- Mutations: pathogenic muts from Minikel paper, M129V
- If you use the Stackexchange approach, you can use your old code to create a VCF file for every mutation
- Combine all VCF for a given mutation to one file (there’s probably a GATK tool that can do this)
- Create a table (e.g. as TSV file) that contains sample name, CHROM, POS, ID, REF, ALT, QUAL, FILTER, GT, PL and the AD for REF and ALT. Chatgpt created a Python approach for this (C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2024-03-03_python_VCF). This needs understanding and testing.
- GT (genotype) isn’t really needed, but you can check whether the zygosity was interpreted correctly.
- Then use R to do statistics on the AD. Is it higher in CJD? If not, probably no more SNV analysis required.
- Then proceed to Mosaic Forecast.
- Show Christos the analyses.

## 04.03.2024

- Following this: https://bioinformatics.stackexchange.com/questions/22237/how-do-i-quantify-a-specific-somatic-variant I created "/add/2024-03-04_force_mutect/find_E200K_new.sh" which forces Mutect2 to call my variant of interest.
- It works, but seems to filter the reads so that the AD values are lower than the values produced by the approach on 03.03.2024.
- Idea: make stats with both the high-quality (Mutect2) and low-quality (bcftools with high -d parameter) and compare the two.
- But perhaps this isn’t worth doing, as the reads that aren’t counted are skipped by default.
- See also here: https://gatk.broadinstitute.org/hc/en-us/articles/360035532252-Allele-Depth-AD-is-lower-than-expected.

## 05.03.2024

- I created find_E200K_loop.sh. This takes the file names from samples.tsv and runs the E200K command on all of them.
- Now I need to find a way to summarise the allele depths in the VCFs (https://bioinformatics.stackexchange.com/questions/22242/how-do-i-summarise-vcf-files)
- In the meantime, I can create similar loops for the other variants.

### Variants:

- Created list of variants with coordinates, based on the Minikel paper: variant_list.xlsx (/CJD PRNP/Experiments/SureSelect-sequencing/Analysis/2024-04-05_compile_VCF_for_all_variants).
- Next, run the mutation_finder script (with bcftools)

## 06.03.2024

- Summarised the VCFs for E200K (according to the SE answer above)  E200K_summary.tsv
- Some CJDs appear to have about 1/1000 ALT.
- Also a few Ctrls.
- Several samples have poor read counts (why? Batch effect?)
- Create some R stats from the summary file: C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2024-04-06_R_read_files
- The AAF of E200K seems to be in a similar range as the ddPCR results. The patients with higher ddPCR AAF don’t seem to have higher Seq AAF, but I need to check this in detail (also, Seq is only from frontal lobe, whereas the ddPCR values are averaged across 6 samples). 
- Now I want to do the same for the other variants. The trouble is that I haven’t got VCF files for the other variants. I made a custom one with data from Ensembl (http://www.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=20:4699252-4700252;v=rs74315403;vdb=variation;vf=1054935347), but Mutect gives me errors. I’ve posted forum questions: 
- https://bioinformatics.stackexchange.com/questions/22237/how-do-i-quantify-a-specific-somatic-variant/22238?noredirect=1#comment34461_22238
- https://www.biostars.org/p/9589280/

## 08.03.2024

Create VCF using bcftools
- I continued my attempts to create VCF files for my variants of interest. The problem is that I can’t call e.g. D178N (chr20: 4699752 G>A) if it isn’t actually present. So, if I use bcftools as in /add/code/old_code/2024-03-08_bcftoolsmpileup_test.sh, it won’t actually produce an output file.
- Conclusion: It’s probably best to create a VCF manually using data from Ensembl.

Create VCF manually
- Finally succeeded in doing this. They key it to make sure that the values are separated by tabs. You can check this with: zgrep '^#C' my.vcf.gz | od -c

Call D178N mutation
- My script just creates empty VCF files. Maybe because the variant is not present? But it should produce an entry like 2000;0 (many REF reads, 0 ALT)
M129V polymorphism
- Indeed, the script only works if the samples are heterozygous. I asked an additional question: https://bioinformatics.stackexchange.com/questions/22237/how-do-i-quantify-a-specific-somatic-variant
- Test: see if it works if I use a non-manually created VCF  it does
- Manual: M129V.vcf.gz
- Generated: CJD6_M129V_frequency.vcf

## 10.03.2024

- I tried editing the manual input VCF (here, for E200K) in various ways (change sample, FORMAT, add headers) -> then ran "/add/code/mutect_force_variant/E200K_test.sh" with "/add/results/2024-03-04_force_mutect/2024-03-10_E200K_VCFtest/" as output -> didn’t work
- Asked question in GATK forum: https://gatk.broadinstitute.org/hc/en-us/community/posts/23669154192155-Mutect2-forced-calling-via-alleles-does-not-always-generate-an-output
- Apply FilterMutectCalls to my forced VCFs? According to the post above, variants with a high probability of being germline are filtered, so maybe that wouldn’t help (reasoning: very low AAF as in my case -> probably germline)
- Next step: use the “hack” version of the VCF (E200K VCF modified with the D178N variant) for all variants.
- See question above: to troubleshoot, try running it with the E200K header, INFO, FORMAT field (1 option at a time).
- Discussion with Nils:
- GenomicRanges, MutationalPatterns
- eQTL: http://bioinfo.szbl.ac.cn/scQTLbase/Search/ (if the variants are annotated: rs…)

## 13.10.2024 

### Replicate E200K allele counting process

To refresh my memory, I replicated the process I used to find reads with my mutation of interest. I succeeded in doing so for E200K.
- the code that produces the E200K table is on Ubuntu 

```
"/add/code/2024-10-14_E200K_allele_count/"
```

On Windows, I have created a folder 

```
\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\mutation_read_results 
```

where I deposit the result of the read counting procedure. Here, there is also an R script that creates a boxplot with the frequencies. I considered using Fisher’s exact test to compare the frequency of “mutated” status, but the sample sizes are probably too small for this to be legitimate.

### Proceed with D178N allele counting

I did the same for D178N, using the “hack” version of the mutation VCF

### Created loop that counts all variants

- a list of variants is found here:

```
 \CJD PRNP\Experiments\SureSelect-sequencing\Analysis\variant_list.tsv
```

I placed it here:

```
/add/variant_VCF/variant_list_short.tsv
```

- also, a short version with only the 4 highly pathogenic mutations: variant_list_short.tsv
- In add/code, I created code that counts E200K and D178N alleles
- I also created a loop that counts all desired variants (2024-10-13_all_mutation_counter)
- I added M129V, as that had been omitted from the loop
The results from the variant counting are found here: 

```
Directory: /add/results/
2024-10-13_replicate_E200K
2024-10-13_count_D178N
2024-10-13_countAll
2024-10-13_M129V
```

### Summarise counts of individual variants to single file

Code:

```
/add/code/2024-10-13_all_mutation_counter/2024-10-13_all_mutation_counter.sh
```

Results: 

```
/add/results/2024-10-13_VCF_summary/all_summaries/
```

## 19.10.2024

I created an R pipeline that creates a summary table of all reads supporting mutations (\Analysis\2024-10-19_read_count_summaries).

Samples with slightly elevated AAF:
CJD29: P102L (0.52%)
CJD25: E200K (0.11%)

### 21.10.2024

- I’m not sure whether the UMIs have been accounted for.
- Backup of old pipeline: \CJD PRNP\Experiments\SureSelect-sequencing\Analysis\code_backup
- Pipeline v1 is saved here now
- I checked whether UMIs are present in my reads (10 bp at the beginning and at the end). These regions both had only 12% uniqueness, suggesting that they are not UMIs (we’d expect very high uniqueness.
- There is no information in the header that UMIs were used to remove duplicates. However, I entered the 10x N (NNNNNNNNNN) when ordering the sequencing, so I assume that they were processed.
- I also examined the Picard .metrics file. I have a high duplicate rate, which is probably due to the target sequence being small and too many PCR cycles (I had to use many cycles to obtain enough DNA to work with)
- There are also many supplementary alignments. Maybe this is because a large part of the target (the intron) contains repetitive sequences.

Next steps:
Create a table with all the mutations for which reads were analysed, with Position, REF, ALT.  This info should be present in the input for the R pipeline above. Otherwise, you can find the VCFs in /add/variant_VCF/

What shall I do about the samples with low coverage?
- Inspect the quality, see script here: /add/2024-10-21_examine_CJD5_quality/
- Exclude the poor quality files from the read counting section. Poor read Per base sequence quality, possibly due to residual phenol
- Set a limit at 100-1000x coverage, based on limit of detection in spike-in exp
- I think I read that BWA-Mem is better at low-quality bases and/or repetitive regions: compare with Minimap2

## 27.10.2024

I would like to also run MosaicForecast. I executed the commands on the GitHub website: https://github.com/parklab/MosaicForecast

This explains how to “pull” MF from docker. It then unpacks the human genome hs37 and runs a demo code (last line):

```
docker image pull yanmei/mosaicforecast:0.0.1
docker run -v ${your_local_directory}:/MF --rm -it yanmei/mosaicforecast:0.0.1 /bin/bash
gunzip hs37d5.fa.gz

Phase.py /MF/demo/ /MF/demo/phasing hs37d5.fa /MF/demo/test.input 20 k24.umap.wg.bw 4
The last line gives me this output:

root@2482a8fe7422:/usr/local/bin# Phase.py /MF/demo/ /MF/demo/phasing hs37d5.fa /MF/demo/test.input 20 k24.umap.wg.bw 4 Usage: python Phase.py bam_dir output_dir ref_fasta input_positions(file format:chr pos-1 pos ref alt sample, sep=\t) min_dp_inforSNPs(int) Umap_mappability(bigWig file,k=24) n_threads_parallel sequencing_file_format(bam/cram) 
Note: 
1. Name of bam files should be "sample.bam" under the bam_dir, and there should be corresponding index files. 
2. There should be a fai file under the same dir of the fasta file (samtools faidx input.fa). 
3. The "min_dp_inforSNPs" is the minimum depth of coverage of trustworthy neaby het SNPs. 
4. Bam file is preferred than cram file, as the program would run much more slowly if using cram format.
```

I need to build my own command following the syntax given at the top.


Next, I started to obtain the resources required for MF to run.

These commands download, extract, and convert mappability data into a bigWig file (k24.umap.wg.bw). This file can then be used in genomic analysis tools that need mappability information to identify unique versus repetitive regions, particularly useful in tasks like phasing.



```
wget https://bismap.hoffmanlab.org/raw/hg38.umap.tar.gz  
tar -zxvf hg38.umap.tar.gz  
cd hg38  
fetchChromSizes hg38> hg38.chrom.sizes  
wigToBigWig <(zcat k24.umap.wg.gz) hg38.chrom.sizes k24.umap.wg.bw  
```

I had to install fetchChromSizes and wigToBigWig using Conda.
Problem: BigWig apparently needs 30GB RAM and I’ve only got 16GB. 

### Workaround for lack of RAM:

Since the process previously used around 7 GB of RAM before being killed (I saw this in the error message), creating a swap file of 16 GB or larger should provide a sufficient buffer. Run the following commands to create a swap file:

```
sudo fallocate -l 16G /swapfile   # Creates a 16GB swap file
sudo chmod 600 /swapfile          # Sets the correct permissions
sudo mkswap /swapfile		# Format the new file as swap space
sudo swapon /swapfile	# Activate the swap file so the system can start using it immediately
swapon –show	#Check to confirm that the swap has been added successfully
free -h
```

In the output, you should see the additional swap space listed. The free -h command will display both your RAM and swap space, showing the total memory available for processes.

![Embedded image](analysis_journal_media/image5.png)

## 03.11.2024

The bigwig command seems to have worked using the RAM workaround.

Fetching other files I need (GRCh38/hg38):
Segmental Duplication regions (should be removed before calling all kinds of mosaics):

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
```

Regions enriched for SNVs with >=3 haplotypes (should be removed before calling all kinds of mosaics): This command did not work  404 not found

```
wget https://raw.githubusercontent.com/parklab/MosaicForecast/master/resources/predictedhap3ormore_cluster.GRCh37.bed
```

ChatGPT: The file predictedhap3ormore_cluster.GRCh38.bed is not available in the MosaicForecast repository. However, you can use the SegDup_and_clustered.GRCh38.bed file, which combines segmental duplication regions with clustered single nucleotide variant (SNV) regions for the GRCh38 reference genome. This file is suitable for excluding regions prone to artifacts when calling mosaic variants.

```
wget https://raw.githubusercontent.com/parklab/MosaicForecast/master/resources/SegDup_and_clustered.GRCh38.bed
```

This could indeed be downloaded.

Simple repeats (should be removed before calling mosaic INDELS):

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
```

### How to run MF

Usage:
python Phase.py bam_dir output_dir ref_fasta input_positions min_dp_inforSNPs Umap_mappability(bigWig file,k=24) n_threads_parallel sequencing_file_format(bam/cram)
Note:
Name of bam files should be "sample.bam" under the bam_dir, and there should be index files under the same directory (samtools index sample.bam).
There should be a fai file under the same dir of the fasta file (samtools faidx input.fa).
File format of the input_positions: chr pos-1 pos ref alt sample, sep=\t
The "min_dp_inforSNPs" is the minimum depth of coverage of trustworthy neaby het SNPs, can be set to 20.
The program to extract mappability score: "bigWigAverageOverBed" should be downloaded and installed, and its path should be added to the PATH environment variable.

I’ve started contructing a command in add/MosaicForecast. I just realised that this tool doesn’t look for new SNPs (as far as I can tell), but examines SNPs which I need to define myself. 
Idea: run MF with the 4 highly pathogenic variants.

Finish setup:
Install all the python packages listed at the top of the manual.
I think I’ve succeeded, with the exception of “pysamstats”.
Checked on 26.01.2025: pysamstats is installed within the Conda environment.

## 26.01.2025

Test run of MF
I verified that the directories specified at the top of the script are correct.
Added error messages for missing directories.
I tried to run the script ("/add/MosaicForecast/runMF_errorcheck.sh"), however, this generated a recursion error. According to ChatGPT, this could be caused by old versions of the required Python packages.
Solution: Use Docker image in the documentation that has got the correct versions

- I mounted /add/MosaicForecast from the host to /MF in the container.

To run the analysis, you need the following:
- /add/MosaicForecast/input/ must contain
- chr2_chr4_chr20.fasta
- chr2_chr4_chr20.fasta.fai
- sample.bam (can be renamed e.g. CJD1.bam, as long as it’s referenced in input_positions.tsv)
- sample.bam.bai (likewise)
- input_positions.tsv
- input_positions.tsv
- no header, columns are chr pos-1 pos ref alt sample
- make sure that the columns are separated by tabs (view -> special chars)

```
chr20	4699817	4699818	G	A	CJD1
chr20	4699524	4699525	C	T	CJD1
chr20	4699569	4699570	C	T	CJD1
chr20	4699751	4699752	G	A	CJD1
```

- Command to start the Docker container (saved in "/add/MosaicForecast/2025-01-26_MF_using_docker"):

```
sudo docker run \
    -v /add/MosaicForecast:/MF \
    --rm -it \
    yanmei/mosaicforecast:0.0.1 /bin/bash
```

Now the user is root@8818837e2f6c:/usr/local/bin#
- Correct directory:

```
cd /MF
```

- Now execute MosaicForecast:

```
python /usr/local/bin/Phase.py input/ output/ input/chr2_chr4_chr20.fasta input/input_positions.tsv 20 k24.umap.wg.bw 4 bam
```

This generated a number of files in output:
- all_2x2table: Contains 2x2 contingency tables for the phasing process.
- all_candidates: Lists candidate mosaic variants identified during the analysis.
- all.merged.inforSNPs.pos: Consolidated positions of nearby informative SNPs.
- all.phasing: The main phasing output file with results for all variants processed.
- all.phasing_2by2: Phasing results tied to the 2x2 contingency tables.
- multiple_inforSNPs.log: A log file detailing results for sites with multiple nearby informative SNPs.

For CJD1, these files were empty, likely indicating that no mosaic variants are present.

## 02.02.2025

Researching options to improve quality: o3 suggestions
  - more aggressive trimming with Trimmomatic
  - bwa-mem instead of Minimap2
Now, I’ll run the QC step on all samples. Then I can try the new pipeline on selected samples and evaluate the effect.

### QC of all samples

- In /add/2025-02-02_examine_read_quality/ I created examine_sample_quality.sh, which loops through all markedDup.recal.bam files to create QC reports

## 03.02.2025

- I also need QC reports on the spike-in samples. Done with "/add/2025-02-02_examine_read_quality/examine_sample_quality_v2.sh"

Summarise data with MultiQC
- Installed MultiQC

```
multiqc /add/2025-02-02_examine_read_quality/
```

- Saved output: multiqc_report_original_pipeline.html in here

## 07.02.2025

- I’ve run the more aggressive pipeline on all samples overnight
- Code: /add/code/2025-02-07_aggressive_pipeline/
- Output: 
- "/add/seq_data/2025-02-05_dilutions_bwa/"
- "/add/seq_data/2025-02-05_samples_bwa/"

## 08.02.2025

- I need a strategy for comparing the performance of these pipelines
- To do: run Mutect and MF on the dilution samples from old and new pipelines

Mutect
- I’ve run Mutect on the dilutions and samples from the new pipeline
- Code: /add/code/mutect_bwa/
- Results: /add/results/mutect_output/dilutions_bwa and samples_bwa

### Evaluate Statistics and Read Quality

QC strategy

When comparing two pipelines—one that uses your original workflow and one that uses your aggressive preprocessing—you want to assess which one provides a more accurate and reliable set of somatic variant calls. A systematic evaluation generally involves several steps and metrics. Here’s a detailed strategy:

1. Evaluate Alignment and BAM File Quality
Mapping Statistics:
Coverage and Uniformity: Use tools like samtools stats or Picard’s CollectAlignmentSummaryMetrics to assess overall coverage, coverage uniformity, and duplication rates. Higher, more uniform coverage is typically desirable.
Mapping Quality Distribution: Compare the distribution of mapping quality scores between the two pipelines. Aggressive trimming might remove low-quality bases that would otherwise reduce mapping quality, but it can also lead to shorter reads and potential loss of coverage.
Insert Size Metrics: Use Picard’s CollectInsertSizeMetrics. Changes in insert size distribution may indicate differences in how reads are trimmed and aligned.

2. Compare Variant Calling Metrics in the VCFs
Variant Statistics:
Total Number of Calls: Compare the total number of somatic variants called by Mutect2 in each pipeline. An excessive number of calls might indicate an elevated false positive rate.
Quality Scores: Examine Mutect2 quality metrics (e.g., TLOD scores, NLOD scores, allele frequencies) in each VCF. Higher confidence metrics and a clear separation between likely true calls and background noise can be an indicator of better preprocessing.
Variant Allele Frequency (VAF) Distribution: Plot and compare the VAF distributions. With aggressive trimming, you might see a shift in VAFs because lower-quality regions (which might create false positives) are removed.
Filter Metrics:
Filter Pass Rates: Look at how many variants pass all filters in each VCF. A lower pass rate might indicate that many calls are borderline, whereas a pipeline with clear separation between true positives and noise might have more calls passing stringent filters.

3. Use Spike-In Controls or Known Variants
If you have spike-in samples (or other controls with known somatic mutations), these provide a gold standard for comparison:
Sensitivity (Recall): Evaluate how many of the known spike-in variants are detected by each pipeline.
Specificity (Precision): Check for false positive calls. A pipeline that produces fewer spurious calls in regions where no mutations are expected is generally preferable.
Correlation Analysis: Compare the measured allele frequencies in the spike-in experiments with the expected frequencies. More accurate VAFs suggest that preprocessing has effectively removed errors while retaining true signals.

4. Qualitative Review Using Genome Browsers
Visual Inspection: Load representative BAM files and VCFs into a genome browser (e.g., IGV). Look at regions with variant calls from both pipelines. Determine whether the aggressive trimming leads to better-resolved alignments and more convincing variant evidence (e.g., clear split-read or discordant pair support) compared to the original pipeline.
Edge Cases: Particularly check regions that are problematic in the original pipeline—such as those with lower base quality—and see if aggressive preprocessing improves clarity.

5. Downstream Impact
Consider the broader biological context:
Reproducibility: Evaluate whether the variant calls are reproducible across technical replicates or similar sample types.
Functional Annotation: Run both sets of VCFs through your downstream annotation pipelines (e.g., using tools like ANNOVAR or VEP) and see if the biological interpretation (e.g., predicted impact on protein function) is more consistent or plausible with one approach.

6. Statistical and Computational Comparison
Comparison Tools: Use tools such as VCFtools, bcftools stats, or even benchmarking tools like hap.py to generate summary statistics and compare variant call sets quantitatively.
Concordance Analysis: Measure the overlap between the two call sets. Are there systematic differences (e.g., one pipeline consistently calls more low-frequency variants that are later filtered out)?

Final Considerations
Trade-Offs: Aggressive trimming can remove low-quality bases and reduce false positives, but it may also lead to loss of coverage in critical regions, thereby missing true low-frequency mutations. On the other hand, the original pipeline might retain more reads but at the cost of higher noise.
Holistic Assessment: No single metric tells the whole story. You should combine alignment quality metrics, variant calling statistics, control-based evaluations, and manual review to decide which pipeline performs best for your biological question.
By following this multifaceted evaluation strategy, you can determine which pipeline—original or aggressive preprocessing—provides the most reliable, reproducible, and biologically meaningful results for your somatic mutation analysis.

Samtools Flagstat Analysis

```
samtools flagstat your_sample.bam > your_sample.flagstat.txt
```

I created loops to run flagstat on the dilution BAM files (BAM that resulted from the 2 pipelines)

code: /add/code/pipeline_comparison/2025-02-08_flagstat_dilutions/
Results: /add/results/pipeline_comparison/2025-02-08_flagstat_dilutions/
The flagstat results are saved with “orig” for the original pipeline and “bwa” for the new aggressive pipeline.
The AWK script parse_flagstat.awk summarises the parameters in the flagstat result files (run_parse_flagstat.sh):

```
cd /add/results/pipeline_comparison/2025-02-08_flagstat_dilutions/
awk -f "/add/code/pipeline_comparison/2025-02-08_flagstat_dilutions/parse_flagstat.awk" *.flagstat.txt > summary_flagstat.tsv
```

Samtools Flatstat Results
Alignment Quality:
The aggressive pipeline (bwa) appears to improve mapping efficiency and reduce the number of singletons and duplicates. These are positive indicators for overall alignment quality.
Trade-Off:
The cost of these improvements is that many more reads are discarded (leading to fewer total reads) and there is a higher incidence of discordant mappings. You’ll need to determine whether the loss of some reads is acceptable in exchange for cleaner alignments.
Next Steps:
Visual Inspection: Load representative regions in a genome browser (like IGV) to check if the increased discordant mappings in the aggressive pipeline represent true structural variation or mapping artifacts.
Variant Calling Impact: Compare the downstream variant calling results. If the aggressive pipeline’s improvements in mapping and duplicate reduction lead to more accurate somatic mutation detection (for example, validated via spike-ins or controls), it might be the better option.
Reproducibility & Sensitivity: Assess sensitivity (true positive rate) and specificity (false positive rate) using known mutations or spike-in controls.
In summary, the aggressive pipeline improves mapping percentages and reduces problematic singletons/duplicates—but it also introduces more discordant mappings. These differences should be evaluated in the context of your variant calling performance and the biological plausibility of the discordant events.
Summary of what I’ve done with Flagstat
- Code: /add/code/pipeline_comparison/
- Results: /add/results/pipeline_comparison/
- I’ve run flagstat on all dilutions and samples, old and new pipelines
- Summarised the results with the AWK script  summary_flatgstat.tsv

### Interpretation (o3)

Improved Mapping Quality:
The aggressive pipeline (bwa) consistently produces a higher proportion of mapped and properly paired reads relative to its total. Even though the total read count is lower, the reads that remain appear to be of higher quality, which is likely beneficial for downstream variant calling.
Reduction in Duplicates and Singletons:
Lower duplicate and singleton counts in the aggressive pipeline indicate that it removes a significant amount of ambiguous or low-quality data, thereby potentially reducing false positives in variant calling.
Trade-offs:
The aggressive trimming comes at the cost of losing some total read depth. In many cases, however, the benefit of having a higher quality, more confidently mapped set of reads can outweigh the loss in coverage—especially if your downstream analyses (such as somatic variant calling) rely on accurate alignment rather than sheer depth.
– The higher numbers of discordant mappings in some “bwa” samples need additional evaluation. If these reflect true biological signals (e.g., structural variants) or if they are confidently mapped (as suggested by higher mapQ counts), then this may be an advantage. If they are artifacts of over-trimming, they could be a concern.
Next Steps
Visual Review:
Use a genome browser (e.g., IGV) to inspect regions with discordant mappings, especially in samples where “bwa” shows a marked increase in reads with mates on different chromosomes.
Downstream Impact:
Ultimately, compare variant calling results (e.g., sensitivity and specificity using known spike-ins) between the two pipelines. Higher mapping quality and lower duplication rates should, in principle, translate into more reliable mutation detection.
Statistical Analysis:
Consider normalizing metrics (e.g., percentage mapped, percentage properly paired) and plotting these values for each sample to visually compare the performance of the two pipelines across the board.

![Embedded image](analysis_journal_media/image6.png)

Conclusion: Percentage mapped and pct properly paired are generally better in bwa (or rather: the aggressive trimming helped to achieve this). Duplicates are similar, but I expect many duplicates due to the many PCR cycles, so this is less important.

## 09.02.2025

FastQC of all samples
- Code: "/add/code/pipeline_comparison/2025-02-09_fastQC_all_samples/run_FastQC_comparison.sh"
- Results: /add/results/pipeline_comparison/2025-02-09_FastQC_samples/
MultiQC
- Run separately on bwa and original outputs (same dir)

## 10.02.2025

### Run Mutect forced-call script on bwa output

- This is the script that looks for specific mutations and force-calls them.
- Code: "/add/code/2025-02-10_bwa_force_mutect/runForceMutect_bwa.sh"
- Results: /add/results/2025-02-10_force_mutect_bwa/

### Run MF as loop (BWA samples)

- Code: "/add/code/2025-02-10_bwa_MF_loop/"
- Setup:  run move_files_for_MF.sh to move the BAM and BAI to be analysed to the input directory used by the docker container
- Run analysis: run_MF_loop.sh. This loops through every sample by changing the tsv
- retrieve sample names (e.g. CJD1_bwa, CJD14_bwa, Ctrl7_bwa) from the input directory, then start loop for first sample
- modify "/add/MosaicForecast/input/input_positions.tsv" so that the name column contains the sample name (e.g., if the script is working on CJD1_bwa, change the name column to "CJD1_bwa", where there should be an entry in every row, in this case 4 rows)
- run MF on that particular sample
- copy the six output files to "/add/results/2025-02-10_MF_loop_bwa/", inside a directory for every sample (e.g., for CJD1_bwa, the directory is "/add/results/2025-02-10_MF_loop_bwa/CJD1_bwa". The six output files should have the sample name in their filename (e.g. CJD1_bwa_all_2x2table and CJD1_bwa_all_candidates) 
- proceed with next sample
- Results: /add/results/2025-02-10_samples_bwa_MF_loop
- Clear-up: remove BAM and BAI from the input dir

### Run MF as loop (original samples)

- Similar procedure as above
- Setup: "/add/code/2025-02-10_samples_original_MF_loop/move_original_files_for_MF.sh"
- Code: "/add/code/2025-02-10_samples_original_MF_loop/run_MF_loop_original.sh"
- Results: "/add/results/2025-02-10_MF_loop_original/"

### Run MF on dilutions (bwa)

- Code: /add/code/2025-02-10_dilutions_bwa_MF_loop/
- Results: /add/results/2025-02-10_MF_loop_bwa_dilutions/

### Run MF on dilutions (original)

- Code: /add/code/2025-02-10_dilutions_original_MF_loop/
- Results: /add/results/2025-02-10_MF_loop_original_dilutions/

## 23.02.2025

### Summarise forced call VCFs

- Adapted summarising script created on 13.10.2024
- Code: "/add/code/2025-02-23_summarise_VCF_bwa/2025-02-23_summarise_VCF_bwa.sh"

- Input Directories:
- The input files are stored under /add/results/2025-02-10_force_mutect_bwa, with each subdirectory named after a mutation (e.g., “E200K”, “F198V”, etc.).
- Output Directories: 
- Two main output locations are defined:
One for storing summarized, compressed, and indexed VCF files at /add/results/2025-02-23_bwa_forced_calls_VCF_summary
One for merged VCF files at /add/results/2025-02-23_bwa_forced_calls_merged_VCF
Additionally, a subdirectory all_summaries is created to store summary TSV files.
- Variant List File:
The script reads a variant list from /add/variant_VCF/variant_list_short.tsv. This file (after skipping the header) provides details for each mutation (e.g., variant name, dbSNP ID, position, reference, and alternate allele).

Based on 19.10.2024, I created an identical R script to summarise the results:
- Dir: C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-02-23_bwa_read_count_summaries
- Script: summarise_reads_stats_bwa.R

Next step: compare the outputs

## 04.03.2025

### 
Compare forced call VCF summary CSVs read counts

I need to compare the old and the “aggressive” bwa pipelines. For both, I have a large CSV with read counts for the REF and ALT allele for every mutation and sample.
To compare, subtract values: bwa – old  see where counts are better, particularly in the samples with lower coverage. I created an R script for this in this dir:

Dir: C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-03-04_compare_read_counts
Result: bwa has got higher read counts in most samples

### 
Find samples where AAF is over detection limit

The detection limit should be based on the detection limit from your dilution samples.
Added script to the R code above. There seem to be almost no reads supporting mutations (perhaps 0.1% or so).
Next: run force call script on dilution samples

## 15.03.2025

### Ensure that all mutation algorithms have been run on the dilution samples

- Bwa + Mutect2 regular: "/add/results/mutect_output/2025-02-08_dilutions_bwa/"
- Original + Mutect2 regular: "/add/results/mutect_output/2021-02-03_sequencing_of_dilutions/"
- Bwa + Mutect2 force call: "/add/results/2025-03-15_bwa_force_call_dilutions/"
- Original + Mutect2 force call: "/add/results/2025-03-15_original_force_call_dilutions/
- Bwa + MosaicForecast: /add/results/2025-02-10_MF_loop_bwa_dilutions/
- Original + MosaicForecast: /add/results/2025-02-10_MF_loop_original/

### Summarise force-calls of dilutions

- Original:
- Code: "/add/code/2025-03-15_original_force_call_dilutions/2025-03-15_original_force_call_dilutions_summarise.sh"
- Results: /add/results/2025-03-15_original_forced_calls_merged_VCF_dilutions/
- Bwa: 
- Code: "/add/code/2025-03-15_bwa_force_call_dilutions/2025-03-15_bwa_force_call_dilutions_summarise.sh"
- Results: "/add/results/2025-03-15_original_forced_calls_VCF_summary_dilutions/"

As on 23.03.2025, I created R scripts to create summary tables (C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-03-15_dilution_force_call_summaries)
I created summary tables for both original and bwa. They detect the A mutation spike-in in the respective dilution sample.

## 22.03.2025

#### Comparison of force call results

- BWA has better read counts
- Detection limit: 0.006  0.6%, based on NA99A1_1to5 (was supposed to be 1% spike-in, approximately correct)
- NA995A05: was supposed to have 0.5% spike-in, but 0.014 AAF was found (probably inaccurate spike-in)

#### FilterMutectCalls

The VCFs produced from the dilution samples contain many variants that are probably false positives. I should run FilterMutectCalls
Code: /add/code/2025-03-22_FilterMutectCalls/
results: File containing “_filtered”, saved in the same directory as the zipped input VCF

#### Comparison of Mutect2

I’ll start by comparing the filtered NA99A1_1to5 VCF in IGV for both the bwa and original pipelines to see whether the Axx mutation was detected correctly and whether there are additional variants that might be false positives. If there are, I’ll need to see whether they can also be detected in the NA sample and in the pure spike-in. If not, they are probably false positive. 

- In the bwa version, the A117V spike-in is detected by Mutect and is not filtered out by FilterMutectCalls  success
- In the original version, it is not (which is strange)
- Therefore, I’ll proceed with the bwa version
- I can’t see any potential false positives

![Embedded image](analysis_journal_media/image7.png)

Figure: In the bwa variant, the A117V spike-in is detected by Mutect2 and survives filtering (bottom track, the MV polymorphism is filtered out)

#### Compare MF output

When I ran MF, I made it specifically test for the four highly pathogenic variants

For NA99A1_1to5 bwa, the spike-in was detected:

```
NA99A1_undil_bwa chr20 4699818 G A 3180 2 0 SNP
NA99A1_undil_bwa chr20 4699525 C T 2643 0 0 SNP
NA99A1_undil_bwa chr20 4699570 C T 3476 24 0 SNP
NA99A1_undil_bwa chr20 4699752 G A 3111 1 0 SNP
```

A117V

```
#CHROM	POS	ID	REF	ALT
chr20	4699570	rs74315402	C	T
```


However, the phasing and 2x2table files are empty. Probable reasons:


Would one expect empty outputs from MosaicForecast here?
Empty *_2x2table file:
At very low VAF (≤ ~1%), MosaicForecast often does not detect sufficient statistical evidence to confidently report allelic imbalance. The test requires clear deviation from expected allele counts to populate this table. With only ~24 alternate reads out of 3476 total reads (0.69% VAF), statistical power is limited, and MosaicForecast may not output entries here.
Empty *_phasing file:
MosaicForecast requires informative SNPs nearby to phase variants. An artificial spike-in sample often lacks informative heterozygous variants nearby, resulting in no phasing information. Hence, empty phasing outputs are entirely expected.

Interpreting your scenario clearly:
Your main variant (A117V) is indeed detected in the *_all_candidates file (24 reads, ~0.69% VAF), exactly matching your expectations.
Absence of additional statistical and phasing evidence is consistent with the very low allele frequency and lack of informative SNPs.


Possible conclusion: Limited sensitivity for MF in detecting my variant of interest.

NA100_1to10 (no spike-in control)

```
chr20 4699818 G>A | DP=5913 | ALT=1 read | VAF=0.017%
chr20 4699525 C>T | DP=4558 | ALT=2 reads | VAF=0.044%
chr20 4699570 C>T | DP=5227 | ALT=1 read | VAF=0.019%
chr20 4699752 G>A | DP=5780 | ALT=1 read | VAF=0.017%
```


Conclusion:
- ~0.02–0.05% background sequencing noise (1–2 ALT reads per ~5000 depth).
- 0.6% can be clearly detected

#### Summarise MF output

- Start by inspecting the MF output files to understand their structure.
- Created a script, but it calculates VAF incorrectly, you need to check it: "/add/code/2025-03-22_summarise_MF/summarise_MF_dilutions_bwa.sh"

## 29.03.2025

- I have analysed the MF code for phase.py to try to understand the output that is written into the all_candidates files (https://github.com/parklab/MosaicForecast/blob/master/Phase.py)
- The columns in the all_candidates files appear to be the following: 
- Sample
- Chromosome
- Position
- Major allele (often representing the reference allele)
- Minor allele (typically the alternate allele)
- Major read count (number of reads supporting the major allele, akin to REF_count)
- Minor read count (number of reads supporting the minor allele, similar to ALT_count)
- Conflict count (the number of reads that support both alleles)
- Variant type (e.g., SNP, MNP, DEL, INS)

- I have edited the code to create a summary of all MF all_candidate files, with the following columns:
- Sample, chr, position, REF, ALT, REF_reads, ALT_reads, conflict_count, variant_type

### Mutation algorithms: summarised results for dilution samples

- Bwa + Mutect2 force call: "/add/results/2025-03-15_bwa_forced_calls_merged_VCF_dilutions/"
- Original + Mutect2 force call: "/add/results/2025-03-15_original_forced_calls_merged_VCF_dilutions/"
- Bwa + MosaicForecast: /add/results/2025-03-22_MF_dilution_summary/
- Original + MosaicForecast: /add/results/2025-03-22_MF_dilution_summary/
For the regular Mutect output, I haven’t got a summary, but the FilterMutectCalls results that can be evaluated in IGV (e.g. o23928_1_7-NA995A05_undil_G01.Mutect_filtered.vcf.gz)

- Bwa + Mutect2 regular: "/add/results/mutect_output/2025-02-08_dilutions_bwa/"
- Original + Mutect2 regular: "/add/results/mutect_output/2021-02-03_sequencing_of_dilutions/"
The summaries (without Mutect2 regular) have been added to:

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\final results\dilution samples
```

### Choose bwa or original


- MosaicForecast:
- The depth seems to be similar for both. No clear difference. 12 rows have higher depth in bwa, 12 have higher in original. Bwa might have slightly higher counts.
- Mutect2 original + FilterMutectCalls: Bwa seemed to perform better, see 22.03.2025
- Mutect2 force call: Bwa has higher read counts 
- see R script here: C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-03-15_dilution_force_call_summaries

-  Conclusion: Work with bwa

Summary of other MF outputs
A script summarises all additional files: /add/code/2025-03-29_summarise_MF_files/
However, with the exception of this file

```
A100_1to2_bwa_all.merged.inforSNPs.pos"
```

Contents:

```
A100_1to2_bwa chr20 4699570 C T chr20 4699570 T C 0 SNP
```

They are empty, so there isn’t anything to summarise.

### Summarise MF for all bwa samples


all_candidates files:

```
"/add/results/2025-03-25_MF_all_bwa_samples/"

```


other files:

```
"/add/results/2025-03-29_MF_files_summary_bwa_samples/"
```

### Interpretation of MF summaries

Empty Files (Phasing, Phasing_2by2, and inforSNPs_log):
The summary files for the phasing outputs—MF_files_summary_phasing.tsv, MF_files_summary_phasing2by2.tsv, and MF_files_summary_inforSNPs_log.tsv—are empty.
This may indicate that for these particular samples (or at least in the BWA dataset analyzed), no additional phasing information was generated or no informative heterozygous SNPs were found that met the criteria for phasing. In other words, the pipeline may not have required phasing corrections, or the data did not produce any output in these steps.
MF_files_summary_merged_inforSNPs.tsv:
For sample CJD29_bwa, there is one entry at candidate position 4699525 on chr20.
The candidate shows a reference allele of C and an alternate allele of T.
The associated informative SNP is reported at position 4699605 with alleles G (inforSNP_ref) and A (inforSNP_alt).
The conflict count is 0, and the variant is classified as a SNP.
This suggests that, at this candidate site, the nearby informative SNP data is available and does not indicate conflicting evidence (i.e. no reads supporting both alleles).
MF_files_summary_2x2table.tsv:
There are two entries in this file.
For CJD29_bwa at candidate position 4699525:
The informative SNP is again at 4699605 with alleles G and A.
The 2×2 table shows very high counts for the major (reference) allele (416 and 395) and zero counts for the minor allele.
This strong and exclusive support for the reference allele at the informative SNP suggests that there is no evidence for a mosaic event at this candidate, reinforcing its classification as a germline SNP.
For CJD30_bwa at candidate position 4699818:
The candidate position and the informative SNP position are the same (4699818) with both alleles matching (G and A).
The absence of additional detailed counts in the summary implies that the signal at this site is consistent with a SNP rather than a mosaic event.
Overall Interpretation:
Phasing Data Absence:
The lack of output from the phasing-related files might indicate that either the candidates did not have sufficient nearby heterozygous sites for a robust phasing analysis or that the pipeline determined phasing corrections were unnecessary for these samples.
Merged Informative SNPs and 2×2 Tables:
For the candidate sites present, the merged informative SNP and 2×2 table outputs suggest that the evidence strongly supports a germline SNP call rather than a mosaic variant. In particular, for sample CJD29_bwa, the high reference allele counts with no evidence of alternative support at the informative site indicate a consistent signal, while for CJD30_bwa, the matching candidate and informative SNP positions further reinforce the SNP classification.
In summary, the data available from the merged informative SNPs and 2×2 tables support the interpretation that the candidates detected in these samples are likely germline SNPs rather than mosaic variants, and the phasing outputs being empty further corroborates that no additional mosaic-specific phasing signal was detected in this dataset.

## 01.04.2025

I need to evaluate the filtered Mutect2 calls. I’ve created a script that automatically extracts the PRNP variants ("/add/code/2025-04-01_extract_filtered_prnp_variants/") to this directory: "/add/results/2025-04-01_filtered_prnp_vcf/"
The script worked. I’d like to create a convenient table. Instead of creating a new script, find a tool that already exists. Gatk VariantsToTable, VCF tools or bcf tools could be good.
Also implement the corrections by chatGPT

05.04.2025
I created boxplots with t-tests for the read counts (Mutect2 force-call) here: 

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\final results\plots
```

I’ve run VariantsToTable to create a separate table for every sample, which are saved in 

```
/add/results/2025-04-05_combined_variant_tables/
```

Next, I want to combine them to a more convenient table. I’ve saved them to Windows:

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-04-05_combined_variant_tables_dilutions
```

R script to combine them. I’ve coded it so that it adds e.g. “E200K position” if the respective position is affected. If it is, you’ll have to check manually in IGV whether it changes the amino acid. This solution avoids having to create a command for every possible nucleotide change.

## 12.04.2025

Instead of the manual solution in the last step, it might be good to directly run a pathogenicity analysis on the entire output of Mutect.
https://chatgpt.com/c/67fa3569-6ef8-8002-a978-13613b6313ce
Pathogenicity analysis pipeline
Preprocess VCF files:
Normalize using bcftools/vt
Apply additional filtering criteria
Annotate Variants:
Use VEP/ANNOVAR/snpEff with databases (ClinVar, COSMIC, gnomAD, etc.)
Predict Pathogenicity:
Run in-silico predictors (SIFT, PolyPhen-2, CADD, etc.)
Integrate cancer-specific annotations
Filter and Prioritize:
Filter by population frequency, impact, and known clinical annotations
Tier the variants
Review and Validate:
Visualize in IGV and perform manual curation
Validate key variants experimentally if required

### Organise filtered calls into single directory (mutect_filtered_collection) to simplify pipeline design

```
cp /add/results/mutect_output/2025-02-08_dilutions_bwa/*filtered* /add/results/mutect_filtered_collection/dilutions/
cp /add/results/mutect_output/2025-02-08_samples_bwa/*filtered* /add/results/mutect_filtered_collection/samples/
```

### Annotation pipeline

I’m building a pipeline that executes all annotation steps. The steps are saved here:

```
/add/code/annotation_pipeline/
```

### Normalisation with bcftools

1_annotation_preprocessing.sh executes bcftools:
Filtered VCF files generated from the Mutect2 and FilterMutectCalls pipeline were normalized using a custom bash script that applied bcftools. The script iterated over two datasets—samples (comprising CJD and control samples) and dilutions (control samples)—located in separate directories. For each file, multi-allelic sites were decomposed using the bcftools norm -m -any option, and the resulting records were then left-aligned based on the reference genome (chr2, chr4, and chr20 from /home/mcarta/databases/chr2_chr4_chr20.fasta) by invoking bcftools norm -f. The normalized variants were compressed, indexed, and stored into dedicated output directories, ensuring uniform variant representation for subsequent functional annotation steps.

### Funcotator (preparation)

https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial
- I have set up the necessary reference files here by downloading them using FuncotatorDataSourceDownloader: "/add/funcotator_data_somatic/hg38"
- I have indexed the vcf files to create a tbi for every sample and dilution using bcftools index

```
for vcf in /add/results/annotation_output/1_normalisation/samples/*.vcf.gz; do bcftools index -t "$vcf"; done
```

```
for vcf in /add/results/annotation_output/1_normalisation/dilutions/*.vcf.gz; do bcftools index -t "$vcf"; done
```

### Funcotator (execute)

Funcotator script: "/add/code/annotation_pipeline/2_funcotator_dilutions.sh" and “2_funcotator_samples.sh”

I’ve run funcotator overnight on both dilutions and samples and it seems to have worked. I still need to inspect the results. However, it seems like it didn’t use gnomAD. I’ve extracted the two tar.gz files related to gnomAD to /hg38 and I’ll rerun the tool.
Including gnomAD seems to greatly lengthen the time required. The manual mentions that this might depend on the speed of the internet connection.

## 13.04.2025

I’ve saved the funcotator output in the directories dilutions_1st and samples_1st for now. 
I previously created a script to convert the VCFs to a table. See 01 and 05.04.2025.
I can also use the R script to check whether the 4 pathogenic variants are found.
Next, I’ll evaluate the funcotator analysis. I’ll need to define filtering parameters to remove variants that are likely to be noise. On the other hand, I’ve already used FilterMutectCalls, maybe that’s enough.
The files generated by funcotator are very large. It’ll be difficult to collate them to a large dataframe, as I originally intended. 
Next steps:
- I’ll start by using my R workflow to search for known pathogenic mutations.
- Next, I’ll have to apply filtering to reduce the amount of data. Probably based on coverage, AAF (I want it be low enough to have a chance of finding results, but I don’t want too much junk either), other quality scores which I’ll have to learn about

### My R workflow to create an easy dataframe (without funcotator annotations)

I’m doing what I did on 05.04.25 with the dilutions, but with the samples.
- Run "/add/code/2025-04-01_extract_filtered_prnp_variants/1_extract_filtered_prnp_variants.sh"
- This extracts PRNP variants with a PASS score only
- Run "/add/code/2025-04-01_extract_filtered_prnp_variants/VariantsToTable.sh"
- Creates convenient tables from the PRNP variants I just extracted
- I moved the tables to C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-04-05_combined_variant_tables\samples
- Summary table with R script in the same directory  no pathogenic variants (other than the heterozygous E200K in CJD30). Otherwise, we have 19 variants in the protein coding region.
- Added the other variants (other than the four highly pathogenic ones) to the R script  still no pathogenic variants
- Some might be heterozygous variants (with AAF of e.g. 0.4) and not somatic. There are about 8 with AAF in the 1-20% range.
We’ll have to run a similar workflow with the funcotator-annotated variants.

## 18.04.2025

### Review of PRNP mutation table

I’ve reviewed my summary table of mutations found by Mutect2 that carry the PASS tag (C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2024-10-19_read_count_summaries). 
None are within the PRNP protein coding region (except for E200K in the heterozygous sample).
Other than an intronic mutation in Ctrl5, no controls carry mutations. A point of interest and potential result could be that CJDs could have a higher mutation rate. Another finding pointing is this direction is that CJDs and ddPCR had a trend towards higher read counts supporting pathogenic mutations (was it only E200K?). This could also point towards a higher mutation rate in PRNP.
I also need to see whether there are more variants in TET2 and TTN. The funcotator output could be useful for this, as it might help me find protein coding regions mutations more easily.
I could modify the pipeline step upstream of the R script above (VCF to table) to include TET2 and TTN as well and inspect any variants manually.

### Using a panel of normals (PoN)

In Mutect2, the best practice is to use a PoN to filter artifacts more effectively. Prior to this project, I didn’t know whether the controls might contain somatic PRNP mutations or not. The results mentioned above (only 1 mutation in Ctrl5) suggest that this is not the case. Based on this finding, I can assume that they are indeed “normals”. Using a PoN would help me refine the PRNP mutation table, which are all intronic.
I can create a PoN with the controls (NOT with the commercial gDNA, because it comes from multiple donors).
Then I run Mutect2 with PoN mode on the CJDs, then FilterMutectCalls, then I review the table again.

## Create PoN

Chat on this part

### Sequence dictionary

Picard interval lists (used so that Mutect2 focuses on the target regions only) must reference a “.dict” file that matches your FASTA.

```
# Only run this once
gatk CreateSequenceDictionary \
  -R /home/mcarta/databases/chr2_chr4_chr20.fasta \
  -O /home/mcarta/databases/chr2_chr4_chr20.dict
```

### Write the three capture regions to a BED file

In the databases folder:

```
chr20   4686134   4701605
chr4   105233897 105237351
chr2   178739141 178741921   # note: start < end
```

### Convert BED → interval list

```
picard BedToIntervalList \
    I=targets.bed \
    O=/home/mcarta/databases/capture_targets.interval_list \
    SD=/home/mcarta/databases/chr2_chr4_chr20.dict
```

This created the file capture_targets.interval_list

### Create raw files for PoN

Code: /add/code/2025-04-18_create_pon/run_mutect2_controls_for_pon.sh
Output: /add/results/PoN/controls_tumor_only

### Merge VCFs

Merged six single‑sample Mutect2 VCFs (Ctrl1‑5 & Ctrl7) into one multi‑sample VCF suitable for CreateSomaticPanelOfNormals:

```
/add/code/2025-04-18_create_pon/merge_vcfs.sh
```

I also indexed the merged file with IndexFeatureFile.

### Created PoN

This PoN can be used to re-run the Mutect analysis of the CJDs.
Code: "/add/code/2025-04-18_create_pon/create_pon.sh"
Result: "/add/results/PoN/CJD_controls_PoN.vcf.gz"

### Mutect, processing and annotation

```
│/add/code/2025-04-18_mutect_with_pon/ contains 
│── 1_mutect2_cjd_with_pon.sh
│── 2_learn_orientation_model.sh
│── 3_FilterMutectCalls.sh
│── 4_normalise_vcfs.sh
└── 5_funcotator.sh
```

The intermediate files are saved in 

```
/add/results/mutect_output/2025-04-18_CJD_with_PoN/
```

Funcotator:

```
/add/results/annotation_output/funcotator_PoN/
```

## 19.04.2025

### All in one pipeline

Created all-in-one pipeline to analyse CJD, Ctrls and dilutions without PoN:

```
"/add/code/mutect_without_pon/run_no_pon_pipeline.sh"
```

Pipeline steps:
- Stage 0 : BAM indexing
- 1) Mutect2 tumour‑only  (no PoN)
- 2) LearnReadOrientationModel
- 3) FilterMutectCalls  (uses orientation priors)
- 4) bcftools normalise (split + left‑align)
- 5) Funcotator annotation

The results are stored in this directory tree:

```
/add/results/no_PoN/                     (rootDir)
  │─ mutect_raw/          (CJDs + Ctrls)        <-- stage 1
  │─ mutect_raw_dil/      (Dilutions)           <-- stage 1
  │─ orientation/                               <-- stage 2
  │─ orientation_dil/                           <-- stage 2
  │─ filtered/                                  <-- stage 3
  │─ filtered_dil/                              <-- stage 3
  │─ norm/                                      <-- stage 4
  │─ norm_dil/                                  <-- stage 4
  └─ annot/                                     <-- stage 5
      └─ dil/                                   <-- stage 5
```

While this runs, I’ll inspect the results of the PoN pipeline.

### VariantsToTable on PoN pipeline results

The results are saved in this dir:

```
/add/results/annotation_output/funcotator_PoN/
```

I’ll export these results as a table using this script:

```
"/add/code/2025-04-01_extract_filtered_prnp_variants/VariantsToTable.sh"
```

### No-PoN pipeline

I’ve re-run the CJDs, Ctrls and dilutions on this pipeline, which is the same as the PoN pipeline (i.e. it includes the orientation bias model, which I hadn’t included before).

## 21.04.2025

### Bcftools to annotate gnomAD variants

I copied the script from the PoN pipeline to the non-PoN code directory:
"/add/code/mutect_without_pon/run_no_pon_pipeline.sh"
The output is saved in
/add/results/no_PoN/annot_with_gnomAD/

## 26.04.2025

I’ve run the data from /add/results/no_PoN/annot_with_gnomAD/ through VariantsToTable:

```
/add/results/no_PoN/variant_tables/
```

Then I transferred the VCF to windows:

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-04-26_combined_variant_table_noPoN
```

And ran the R script to create a table. I noticed that the VCFs were huge, the reason being that they contained variants without the PASS tag (almost certainly artifacts). I need to apply a filtering step to remove non-PASS variants (I thought I already had this?). 

I added SelectVariants to the script. I also noticed that read counts weren’t being extracted, added that as well.
You’ll also have to look at the AAF and the gnomAD AAF to see whether the variants might be germline (rather than somatic).

## 29.04.2025

### VariantsToTable code

Look at the variants generated with the updated no-PoN pipeline (includes SelectVariants).
I saw that read counts were missing, so I edited the code:

```
"/add/code/mutect_without_pon/VariantsToTable_noPoN.sh"
```

Copied to directory with R formatting script:

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-04-26_combined_variant_table_noPoN
```

### Creating table of non-PoN variants in R (continued)

Edited the R script to include read counts. They are in the AD column, with this annotation:

```
AD: 20,15 → 20 reads support the REF allele, 15 support the first ALT allele
```

I edited the script to separate them into REF_count and ALT_count and calculate AAF.
I also added to code to extract various important data from the FUNCOTATION column.
- Gene symbol, gene region, sbSNP ID, mutation type

### Population frequencies

Extracting population frequency from gnomAD doesn’t seem to have worked, which is why I looked up the rs IDs manually on gnomad (https://gnomad.broadinstitute.org/) and dbSNP (there’s a link from the gnomAD page for every variant). I added the variants to the R code to be added manually. Some variants didn’t have rs IDs, there I looked the variants up manually (e.g. “2-178581749-C-T”) but didn’t find population data on any of these.

### Filtering parameters

I discovered that the brain somatic mosaicism network published an analysis best-practice workflow:
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02285-3?utm
It probably can’t be directly applied to my samples and I might have already implemented many of the recommendations. I want to add minimum values for DP, ALT_count and MQ to filter out junk.

BSMN “common pipeline” reference WGS (neurotypical cortex; 85–245 ×)
- Mean depth (brain): ≥ 50 × at the locus
- Keep a site only if
- ≥ 5 ALT reads and <3 in panel-of-normals
- VAF ≥ 1 %
- MQ ≥ 20; at least 1 read/strand; strand-bias P < 0.05

Ultra-deep epilepsy-gene panel
- Mean depth (brain): 2–5 k ×
- Keep a site only if
- ≥ 2 ALT reads per strand in two independent libraries (replicate rule)
The epilepsy methods can’t be directly applied to my data, as I didn’t use the two replicates per library strategy.

## 30.04.2025

### Issue with QUAL field

Mutect2 does not calculate the VCF QUAL field.
In raw Mutect2 output the column is written as “.” (missing). When you open the file with tools that expect a numeric value (e.g. bcftools, some R packages or Excel), the missing “.” is interpreted as -10 — a sentinel that means “no value”, not a real Phred score. Therefore -10 does not say anything about read depth or variant support (source)

### Implemented filtering for the non-PoN pipeline

In summary, variants were retained only if they fulfilled the following selection requirements, which include prior selection steps: (i) location within genomic regions covered by the hybridisation-capture probes used for sequencing; (ii) presence of a Mutect2 PASS tag, indicating a high-confidence somatic call; (iii) the variant was supported by ≥10 alternate reads, with  ≥3 reads on each DNA strand, at a total read depth of ≥100; (iv) the variant was either absent or recorded at an allele frequency of <0.001 in gnomAD and dbSNP; (v) the variant-allele fraction (VAF) differed significantly from 0.5 (consistent with a non-germline origin), according to a two-sided binomial test with p < 10⁻⁶; and (vi) VAF was >0.006, based on our limit-of-detection analysis.

### MQ filter

I also want to add base and mapping quality as they were used in the BSMN paper. 
(vi) mandate candidate mosaic SNVs have at least five independent non-duplicated supporting reads that have minimum values of 20 for mapping and base quality

### Re-run no-PoN pipeline

As these values were not all output by Mutect2, I had to edit the pipeline and re-run it. I also told the tools to only evaluate the intervals (i.e. the genomic regions targeted by our probes); this significantly sped up the pipeline (3h instead of 1d, or so)

## 01.05.2025

### Implement mapping and base quality filtering

I’ve created a new script:

```
"/add/code/mutect_without_pon/readcount_QC_pipeline.sh"
```

This script implements a streamlined pipeline to generate per-allele readcount metrics for variant positions in a no-PoN (Panel of Normals-free) analysis. It begins by extracting variant coordinates from gzipped VCF files into BED format using bcftools, preparing one BED file per sample or dilution. Next, recalibrated BAM files (which already contain proper read group and sample tags) are symlinked into a working directory to preserve the originals. PCR duplicates are then removed using samtools view with flag 0x400, and the resulting deduplicated BAMs are indexed. bam-readcount is then run on each deduplicated BAM at the corresponding variant positions, using a chromosome-limited FASTA reference. The output consists of tab-delimited readcount files, one per sample, which include per-base metrics such as counts, mapping quality, and base quality. Each step is built to skip existing output files, allowing the pipeline to resume cleanly if interrupted.

## Parse readcount file to readable TSV

AI failed to create a readable TSV from the readcount file using awk and python. Tomorrow, create a stand-alone python script with the objective of processing a single file successfully. If that can be done, implement a loop through all the files.

## 02.05.2025

### Python parsing to readable TSV, then R summary to CSV

Anthropic Claude did this easily. The code:

```
"/add/code/mutect_without_pon/readcount_to_TSV.py"
```

I’ve included it in the bash pipeline:

```
"/add/code/mutect_without_pon/readcount_QC_pipeline.sh"
```

Next, I moved the output to Windows:

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-05-02_readcount_TSV_noPoN
```

and created an R pipeline (readcount_dataframe.R) to summarise them to a large CSV each. I also added it to the bash pipeline.


### Integrate into R variant filtering pipeline

The relevant columns in the basecount summary CSV are:
- MEAN_BQ	Mean base-quality (Phred) of reads at this position.
- MEAN_MQ	Mean mapping-quality (Phred) of reads.
I left-joined these values with the summary_table rows, based on sample, genomic coordinates and variant and implemented the filter (removes 2 / 4 of the remaining variants, the other two were slightly below the MQ threshold of 20).

Next: Run this analysis on the dilutions, but without VAF limit.

### Where are the UMIs?

- I previously checked whether the reads in the R1 and R2 FASTQ files contain barcodes — they do not. Based on the o3 assessment, their length corresponds to to reads without barcodes
- This means that the barcodes were removed upon FASTQ creation. I should have received files containing the UMIs (I1 and I2), but I didn’t.
- Processing was for sureselect XT HS version 1 (not 2)
- I have contacted FGCZ to see whether the UMIs are either saved or can be recovered from the BCL files.

## 03.05.2025

### Search for UMIs

I found Stats_i1-8_i2-0.standard.json on the FGCZ server. (https://fgcz-gstore.uzh.ch/projects/p3111/, log in with b-fabric login). It seems like only the sample indexes were read, whereas the UMIs (i2) were not read. If this is confirmed, the UMIs weren’t read in the first place and can’t be recovered.
I’ll wait for the facility to answer. If it turns out that the sequencing settings were incorrect, I could ask for a free re-sequencing. Based on my notes (Excel sheets with library volume calculations), I might still have enough library left for a second run. I’d have to find the libraries — I recall putting them in tidy boxes in the freezer in the main lab (hopefully they haven’t been trashed).
In the meantime, I’ll continue work without UMI analysis.

### Examine result of QC pipeline on dilution samples

Goals: 
- demonstrate that the analysis and QC pipeline can detect the spike-in
- quantify the spike-in  define as limit of detection for the QC pipeline

The pipeline has already been run on the dilution samples. The respective files are stored in 

```
/add/results/no_PoN/
```

with the relevant directories ending in “_dil”.
A summary of the dilution readcounts was created yesterday:

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-05-02_readcount_TSV_noPoN_samples\output\readcounts_dilutions.csv
```

The TSV files (generated by the Mutect + annotation pipeline) are here:

```
/add/results/no_PoN/variant_tables/dilutions/ copied to C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-05-03_dilution_QC_analysis
```

I’ve created a table of the samples in which A117V was detected:

## 04.05.2025

### A117V table: R script and Latex

I made some changes to format the table nicely for publication, with the columns that are relevant. I tried to export the table with the gt package, but while the table was nice, the output formatting never worked properly. In the end, I created a Latex script that builds a booktabs style table with the values being defined explicitly (importing the values from the CSV didn’t work)

### Create PoN


- Code: /add/code/mutect_with_pon/PoN_creation/
- Result: /add/results/PoN/panel_of_normals/

### 
Run PoN pipeline

Re-analysed the CJD samples using the latest pipeline parameters. Adapted the no-PoN pipeline to use the newly created PoN. To speed it up, Mutect2, VariantAnnotator and Funcotator use the intervals tag to focus only on the 3 genes.

### Run QC pipeline

The QC pipeline probably doesn’t need editing.

## 10.05.2025

The PoN pipeline concluded successfully.

### VariantsToTable

"/add/code/mutect_with_pon/VariantsToTable_withPoN.sh"

### Readcounts to TSV

"/add/code/mutect_with_pon/readcount_to_TSV_CJD.py"

### Perform QC

C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-05-10_CJD_PoN_QC_analysis

For some reason, I don’t obtain base and mapping quality scores for all variants. I need to see whether the pipeline that produces these values processes all the variants I need, as some variants have scores, while others haven’t.
Actually, the issue might lie with the join, resulting in values being lost -> NA

## 20.05.2025

In createTable_CJD.R, located in 

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-05-10_CJD_PoN_QC_analysis
```

I’m working on the step where the BQ and MQ data are joined to my CJD variant data. For some reason, I get NA values for many variants, even though the dataset containing the quality data should be designed to contain values for all my variants of interest.
Solution: in the basecount dataset, some REF bases were in lower case, therefore the join didn’t work.
After applying QC, I am left with 4 variants. Only one is narrowly above the limit of detection of 0.8%, but others are around 0.5-0.6% and are backed by many ALT reads, so I think they can be considered legitimate.

## 02.08.2025

Create final table of SNVs for the paper

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Manuscript\Tables\SNV_summary_table
```

## 02.08.2025

### Create Lollipop plot

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Manuscript\Figures\SNV_lollipop
```

The design was polished in Inkscape (legend for samples, horizontal line).
I just noticed that three samples contain the same intronic mutation: 4691920 G>A. They were sequenced in the same batch, which could indicate contamination. Also, they are below my LOD of 0.8 (I hadn’t noticed this). Their low pravalence could suggest index hopping.
I think I will keep them in the paper, as they are supported by the number of reads in the best practice guidelines, but clearly state that they are less certain.
https://chatgpt.com/g/g-p-679f4d39752c81918971a06f142cf27e-cjd-mutation-analysis/project

Next:
- Make x-axis scale wider to cover entire gene track
- Re-generate ddPCR plots

## 18.09.2025

### Revise Lollipop plot

Extended x axis with genomic coordinates.
Tweaked image: improved visibility of elements, alignment. Image now complete.

### Review ddPCR analysis

- Find directories with data and understand the structure
- Review R script
- Are all the data being processed?
- Review the graph. Adding more data might disrupt formatting, if so, fix
- Edit in inkscape
- Also send to Mirka

#### Directory structure

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\ddPCR\ddPCR-files
```

Contains folders, one for each run (e.g. 2020-10-12), which in turn contain ddpcr files (QuantaSoft, contains graphs) and an output CSV.

#### CSV structure

https://chatgpt.com/c/68cbf5bb-8d88-8322-8fe7-9c85a6e1a7a3
I have checked that we have CSVs for every ddPCR
https://chatgpt.com/c/68cbfa12-f5d4-8324-a212-caa5c1c91fd3

#### Create big table of ddPCR data

Here:

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\ddPCR\CSV files pool
```

I’m reviewing the script create_snv_dataframe.R
It’s mostly fine, I’ve just tweaked some code to make it clearer. We might also need some extra QC analyses, which I’ve added on the end.

## 19.09.2025

### Continue with ddPCR dataframe script

I have reviewed the entire code. The next steps:
- Combine droplet counts for the samples with multiple analyses
- sum positives and sum accepted droplets across runs, then compute FA = Σ pos / Σ drops and derive Poisson CI from the pooled counts
- Define a QC strategy
- How to clearly distinguish from noise?

![Embedded image](analysis_journal_media/image8.png)

Next step: Continue reviewing the R script

## 25.09.2025

### ddPCR

- Completed code for pooling.
- Updated methods section with section on analysis.

### LoD, LoB, LoQ

- LoD: determined with my dilution experiment, already done
- LoB (limit of blank): calculate with blanks
- LoQ: I’m not sure I need this.

### LoB

- I’ve inserted code that calculates LoB for every sample based on the controls on the respective plate.
- Added methods section.

Positive samples: Samples exceeding both LoB and LoD.

### ddPCR Graphs

```
C:\...\Manuscript\Figures\ddPCR_fractional_abundance
```

In the dir above, there’s the R code to create the new graphs. The graphs need some tweaking as they look a bit messy.

## 26.09.2025

### Pooled droplets from whole brain

In the R script used for analysis yesterday:

```
create_snv_dataframe.R
```

I pooled the droplet counts for the whole brain and calculated LoB, based on the max p0 of the plates that were used (“max-of-plates” approach). Also created figures:

```
C:\...\Manuscript\Figures\ddPCR_fractional_abundance_pooled
```

The FAs are extremely low, most of them below LoB, all of them well below LoD.

### Limit of detection experiment

- Checked code for LoD experiment graphs
- Calculated LoB based on plate (the limitation is that it’s based on two wells only, but I calculated LoB based on the plate of the experiment for the samples). The value is not so different from the overall LoB (made from all control wells)
- Added LoB bar to LoD experiment result graph (“fraction)
- Value has to be above LoD and LoB. D178N lower CI intersects with LoB, but I think I’m fine. The detection limit is simply so low that it’s close to the level where noise predominates. Stress this point in the paper.

LoD Figures:

```
C:\...\Manuscript\Figures\ddPCR_LoD
```

LoB calculations for LoD experiment:

```
C:\...\Manuscript\Tables\LoD_calculations
```

To do: 
- In the fraction plots there’s no error bar for NTC. This is because the relevant values are NA. Investigate why. If you can obtain values, insert them. If not, remove NTC from the graph as it’s empty.

## 08.10.2025

### Complete ddPCR plots

I’ve added a “LoB” annotation to the horizontal line, tweaked the text sizes and used Helvetica font. The parameters can be adapted to the journal’s style.

## 09.10.2025

Currently going through the figures in the thesis, tweaking the new figures, writing legends. They are currently collected in a .tex

```
...\CJD PRNP\Manuscript\Figures\all figures with legends
```

## 20.10.2025

Today I’ll review the methods and results and make sure that they are publication ready.

### Table of ddPCR LoD results for the supplement

Code

```
…\Manuscript\Tables\LoD_calculations\LoD_table.R
```

Source material (equivalent for other mutations)

```
...\Experiments\ddPCR\LOD graphs\D178N\Manfredi_D178N_LOD.csv
```

Latex

```
LoD_table.tex
```

Latex with all tables

```
…\Manuscript\Tables\all tables with legends
```

## 21.10.2025

I have created an Excel of all ddPCR sample data for the supplement (Manuscript\Supplement)
Working on nice table of data for paper here:

```
…\Manuscript\Tables\ddPCR_sample_results\ddPCR_samples_results_tbl.R
```

## 10.11.2025

### Finish ddPCR tables

Excels with the results of the ddPCR analysis of CJDs and Controls were created with:

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\ddPCR\CSV files pool\create_snv_dataframe.R
```

The R creates a concise table for the paper that is saved as a string of Latex code in the R object latex_code

I copied the Latex code into:

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Manuscript\Tables\all tables with legends\all_tables.tex
```

The Latex table itself is huge and I’d include it in the Supplement. Maybe shade the LoD + LoB fulfilled cells blue to make them more visible.
Full Excel table for supplement is saved here:

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Manuscript\Supplementary data\ddPCR_results_table.xlsx
```

## 11.11.2025

### Table for sequencing limit of detection

I’ve integrated this table into the PDF with all tables and written the results section

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\SureSelect-sequencing\Analysis\2025-05-03_dilution_QC_analysis
```

## 12.11.2025

### CNV workflow

Brainstorming on which tools to use:
https://chatgpt.com/g/g-p-68e6cd7050ec81918984fe18cbd6143e-manfredi-s-research/c/69132c28-e3e0-832d-8b4a-c89cdaec18f7

### Setup

I’ve built a directory tree for this analysis in /add/prnp-cnv/
I’ve created a list of all samples to be processed in

```
/add/prnp-cnv/data/metadata/samples.txt
```

### Metric calculation calculation

```
/add/prnp-cnv/src/2025-11-12_GC_coverage_mapping.sh
```

- bedtools nuc: GC content
- bedcov: how many reads overlap each interval in the BED (i.e. the exons)
- bedcov: same for entire target

### Github

I created a repository called prnp-cnv, which I will use to sync the src folder. I’m not syncing the files as repository disk space might not be sufficient.

### 13.11.2025

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Experiments\prnp-cnv\src
```

In this exploratory analysis, I see that the PRNP exon has higher coverage in most samples, both CJD and controls. This is probably due to the denser probe coverage. As probe coverage creates an imbalance, I have to conclude that a CNV analysis is not feasible.

## 14.11.2025

I have abandoned the CNV analysis. Something I haven’t looked at yet is whether the probes targeting fused exons (targeting spliced PRNP, this idea was based on the spliced APP Nature paper) actually captured anything. I will see whether spliced PRNP exons can be found in my reads. This will be the last analysis.

## 15.11.2025 – 16.11.2025

### Exon-exon junction analysis

- Create short artificial contigs that encode the exon1→exon2 junction for my three PRNP transcripts, save as FASTA
- Set-up Github repo prnp-junctions
I've built the exon-exon junction FASTA 

```
…\prnp-junctions\src\01_build_prnp_junction_fasta.R
…\prnp-junctions\resources\prnp_junctions.fa
```

And a BED file for PRNP +/- 1kb:

```
…\src\02_make_prnp_bed.R
…\resources\PRNP.pad50kb.hg38.bed
```

Next: working on "/add/prnp-junctions/src/03_process_bam.sh"

## 17.11.2025

I have run 03_process_bam.sh:
- Input: Whole-panel BAM (CJD1.bam, Ctrl5.bam, …).
- Subset: Keep only reads overlapping PRNP (based on BED); remove duplicates; name-sort.
- Convert: Turn that subset into paired FASTQ.
- Realign: Align those reads to a custom PRNP junction FASTA.
- Output: A PRNP-specific, junction-centric BAM (${SAMPLE}.PRNP.toJunc.bam), one per sample.
Output is saved in:

```
/add/prnp-junctions/results/junction_align/
```

## 31.01.2026

I need to provide techinical details of the sequencing quality metrics in the paper supplement. I will place the relevant scripts in 

```
/add/authoritative_files/
```

On windows, authoritative references and input files are saved in

```
C:\Users\Manfredi\USZ\Neuropathologie - Carta Manfredi\CJD PRNP\Manuscript\Authoritative Reference Files
```

Authoritative_Refs_and_Inputs.docx: contains a list of the source and reference files.
manifest.tsv: List of samples, groups (dilution, CJD, control), batches (sequencing date), input_dir (directory containing the raw FASTQ and recalibrated BAMs)
seq_metrics_config.yaml: contains reference FASTA file, targets BED, filters, groups, naming conventions
sequencing_metrics_per_sample.schema.tsv: header string for QC file extraction

I double-checked that the files mentioned in manifest.tsv are complete with this script:

```
/add/code/2026-01-31_quality_metrics/validate_manifest.sh
```

The file that is output by the script (manifest_qc.tsv) confirms that all files are complete.

### Create quality metrics TSV

```
/add/authoritative_files/compute_sequencing_metrics.py
```

This script turns a manifest of samples into a schema-locked per-sample metrics table by counting reads with samtools view, computing target-depth distribution with samtools depth (MAPQ≥20; duplicates excluded) and extracting duplication/library size from Picard MarkDuplicates metrics.
Quality metrics were output here:

```
/add/results/2026-01-31_quality_metrics/sequencing_metrics_per_sample.tsv
```

Next: checked the TSV  in some samples, such as CJD5, the coverage looks poor. I checked the BAM in IGV and the issue is that the intron has poor coverage / mapping quality, whereas the PRNP Exon 2 is fine.

## 07.02.2026

As the intron often has poor coverage (especially in the 16-sample batch, though not all samples have poor coverage in the intron), but the ORF usually has decent to good coverage, I will extract QC metrics for the ORF, not only for the entire target regions.
Created BED with PRNP coding region: /home/mcarta/databases/prnp_coding.bed

```
chr20	4699220	4699982	PRNP_CODING_hg38_chr20_4699221_4699982
```

Updated locked schema file to add coding-region metrics: /add/authoritative_files/sequencing_metrics_per_sample.schema.tsv

Edited Python3 script to also compute coding-region metrics: /add/authoritative_files/compute_sequencing_metrics.py
There was an issue with mapping quality stats. I initially used samtools -q (yields base quality) instead of samtools -Q (mapping quality). Next, I’ll recreate the QC table using -Q.
The TSV shows metrics that are in line with what I see in IGV, which is good news.
The next step will be to inpect the TSV (probably in R) and create a table for the paper supplement. I will also write a results paragraph on the metrics (coverage etc).

## 08.02.2026

Created the repository skeleton for the PRNP somatic project with a clear separation between authoritative inputs, resources and results.
Defined the initial .gitignore and folder structure for reproducible analyses and future publication.
Legacy Ubuntu scripts were committed to
/add/prnp-somatic/src/legacy/ubuntu
as the historical baseline for the pipeline.
Added a Makefile target for QC metrics with logging and ensured that the metrics TSV is written directly to the outputs together with a stderr log for traceability.
Removed compiled artefacts and refreshed the inventory to keep the repository clean and deterministic.

## 09.02.2026

Added the authoritative files to the repository and enforced LF line endings via .gitattributes to avoid cross-platform inconsistencies between Windows and Ubuntu environments.
Pinned exact tool versions through the Conda environment and lockfile for publication-grade reproducibility.
Documented the toolchain and reference resources in /resources, together with SHA256 checksums.
Updated .gitignore to consistently exclude logs, config files and large resource files while keeping the structure of the resources directory under version control.
Harmonised the authoritative files with the legacy scripts so that the same inputs drive both historical and reproducible workflows.
Implemented a lightweight CI step to validate manifest and resource integrity.
Snapshotting of manifest.tsv and schema.tsv into run outputs with checksums was introduced to guarantee provenance of every analysis run.

## 15.02.2026

Reorganised the legacy Ubuntu scripts by moving them one directory level up to simplify path handling inside the Makefile and future pipeline scripts.
Updated .gitignore to exclude the fastq/ directory and the runs/ folder containing large intermediate outputs.
Created the preprocessing sample sheet
config/preprocessing_samples.tsv
containing batch membership and the paired FASTQ paths.
Added
config/preprocessing.env
for environment-specific settings and excluded it from version control.
Implemented
src/pipelines/preflight_preprocessing.sh
to validate that all required inputs, references and configuration variables are present before starting the preprocessing.
Wrote the first functional version of
src/pipelines/preprocessing.sh
and added the corresponding Makefile block and README documentation.
A test run of the preprocessing workflow is planned next.

## 16.02.2026

Refined preprocessing.sh (v1) and integrated it with the Makefile-based execution model so that the preprocessing can be launched reproducibly with a single target.
The preprocessing workflow now defines:
input discovery via the sample sheet
controlled output locations in runs/preprocessing
thread and memory configuration via the env file
The repository is now structured so that QC, preprocessing and variant calling can be executed as independent, resumable steps.

## 21.02.2026

The pre-processing step seems to have been incomplete. I am re-running it. It is important to run it within conda to ensure that all the correct tools and versions are used.
Afterwards, I’ll run mutect2 on the controls only:

```
/add/prnp-somatic/src/pipelines/1_controls_mutect2_no_pon.sh
```

The goal is to verify that the controls don't contain somatic PRNP mutations and make the PoN next.

## 22.02.2026

Major pipeline restructuring was completed today, with the repository moved to a staged, reproducible workflow for both no-PoN controls and PoN-based CJD+dilution analyses.

Summary of completed work (Git history + today’s execution):

- Standardised legacy path usage from `run/` to `runs/` and aligned defaults/config entries.
- Added `doc/sequencing_methods.md` with formatted methods text.
- Added controls no-PoN staged scripts and docs:
- `src/pipelines/1_controls_mutect2_no_pon.sh`
- `src/pipelines/2_controls_postprocess_no_pon.sh`
- `src/pipelines/3_controls_readcount_qc_no_pon.sh`
- `src/pipelines/4_readcount_to_tsv.py`
- `src/pipelines/5_controls_variant_qc_no_pon.sh`
- `src/pipelines/6_controls_variant_table_qc_no_pon.R`
- Added PoN creation and PoN run stages:
- `src/pipelines/7_controls_create_pon.sh`
- `src/pipelines/8_cjd_dilutions_mutect2_with_pon.sh`
- `src/pipelines/9_cjd_dilutions_postprocess_with_pon.sh`
- `src/pipelines/10_cjd_dilutions_readcount_qc_with_pon.sh`
- `src/pipelines/11_cjd_dilutions_readcount_to_tsv_with_pon.sh`
- `src/pipelines/12_cjd_dilutions_variant_qc_with_pon.sh`
- `src/pipelines/12_cjd_dilutions_variant_table_qc_with_pon.R`
- Updated repo and pipeline documentation (`README.md`, `src/pipelines/README.md`, `config/README.md`) to reflect staged execution and output layout.
- Added/updated config templates for new stage variables (`config/preprocessing.env.example`).
- Added Windows legacy scripts under `src/legacy/windows/` and a dedicated README.
- Added Funcotator resources under `resources/funcotator_data_somatic/` and updated ignore policy so large datasource artefacts remain out of Git.
- Verified that the repo-local Funcotator datasource path is structurally valid:
- `resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38`
- Ran preflight checks for controls post-processing prerequisites:
- tools (`gatk`, `bcftools`, `tabix`)
- reference FASTA + index
- stage-1 outputs (`.raw.vcf.gz`, `.stats`, `.f1r2.tar.gz`) for all controls
- Synced FASTQ storage between `/add/seq_data` and `/add/prnp-somatic/fastq`:
- identified and copied 6 missing files into `fastq/first_CJD_seq`
- verified all 78 mapped `.fastq.gz` files by binary equality (`cmp`)
- removed 78 redundant `.fastq.gz` files from `/add/seq_data` only after exact match confirmation

Outcome:

- Controls and PoN workflows are now represented as explicit reproducible stages in `src/pipelines`.
- Configuration and documentation are aligned with current staged execution.
- FASTQ storage is consolidated under `/add/prnp-somatic/fastq` with verified integrity.
- Conducted full pipeline rune

Next steps:

- Cross check results of the pipeline with tables and graphs that are already made.
- dilutions: remove the min AAF part from the QC, as you used the dilutions to determine the threshold
- add ddPCR to the repo
- add the exon-exon junction part to the repo
- add scripts creating figures to the repo
- sync the repo to Windows properly (in C/Projects)

