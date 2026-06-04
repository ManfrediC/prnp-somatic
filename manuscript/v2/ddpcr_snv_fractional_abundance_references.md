# ddPCR SNV fractional abundance: reference rationale and manuscript wording

## Purpose

This note records the literature basis for the rebuilt ddPCR SNV fractional abundance calculation. The defensible claim is not that a single paper exactly reproduces the full PRNP somatic SNV workflow. Rather, the method is supported by combining:

- dMIQE-aligned dPCR Poisson quantification and reporting.
- Published allele-specific ddPCR mutation fractional abundance definitions.
- A published dPCR concentration-ratio confidence interval method.
- Empirical LoB/LoD handling for rare mutation detection.

## Recommended core references

### 1. dMIQE2020 guidelines

Huggett JF, et al. The Digital MIQE Guidelines Update: Minimum Information for Publication of Quantitative Digital PCR Experiments for 2020. Clinical Chemistry. 2020;66(8):1012-1029. doi: 10.1093/clinchem/hvaa125.

URL: https://academic.oup.com/clinchem/article/66/8/1012/5880117

Why include it:

- It is the community reporting standard for quantitative dPCR experiments.
- It defines the relevant dPCR terms: accepted partitions, positive partitions, negative partitions, Poisson correction, and lambda as the average number of target molecules per partition.
- It states that dPCR should first be analysed using lambda, and that when relative quantities such as fractional abundance are reported, lambda for the respective targets should also be reported.
- It explicitly discusses rare single-nucleotide variant assays, false-positive signal, assay specificity, thresholds, and the need for controls.

How it supports the manuscript:

- Supports the use of Poisson-corrected molecule occupancy as the measurement scale.
- Supports reporting target-level absolute quantities alongside fractional abundance.
- Supports explicit description of thresholds, droplet classification, and statistical methods.

### 2. Allele-specific ddPCR fractional abundance definition

Castellanos-Rizaldos E, Paweletz C, Song C, Oxnard GR, Mamon H, Janne PA, Makrigiorgos GM. Enhanced Ratio of Signals Enables Digital Mutation Scanning for Rare Allele Detection. The Journal of Molecular Diagnostics. 2015;17(3):284-292. doi: 10.1016/j.jmoldx.2014.12.003.

URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC4411249/

Why include it:

- It is an allele-specific ddPCR mutation paper using mutant and wild-type channels.
- It gives the exact fractional abundance definition needed here: mutant copies divided by mutant plus wild-type copies.
- It states that target copies are adjusted by the software to fit a Poisson distribution model.

How it supports the manuscript:

- Supports defining PRNP SNV fractional abundance as:

```text
FA = lambda_mut / (lambda_mut + lambda_wt)
```

- Supports the interpretation of the two-channel ddPCR assay as measuring mutant allele abundance relative to total target allele abundance.

### 3. dPCR concentration-ratio confidence interval

Dube S, Qin J, Ramakrishnan R. Mathematical Analysis of Copy Number Variation in a DNA Sample Using Digital PCR on a Nanofluidic Device. PLoS ONE. 2008;3(8):e2876. doi: 10.1371/journal.pone.0002876.

URL: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0002876

Why include it:

- The biological application in the paper is CNV, but the statistical problem is general dPCR ratio estimation.
- The paper derives concentration estimates from digital positive counts and gives confidence intervals for the ratio of two true target concentrations.
- The ratio problem maps directly to allele-specific SNV ddPCR:

```text
CNV ratio: lambda_target / lambda_reference
SNV ratio: lambda_mut / lambda_wt
FA transform: ratio / (1 + ratio)
```

How it supports the manuscript:

- Supports the confidence interval around the mutant:wild-type concentration ratio.
- Supports transforming the concentration ratio CI to an FA CI.

Important wording constraint:

- Do not describe this as a somatic SNV-specific method.
- Describe it as a published dPCR concentration-ratio method applied to the mutant:wild-type concentration ratio.

### 4. Mutation assay LoB/LoD

Milbury CA, Zhong Q, Lin J, Williams M, Olson J, Link DR, Hutchison B. Determining lower limits of detection of digital PCR assays for cancer-related gene mutations. Biomolecular Detection and Quantification. 2014;1(1):8-22. doi: 10.1016/j.bdq.2014.08.001.

URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC5129438/

Why include it:

- It is directly about digital PCR mutation detection.
- It separates molecule counting from the analytical sensitivity problem.
- It discusses LoB, LoD, false-positive events, false-positive rates, mutant:wild-type ratios, and low-frequency cancer mutation assays.

How it supports the manuscript:

- Supports keeping empirical LoB/LoD as the detection decision layer.
- Supports the argument that rare SNV positivity should not rely on the partition-count confidence interval alone.
- Supports explicit treatment of false-positive mutant droplets and assay background.

## Optional supporting reference

Pablo-Fontecha V, Hernandez-Illan E, Reparaz A, et al. Quantification of rare somatic single nucleotide variants by droplet digital PCR using SuperSelective primers. Scientific Reports. 2023;13:18997. doi: 10.1038/s41598-023-39874-0.

URL: https://www.nature.com/articles/s41598-023-39874-0

Why include it:

- It is explicitly about rare somatic SNV quantification by ddPCR.
- It uses variant allele frequency as the biological quantity and compares ddPCR VAF with sequencing-derived VAF.
- It supports the broader idea that ddPCR can quantify rare somatic SNVs, although its assay design differs from the two-channel mutant/wild-type hydrolysis-probe setup used here.

Suggested use:

- Useful in the introduction or methods discussion if the manuscript needs an explicitly somatic-SNV ddPCR citation.
- Less central than Castellanos-Rizaldos et al. for the exact mutant plus wild-type fractional abundance formula.

## LoB and p0 rule: evidence status

### 1. Detection capability framework

CLSI. Evaluation of Detection Capability for Clinical Laboratory Measurement Procedures; Approved Guideline--Second Edition. CLSI document EP17-A2. Wayne, PA: Clinical and Laboratory Standards Institute; 2012.

URL: https://clsi.org/media/1430/ep17a2_sample.pdf

Armbruster DA, Pry T. Limit of blank, limit of detection and limit of quantitation. Clinical Biochemistry Reviews. 2008;29 Suppl 1:S49-S52.

URL: https://pubmed.ncbi.nlm.nih.gov/18852857/

Why include them:

- CLSI EP17-A2 is the general clinical-laboratory framework for LoB, LoD, and LoQ.
- Armbruster and Pry give a concise, citable explanation of the same framework.
- Both support treating LoB as a high quantile of the blank/no-analyte distribution, not as a biological signal estimate.

How it supports the current code:

- The code's `qbinom(0.95, size = n_tot, prob = p0_use)` is an implementation of the same principle on the droplet-count scale: it asks how many mutant-positive droplets would still be plausible background noise in a blank reaction with the same number of accepted droplets.
- The code's `detected_LoB = x_mut > LoB_count` is therefore a count-scale decision: observed mutant-positive droplets must exceed the blank distribution's 95th percentile.

What it does not prove:

- EP17-A2 and Armbruster/Pry do not prescribe this exact ddPCR implementation.
- They support the concept of a 95th-percentile blank threshold. The exact binomial model, the use of a Clopper-Pearson upper bound for p0, and the max-of-plates rule are conservative implementation choices.

### 2. Exact binomial p0 upper bound

Clopper CJ, Pearson ES. The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika. 1934;26(4):404-413. doi: 10.1093/biomet/26.4.404.

URL: https://doi.org/10.1093/biomet/26.4.404

Why include it:

- The code estimates the blank false-positive probability as the upper 95% exact binomial confidence bound on pooled blank positives.
- This is deliberately conservative: it uses the upper plausible value for the blank-positive probability rather than the observed point estimate.

How it supports the current code:

```text
p0_use = upper 95% exact binomial confidence bound for blank positives / blank droplets
LoB_count = 95th percentile of Binomial(n_sample_droplets, p0_use)
```

Interpretation:

- This is defensible when blank mutant-positive droplets are treated as Bernoulli/binomial background events after droplet classification.
- It is not a named standard ddPCR LoB recipe found in a single paper. It is a transparent conservative construction from standard LoB logic plus exact binomial uncertainty.

### 3. ddPCR rare-mutation background and controls

Milbury CA, Zhong Q, Lin J, Williams M, Olson J, Link DR, Hutchison B. Determining lower limits of detection of digital PCR assays for cancer-related gene mutations. Biomolecular Detection and Quantification. 2014;1(1):8-22. doi: 10.1016/j.bdq.2014.08.001.

URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC5129438/

Why include it:

- It is directly about dPCR mutation assays and low-frequency mutation detection.
- It treats false-positive mutation events as a limiting factor for analytical sensitivity.
- It supports separating molecule counting from the detection decision at very low mutant fractions.

How it supports the current code:

- The PRNP workflow estimates fractional abundance from Poisson-corrected molecule occupancy, then applies empirical LoB/LoD criteria as a separate detection layer.
- This is the right structure for rare-variant ddPCR: a non-zero FA estimate is not automatically a mutation-positive call.

Control-material remark:

- dMIQE2020 explicitly distinguishes no-template negative controls from controls containing the same nucleic-acid background without the targeted variant. For rare variants in mostly wild-type samples, same-background WT controls are the more relevant estimate of assay specificity and low-level false-positive mutant signal.
- NTCs are still useful for contamination monitoring. They are less biologically matched to the test samples because they do not contain wild-type genomic background.
- Therefore the current LoB p0 model uses WT controls only. The previous WT+NTC p0 model is retained only as a sensitivity/QC comparison, because WT controls are the better biological blank for rare SNV false positives.

## Method rationale

The current calculation is literature-backed if described as a modular, dMIQE-aligned dPCR analysis:

1. Accepted droplets are classified by the gating strategy into mutant-positive, wild-type-positive, both-positive, or negative partitions.
2. For each target, the target molecule occupancy is estimated using the standard dPCR Poisson correction:

```text
lambda = -log(negative_partitions / accepted_partitions)
```

3. Mutant fractional abundance is calculated from Poisson-corrected mutant and wild-type target occupancies:

```text
FA = 100 * lambda_mut / (lambda_mut + lambda_wt)
```

4. Confidence intervals are calculated for the mutant:wild-type concentration ratio using a published dPCR concentration-ratio approach and then transformed to fractional abundance:

```text
ratio = lambda_mut / lambda_wt
FA = 100 * ratio / (1 + ratio)
```

5. Detection calls are made using empirical assay-specific LoB/LoD criteria, not by the FA confidence interval alone.

This is scientifically stronger than the earlier droplet-level binomial approach because the earlier approach estimated the fraction of mutant-positive droplets among accepted droplets, whereas the biological question is the fraction of target molecules carrying the mutant allele.

## Suggested manuscript text

> Droplet classifications were reconstructed from the accepted QuantaSoft gate-derived partition clusters. For each target, molecule occupancy per partition was estimated by standard digital PCR Poisson correction, lambda = -log(w/n), where w is the number of target-negative accepted partitions and n is the number of accepted partitions. Mutant fractional abundance was calculated as 100 x lambda_mut/(lambda_mut + lambda_wt), consistent with published allele-specific ddPCR mutation analyses in which mutant copies are expressed relative to mutant plus wild-type copies. Confidence intervals were calculated for the mutant:wild-type concentration ratio using a published digital PCR concentration-ratio approach and transformed to fractional abundance. Mutation positivity was not assigned from the fractional abundance interval alone, but from assay-specific empirical LoB/LoD criteria to account for false-positive mutant partitions and assay background. In line with dMIQE2020 recommendations, target-level lambda values and accepted partition counts were retained for audit and supplementary reporting.

## Suggested citation placement

- Cite dMIQE2020 after the Poisson/lambda reporting sentence.
- Cite Castellanos-Rizaldos et al. after the mutant fractional abundance formula.
- Cite Dube et al. after the concentration-ratio confidence interval sentence.
- Cite Milbury et al. after the LoB/LoD detection sentence.
- Cite Pablo-Fontecha et al. only if an explicitly rare somatic SNV ddPCR reference is useful for context.

## Claims to avoid

- Do not claim that Dube et al. is a somatic SNV ddPCR paper.
- Do not claim that the FA confidence interval alone validates rare mutation detection.
- Do not claim that the method removes the need for empirical false-positive, LoB, or LoD assessment.
- Do not claim that QuantaSoft's internal method is identical to Dube/Fieller; our test implementation gives the same point estimates and unchanged LoB/LoD outcomes, but confidence intervals are not bit-for-bit identical.
