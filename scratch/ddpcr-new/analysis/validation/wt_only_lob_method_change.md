# WT-only LoB method change

## Summary

The main ddPCR SNV analysis now estimates the Limit of Blank (LoB) background rate from `WT_control` wells only. `NTC` wells remain useful contamination controls, but they are no longer pooled into the LoB blank model.

Checkpoint before this change:

```text
b216ee7 Checkpoint ddPCR analysis scripts before WT-only LoB
```

## Code Changed

File:

```text
scratch/ddpcr-new/analysis/scripts/create_snv_dataframe_corrected.R
```

### Before

```r
# --------------------------------------------------------
# >>> >>> LIMIT OF BLANK (LoB) SECTION <<< <<<
# - Use WT + NTC controls from the same plate (Date) and assay
# - Conservative p0: CP upper 95% bound on pooled blank proportion
# - Fallback: assay-wide p0 if a plate lacks blanks
# --------------------------------------------------------

# 1) Build blank table (QC: >=10,000 droplets) from your controls object
blanks <- data.mut.controls %>%
  filter(Target %in% mutation.list, Sample %in% c("WT_control","NTC")) %>% # only blank/control wells
  filter(AcceptedDroplets >= 10000) %>% # only wells with at least 10,000 accepted droplets
  transmute(plate = Date,
            assay = ExperimentType,
            n = AcceptedDroplets,
            x = Positives)
```

### After

```r
# --------------------------------------------------------
# >>> >>> LIMIT OF BLANK (LoB) SECTION <<< <<<
# - Use WT genomic controls from the same plate (Date) and assay
# - NTCs are contamination controls; they are not used to model rare SNV
#   background in wild-type genomic DNA
# - Conservative p0: CP upper 95% bound on pooled blank proportion
# - Fallback: assay-wide p0 if a plate lacks blanks
# --------------------------------------------------------

# WT controls are the biological blanks for the LoB calculation because they
# contain genomic DNA background without the targeted mutation.
lob_blank_samples <- "WT_control"

# 1) Build blank table (QC: >=10,000 droplets) from WT control wells
blanks <- data.mut.controls %>%
  filter(Target %in% mutation.list, Sample %in% lob_blank_samples) %>% # WT genomic blank wells
  filter(AcceptedDroplets >= 10000) %>% # only wells with at least 10,000 accepted droplets
  transmute(plate = Date,
            assay = ExperimentType,
            n = AcceptedDroplets,
            x = Positives)
```

## What Did Not Change

The downstream LoB calculation is unchanged:

```r
LoB_count = qbinom(0.95, size = n_tot, prob = p0_use)
detected_LoB = x_mut > LoB_count
```

The same WT-only `blank_pooled` and `assay_fallback` tables feed both:

- sample-region LoB calls
- participant-pooled LoB calls

The max-of-plates rule is unchanged. If a sample or participant pool includes droplets from multiple plates, the analysis still uses the highest applicable plate-level blank p0.

## Reason

The LoB is meant to model false-positive mutant signal in reactions that contain the same background as the tested biological samples, but not the targeted mutation. `WT_control` wells satisfy that condition because they contain wild-type genomic DNA.

`NTC` wells answer a different QC question: whether there is template contamination or reagent-level background in a no-template reaction. They do not contain wild-type genomic DNA, so they are not the best biological blank for rare SNV false-positive signal.

The literature backing is therefore stronger with WT-only LoB:

- CLSI/Armbruster-Pry support LoB as a high quantile of the blank/no-analyte distribution.
- Clopper-Pearson supports the exact binomial upper confidence bound used for conservative p0 estimation.
- dMIQE2020 and rare-mutation ddPCR literature support using appropriate negative controls and separating analytical detection thresholds from raw fractional abundance estimates.

## Diagnostic Result

The sensitivity diagnostic compares the previous WT+NTC LoB model with the new WT-only LoB model:

```text
scratch/ddpcr-new/analysis/scripts/check_wt_only_lob_sensitivity.R
```

WT-only LoB made every LoB count threshold higher in this dataset, but did not change combined LoD+LoB positive calls:

```text
sample-region rows:          669
sample-region LoB+:          198 -> 38
sample-region LoD+LoB+:        9 -> 9
changed LoD+LoB+ calls:        0

participant-pooled rows:     116
participant LoB+:             31 -> 2
participant LoD+LoB+:          1 -> 1
changed LoD+LoB+ calls:        0
```

After this code change, those same WT-only thresholds are the main analysis output. The WT+NTC results remain available only through the sensitivity diagnostic.
