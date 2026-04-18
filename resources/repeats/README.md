# Repeat Resources

This directory contains committed small reference assets for repeat-analysis
workflows in this repository.

## Files

- `prnp_orr.expansionhunter.json`
  ExpansionHunter variant catalog for the PRNP octapeptide repeat region.
- `prnp_orr.gangstr.bed`
  GangSTR target-locus BED for the mutable middle PRNP ORR repeat block.

## Notes

- The PRNP ORR catalog models the native locus with a fixed `R1` prefix, a
  variable block of `R2`-like repeats, and fixed `R3`/`R4` suffix repeats.
- This representation is intended for targeted short-read screening of OPRI and
  OPRD candidates from the existing BAMs in `results/final_bam`.
- The GangSTR BED targets only the mutable `R2` block, which is `2` reference
  copies of the canonical `24 bp` repeat motif in the hg38 PRNP ORR sequence.
