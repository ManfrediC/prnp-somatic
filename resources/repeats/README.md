# Repeat Resources

This directory contains committed small reference assets for repeat-analysis
workflows in this repository.

## Files

- `prnp_orr.expansionhunter.json`
  ExpansionHunter variant catalog for the PRNP octapeptide repeat region.

## Notes

- The PRNP ORR catalog models the native locus with a fixed `R1` prefix, a
  variable block of `R2`-like repeats, and fixed `R3`/`R4` suffix repeats.
- This representation is intended for targeted short-read screening of OPRI and
  OPRD candidates from the existing BAMs in `results/final_bam`.
