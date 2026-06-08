# Vendored Umbrella Code

This directory contains the upstream 2D Umbrella R scripts used for optional
local smoke testing of the PRNP ddPCR Umbrella input exports.

Source repository: https://github.com/statOmics/umbrella

Vendored from upstream `master` commit:

`afc58372485e94cbda7783b112778faf948aeda2`

Vendored files:

- `2D/Umbrella_2d_V1.R`
- `2D/Graphics_Umbrella_2d_V0.R`

The conversion scripts in `src/ddpcr/umbrella/` create input objects for these
functions. They do not modify the vendored implementation.
