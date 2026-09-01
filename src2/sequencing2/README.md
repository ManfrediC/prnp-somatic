# Sequencing analytical pipeline (adapted draft)

This directory is a self-contained draft of the numbered analytical pipeline
in `src/pipelines`. It excludes FASTQ preprocessing and keeps the existing
inputs, outputs, sample selection, Panel of Normals construction, annotation,
read-count collection and variant-table quality-control behaviour.

## Adapted caller and filter settings

The control and CJD/dilution branches both use the transferable settings
selected in the spike-in pipeline:

```text
Mutect2:
--initial-tumor-lod 0
--max-reads-per-alignment-start 0

FilterMutectCalls:
--max-events-in-region 5
```

All other Mutect2 and FilterMutectCalls settings remain as they were in
`src/pipelines`. In particular, this analytical pipeline does not transfer
`--max-population-af 1.0`: the spike-in parameter review records that setting
as specific to spike-in validation and says the CJD/control analyses should
retain the existing population-AF limit.

The changes are applied consistently to the controls used to build the Panel
of Normals and to the CJD/dilution calls that use it:

- `1_controls_mutect2_no_pon.sh`
- `2_controls_postprocess_no_pon.sh`
- `8_cjd_dilutions_mutect2_with_pon.sh`
- `9_cjd_dilutions_postprocess_with_pon.sh`

The remaining numbered scripts and their Python and R helpers are copied so
this directory does not fall back to helper files in `src/pipelines`.

## Scope and paths

The scripts preserve the established environment-variable names and input
paths, but isolate every generated file from the original pipeline beneath
`results2/sequencing2/`:

- intermediate outputs: `results2/sequencing2/runs/`
- final quality-control tables and the new Panel of Normals:
  `results2/sequencing2/results/`

Each script reads `src2/sequencing2/sequencing2.env` when that optional file
exists. Copy `sequencing2.env.example` to that name to customise the draft.
Set `ENV_FILE` explicitly only when another configuration file is intended;
using the original `config/preprocessing.env` can restore its original output
paths and defeat this isolation. The scripts retain their existing resume/skip
behaviour within the selected output roots.

Create and activate the pipeline's dedicated Conda environment:

```bash
conda env create -f src2/sequencing2/environment.yml
conda activate prnp-sequencing2
```

Run scripts from the repository root through WSL. The existing numbered order
is unchanged: controls calling and post-processing, control quality control and
Panel of Normals construction, then CJD/dilution calling, post-processing,
read-count collection, parsing and final variant-table quality control.

For the publication-path Stage 12 wrapper, run:

```bash
bash src2/sequencing2/run_cjd_dilutions_variant_qc_with_pon.sh
```

This directory contains code only. No sequencing outputs were regenerated as
part of the adaptation.
