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

The final variant-table quality-control stages use the empirical
complete-pipeline LoD derived by the spike-in workflow as their AAF threshold:
`26 / 3891 = 0.006682086867129272`. The existing strict `AAF > threshold`
comparison is unchanged.

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

Run the full analytical pipeline from the repository root through WSL:

```bash
bash src2/sequencing2/run_sequencing2.sh
```

The driver follows the original documented order: controls calling and
post-processing, control quality control and Panel of Normals construction,
then CJD/dilution calling, post-processing, read-count collection, parsing and
final variant-table quality control. It excludes FASTQ preprocessing, as does
the rest of this adapted analytical directory.

After a complete verified run, archive the intermediate `runs/` tree with
Windows-native file handling:

```powershell
powershell -NoProfile -File src2/sequencing2/archive_sequencing2_intermediates.ps1 `
  -RemoveLocal -MakeOnlineOnly
```

The archival script requires the control, Panel of Normals, CJD and dilution
final outputs before it starts. It copies the intermediates to a new run-specific
folder under `CJD intermediates`, excludes the regenerable `bam_work` symlinks,
verifies every copied file by size and SHA-256, then removes the local `runs/`
tree and asks the cloud provider to make the archive online-only. Original BAMs
under `results/final_bam/` and all final outputs under
`results2/sequencing2/results/` remain on device. No live symlink is created,
because it could rehydrate archived files and interfere with resume checks.

For the publication-path Stage 12 wrapper, run:

```bash
bash src2/sequencing2/run_cjd_dilutions_variant_qc_with_pon.sh
```

This directory contains code only. Generated outputs are written under
`results2/sequencing2/`.
