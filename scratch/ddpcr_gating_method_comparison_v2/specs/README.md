# v2 method specifications

Each implementation stage in this scratch experiment must have a reviewed
specification before code is added or extended for that stage.

The expected sequence is:

1. Draft or update the relevant `SPEC_*.md`.
2. Implement the stage in `code/`.
3. Run a method-specific smoke test.
4. Run the stage E2E on the intended scratch inputs.
5. Record results, failures, and fixes in `EXPERIMENTS.md`.

The specifications are based on installed package manuals/vignettes, published
method papers or reviews, and the raw-data contract of this repository.

## Specification Index

- `SPEC_00_shared_import_and_package_inputs.md`
- `SPEC_01_control_geometry.md`
- `SPEC_02_twoddpcr.md`
- `SPEC_03_ddPCRclust.md`
- `SPEC_04_dPCP_native.md`
- `SPEC_05_one_channel_methods.md`
- `SPEC_06_bayesian_mixture.md`
- `SPEC_07_lob_lod_reporting.md`
- `SPEC_08_polygon_gates.md`
