# Backup resources

This directory stores large resource folders that are intentionally excluded from active pipeline execution but kept for traceability and potential restoration.

## Funcotator excluded datasources

The controls post-processing pipeline (`src/pipelines/2_controls_postprocess_no_pon.sh`) uses:

- `resources/funcotator_data_somatic/funcotator_dataSources.v1.8.hg38.20230908s/hg38`

`gnomAD_exome/` and `gnomAD_genome/` were moved out of that active datasource tree to avoid requester-pays access failures during Funcotator initialisation.

Current backup location:

- `resources/backup/funcotator_excluded_datasources/gnomAD_exome/`
- `resources/backup/funcotator_excluded_datasources/gnomAD_genome/`

If access constraints change in future, these directories can be moved back into the active `hg38/` datasource tree.
