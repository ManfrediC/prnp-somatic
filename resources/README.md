# Resources

This directory contains workflow reference files. Small project-specific files are
committed; large public references and generated indexes stay local and are
reconstructed or downloaded as described in
`reproducibility/tooling_and_reference_provenance.md`.

Use `resources/INDEX.tsv` as the canonical file inventory.

## Folder Map

- `adapters/`: small adapter FASTA used by preprocessing.
- `annotations/`: small manual annotation tables.
- `funcotator/`: local Funcotator datasource tree.
- `intervals/`: capture and coding BED/interval files.
- `junctions/`: committed junction reference FASTA/BED plus local BWA indexes.
- `known_sites/`: local dbSNP and Mills/1000G VCF resources.
- `population/`: local gnomAD AF-only and optional 1000G PoN resources.
- `references/hg38/`: local whole-genome hg38 FASTA/GTF plus small sidecars.
- `references/snv/`: local chr2/chr4/chr20 subset FASTA plus small sidecars.
- `repeats/`: committed PRNP ORR repeat workflow resources.

## Git Policy

- Tracked: small hand-authored or project-specific resources, documentation,
  adapter FASTA, BED/interval files, repeat catalogues, `.dict`, and `.fai`.
- Local ignored: large FASTA, GTF, VCF, VCF index, and Funcotator datasource
  files.
- Generated ignored: BWA sidecars (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`).

## Verification

From the repository root:

```bash
make verify_resources
```

This checks `resources/SHA256SUMS.txt` against the files present locally.
