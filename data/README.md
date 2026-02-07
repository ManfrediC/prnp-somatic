# --- Data policy README ---
cat > data/README.md <<'EOF'
# Data policy

This repository does not track raw sequencing data or large derived files.

Expected local layout (example):
- `data/fastq/`   raw FASTQs (not tracked)
- `data/bam/`     aligned BAM (not tracked)
- `data/vcf/`     VCF outputs (not tracked)
- `data/raw/`     any other raw inputs (not tracked)

Provide templates in `config/` (manifest, config) so others can reproduce the workflow.
