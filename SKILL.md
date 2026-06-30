---
name: genomic-baseline-browser
description: Launch an interactive genome browser for a GRCh38 VCF against a gnomAD + ClinVar baseline. Use this skill when the user gives a VCF and wants to browse coverage tracks, genes, ClinVar, and optionally compare their sample to gnomAD.
agent_created: true
---

# Genomic Baseline Browser

## Overview

Launch a local UCSC-style genome browser for a GRCh38 sample VCF. The browser shows pre-computed baseline tracks:

- **Genes** — Ensembl 115 GFF3 exons
- **gnomAD exome** — variant density from gnomAD v4.1 exomes
- **ClinVar** — pathogenic / VUS / benign stacked density

If a local gnomAD exome VCF directory is configured (`GNOMAD_EXOMES_VCF_DIR`), the sample is compared against the baseline and a per-sample variant track is added. If gnomAD is unavailable, a raw sample track is still built from the input VCF so the sample variants are visible.

## Input

- A GRCh38 sample VCF file (`.vcf` or `.vcf.gz`)
- Optional sample ID (inferred from VCF if omitted)

## Entry point

```bash
python scripts/run_skill.py /path/to/sample.vcf.gz
```

## Workflow

1. Detect sample ID from VCF header or filename.
2. Validate the VCF is GRCh38 (best-effort header check).
3. If `GNOMAD_EXOMES_VCF_DIR` points to real `gnomad.exomes.v4.1.sites.chr*.vcf.bgz` files:
   - Run `scripts/02_compare_sample.py` to create `compare/{sample_id}/compare_variants.parquet`.
   - Otherwise, build a raw `compare_variants.parquet` directly from the input VCF (all variants marked as novel).
4. Run `scripts/06_prepare_browser_data.py --skip-gff3` to refresh the browser manifest.
5. Start the FastAPI server on `http://127.0.0.1:8765/`.

## Environment variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `GNOMAD_EXOMES_VCF_DIR` | none | Path to `gnomad.exomes.v4.1.sites.chr*.vcf.bgz` for sample comparison |
| `GNOMAD_GENOMES_VCF_DIR` | none | Path to genome callset VCFs (if used) |
| `CLINVAR_VCF` | `~/.workbuddy/data/clinvar/clinvar.vcf.gz` | ClinVar GRCh38 VCF for live ClinVar track |
| `GFF3_PATH` | `raw/Homo_sapiens.GRCh38.115.gff3` | Ensembl GFF3 (only needed to rebuild genes) |
| `PORT` | `8765` | Browser API port |

## Dependencies

Install into the active Python environment:

```bash
pip install -r browser/requirements.txt
pip install -r scripts/requirements-pipeline.txt
```

Required packages: `fastapi`, `uvicorn`, `pandas`, `pyarrow`, `pysam`.

## Output

- Server URL: `http://127.0.0.1:8765/`
- Manifest API: `http://127.0.0.1:8765/api/manifest`
- Viewport API: `http://127.0.0.1:8765/api/chrom/{chrom}/viewport?start=N&end=N`

## Notes

- Pre-computed baseline bins are in `bins/` and `processed/browser/`. Do not delete them.
- Without `GNOMAD_EXOMES_VCF_DIR`, the sample track is built from the input VCF instead of being compared to gnomAD. Variants are shown as "novel"; gnomAD/ClinVar details for the variant detail pane are still queried on demand where data is available.
- PDF variant reports require the gnomAD VCFs and are generated on demand.
- The server runs in the foreground. Press `Ctrl+C` to stop.
