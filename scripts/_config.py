"""Paths and constants for genome-wide exome baseline.

Environment overrides:
  GNOMAD_EXOMES_VCF_DIR  -- directory containing gnomad.exomes.v4.1.sites.chr*.vcf.bgz
  GNOMAD_GENOMES_VCF_DIR -- directory containing gnomad.genomes.v4.1.sites.chr*.vcf.bgz
  GNOMAD_EXOMES_TSV_DIR  -- directory containing gnomAD exomes TSVs
  GNOMAD_GENOMES_TSV_DIR -- directory containing gnomAD genomes TSVs
  CLINVAR_VCF            -- path to ClinVar GRCh38 VCF.gz
  GFF3_PATH              -- path to Ensembl GFF3 for gene re-extraction
  CONSTRAINT_TSV         -- path to gnomAD constraint metrics TSV
"""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BINS = ROOT / "bins"
BY_CHROM = BINS / "by_chrom"
COMPARE = ROOT / "compare"
SCRIPTS = ROOT / "scripts"
FIGURES = ROOT / "figures"
RAW = ROOT / "raw"
SAMPLES = RAW / "samples"
CYTOBAND = RAW / "cytoBand.txt"
BROWSER = ROOT / "processed" / "browser"
BROWSER_GENES = BROWSER / "genes"


def _env_path(name: str, default: str | Path) -> Path:
    value = os.environ.get(name)
    if value:
        return Path(value)
    return Path(default) if isinstance(default, str) else default


GFF3 = _env_path(
    "GFF3_PATH",
    ROOT / "raw" / "Homo_sapiens.GRCh38.115.gff3",
)

# gnomAD v4.1 release paths (read-only).  If the local VCFs are not present,
# sample comparison and variant report PDF generation gracefully degrade.
GNOMAD_ROOT = _env_path("GNOMAD_ROOT", "/dev/null")
GNOMAD_EXOMES_VCF_DIR = _env_path(
    "GNOMAD_EXOMES_VCF_DIR",
    GNOMAD_ROOT / "vcf" / "exomes",
)
GNOMAD_GENOMES_VCF_DIR = _env_path(
    "GNOMAD_GENOMES_VCF_DIR",
    GNOMAD_ROOT / "vcf" / "genomes",
)
GNOMAD_EXOMES_TSV_DIR = _env_path(
    "GNOMAD_EXOMES_TSV_DIR",
    GNOMAD_ROOT / "tsv" / "tsv" / "exomes",
)
GNOMAD_GENOMES_TSV_DIR = _env_path(
    "GNOMAD_GENOMES_TSV_DIR",
    GNOMAD_ROOT / "tsv" / "tsv" / "genomes",
)

# Default baseline callset (WES); use GNOMAD_GENOMES_* for WGS
GNOMAD_DIR = GNOMAD_EXOMES_VCF_DIR

CLINVAR_VCF = _env_path(
    "CLINVAR_VCF",
    Path.home() / ".workbuddy" / "data" / "clinvar" / "clinvar.vcf.gz",
)
CONSTRAINT_TSV = _env_path(
    "CONSTRAINT_TSV",
    ROOT / "raw" / "gnomad.v4.1.constraint_metrics.tsv",
)

REPORT_SCHEMA_VERSION = "1.1"
GNOMAD_POPULATIONS: tuple[str, ...] = (
    "afr",
    "amr",
    "asj",
    "eas",
    "fin",
    "mid",
    "nfe",
    "sas",
)

BIN_RESOLUTIONS: tuple[int, ...] = (1_000_000, 100_000, 10_000)
RESOLUTION_LABELS: dict[int, str] = {
    1_000_000: "1mb",
    100_000: "100kb",
    10_000: "10kb",
}

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
ALL_CHROMSOMES: tuple[str, ...] = (*AUTOSOMES, "chrX", "chrY")

DATA_RELEASE = {
    "reference_genome": "GRCh38",
    "gnomad": "v4.1",
    "gnomad_callset": "exomes",
    "clinvar": "GRCh38 VCF",
}


def gnomad_vcf_path(chrom: str, callset: str = "exomes") -> Path:
    """Sites VCF for exomes (default) or genomes callset."""
    if callset == "genomes":
        return GNOMAD_GENOMES_VCF_DIR / f"gnomad.genomes.v4.1.sites.{chrom}.vcf.bgz"
    if callset != "exomes":
        raise ValueError(f"unsupported gnomAD callset: {callset!r}")
    return GNOMAD_EXOMES_VCF_DIR / f"gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz"


def gnomad_tsv_dir(callset: str = "exomes") -> Path:
    if callset == "genomes":
        return GNOMAD_GENOMES_TSV_DIR
    if callset == "exomes":
        return GNOMAD_EXOMES_TSV_DIR
    raise ValueError(f"unsupported gnomAD callset: {callset!r}")


def clinvar_contig(chrom: str) -> str:
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    return chrom


def resolution_label(resolution: int) -> str:
    return RESOLUTION_LABELS[resolution]


def gff_contig(chrom: str) -> str:
    if chrom.startswith("chr"):
        return chrom[3:]
    return chrom


# UCSC-style zoom ladder (bp window); first step resolved per-chrom at runtime.
ZOOM_STEPS_BP: tuple[int, ...] = (10_000_000, 1_000_000, 100_000, 10_000)
CLICK_JUMP_BP = 100_000
DEFAULT_BROWSER_CHROM = "chr21"
DEFAULT_BROWSER_SAMPLE = ""
