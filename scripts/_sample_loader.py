"""Load and normalize individual sample VCF records for baseline comparison."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterator

import pysam


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.startswith("chr"):
        return chrom
    if chrom in {"M", "MT"}:
        return "chrM"
    return f"chr{chrom}"


def _resolve_contig(vcf: pysam.VariantFile, chrom: str) -> str:
    candidates = [chrom, chrom.replace("chr", ""), normalize_chrom(chrom)]
    for name in candidates:
        if name in vcf.header.contigs:
            return name
    contigs = list(vcf.header.contigs)
    if not contigs:
        raise ValueError(f"no contigs in VCF header for {chrom}")
    return contigs[0]


def _parse_gt(record: pysam.VariantRecord, sample: str) -> tuple[str | None, str | None]:
    if sample not in record.samples:
        return None, None
    sample_data = record.samples[sample]
    gt = sample_data.get("GT")
    if gt is None or None in gt:
        return None, None
    gt_str = "/".join("." if allele is None else str(allele) for allele in gt)
    alleles = [a for a in gt if a is not None]
    if len(alleles) == 0:
        zygosity = None
    elif len(set(alleles)) == 1 and alleles[0] == 0:
        zygosity = "hom_ref"
    elif len(set(alleles)) == 1:
        zygosity = "hom_alt"
    else:
        zygosity = "het"
    return gt_str, zygosity


def iter_sample_variants(
    vcf_path: Path,
    sample_name: str | None = None,
) -> Iterator[dict[str, Any]]:
    """Yield one row per (record, alt) with normalized chrom and genotype fields."""
    with pysam.VariantFile(str(vcf_path)) as vcf:
        samples = list(vcf.header.samples)
        if not samples:
            raise ValueError(f"no samples in VCF: {vcf_path}")
        sample = sample_name or samples[0]

        for record in vcf:
            chrom = normalize_chrom(record.chrom)
            gt_str, zygosity = _parse_gt(record, sample)
            dp = None
            if sample in record.samples:
                dp_val = record.samples[sample].get("DP")
                if dp_val is not None:
                    dp = int(dp_val)

            for alt in record.alts or []:
                yield {
                    "chrom": chrom,
                    "pos": int(record.pos),
                    "ref": str(record.ref),
                    "alt": str(alt),
                    "variant_id": f"{chrom}-{record.pos}-{record.ref}-{alt}",
                    "sample_gt": gt_str,
                    "sample_zygosity": zygosity,
                    "sample_dp": dp,
                    "sample_qual": float(record.qual) if record.qual is not None else None,
                    "sample_id": sample,
                }


def load_sample_variants_by_chrom(
    vcf_path: Path,
    sample_name: str | None = None,
) -> dict[str, list[dict[str, Any]]]:
    by_chrom: dict[str, list[dict[str, Any]]] = {}
    for row in iter_sample_variants(vcf_path, sample_name=sample_name):
        by_chrom.setdefault(row["chrom"], []).append(row)
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda r: r["pos"])
    return by_chrom


def vcf_sample_name(vcf_path: Path) -> str:
    with pysam.VariantFile(str(vcf_path)) as vcf:
        samples = list(vcf.header.samples)
        if not samples:
            raise ValueError(f"no samples in VCF: {vcf_path}")
        return samples[0]
