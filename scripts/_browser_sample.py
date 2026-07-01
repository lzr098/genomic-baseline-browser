"""Sample variant binned track for the interactive browser (compare parquet)."""

from __future__ import annotations

import json
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

import pysam

from _bin_builder import bin_key, load_clinvar_map
from _browser_gff3 import gene_context_for_position
from _config import COMPARE, ROOT, gnomad_vcf_path
from _variant_report import (
    _allele_index,
    _gnomad_af_for_allele,
    _info_allele_value,
    _resolve_gnomad_contig,
    gnomad_variant_page_url,
    pick_vep_for_alt,
)

_VEP_DISPLAY_LABELS: dict[str, str] = {
    "5_prime_UTR_variant": "5' UTR",
    "3_prime_UTR_variant": "3' UTR",
    "upstream_gene_variant": "upstream",
    "downstream_gene_variant": "downstream",
    "intron_variant": "intron",
    "synonymous_variant": "synonymous",
    "missense_variant": "missense",
    "stop_gained": "stop gained",
    "stop_lost": "stop lost",
    "frameshift_variant": "frameshift",
    "splice_donor_variant": "splice donor",
    "splice_acceptor_variant": "splice acceptor",
    "splice_donor_5th_base_variant": "splice donor",
    "splice_region_variant": "splice region",
    "non_coding_transcript_exon_variant": "nc exon",
    "non_coding_transcript_variant": "nc transcript",
}

_GERMLINE_LABELS: dict[str, str] = {
    "pathogenic": "Pathogenic",
    "likely_pathogenic": "Likely pathogenic",
    "vus": "Uncertain significance",
    "conflicting": "Conflicting interpretations",
    "benign": "Benign",
    "likely_benign": "Likely benign",
}


def compare_variants_path(sample_id: str) -> Path:
    return COMPARE / sample_id / "compare_variants.parquet"


def compare_summary_path(sample_id: str) -> Path:
    return COMPARE / sample_id / "compare_summary.json"


def list_browser_samples(root: Path = ROOT) -> list[dict[str, Any]]:
    """Discover samples with precomputed compare_variants.parquet."""
    samples: list[dict[str, Any]] = []
    compare_dir = root / "compare"
    if not compare_dir.is_dir():
        return samples

    for entry in sorted(compare_dir.iterdir()):
        if not entry.is_dir():
            continue
        parquet = entry / "compare_variants.parquet"
        if not parquet.exists():
            continue
        summary: dict[str, Any] = {}
        summary_path = entry / "compare_summary.json"
        if summary_path.exists():
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        counts = summary.get("counts", {})
        samples.append(
            {
                "id": entry.name,
                "label": entry.name,
                "variant_count": int(counts.get("total_variants", 0)),
                "novel_count": int(counts.get("novel_variants", 0)),
                "compare_parquet": str(parquet.relative_to(root)),
            }
        )
    return samples


@lru_cache(maxsize=4)
def _load_sample_frame(sample_id: str, root_s: str) -> pd.DataFrame:
    path = Path(root_s) / compare_variants_path(sample_id)
    if not path.exists():
        raise FileNotFoundError(f"missing compare parquet for sample {sample_id}: {path}")
    return pd.read_parquet(path)


def load_sample_track_bins(
    sample_id: str,
    chrom: str,
    start: int,
    end: int,
    resolution: int,
    root: Path = ROOT,
) -> dict[str, Any]:
    """Stacked bins: known (gnomAD-matched) + novel counts per resolution window."""
    frame = _load_sample_frame(sample_id, str(root))
    sub = frame[(frame["chrom"] == chrom) & (frame["pos"] >= start) & (frame["pos"] <= end)]

    bin_counters: dict[int, dict[str, int]] = defaultdict(lambda: {"known": 0, "novel": 0})
    for _, row in sub.iterrows():
        bk = bin_key(int(row["pos"]), resolution)
        bucket = "novel" if bool(row.get("is_novel")) else "known"
        bin_counters[bk][bucket] += 1

    bins: list[dict[str, int]] = []
    for bk in sorted(bin_counters.keys()):
        counts = bin_counters[bk]
        known = int(counts["known"])
        novel = int(counts["novel"])
        total = known + novel
        if total <= 0:
            continue
        bins.append(
            {
                "start": int(max(bk, start)),
                "end": int(min(bk + resolution - 1, end)),
                "known": known,
                "novel": novel,
                "total": total,
            }
        )

    novel_n = int(sub["is_novel"].sum()) if not sub.empty else 0
    return {
        "sample_id": sample_id,
        "variant_count": int(len(sub)),
        "novel_count": novel_n,
        "bins": bins,
    }


def _canonical_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom.replace('chr', '')}-{pos}-{ref}-{alt}"


def _format_vep_annotation(consequence: str | None) -> str | None:
    if not consequence:
        return None
    primary = consequence.split("&")[0].strip()
    if primary in _VEP_DISPLAY_LABELS:
        return _VEP_DISPLAY_LABELS[primary]
    return primary.replace("_variant", "").replace("_", " ")


def _vep_category(consequence: str | None) -> str:
    if not consequence:
        return "other"
    primary = consequence.split("&")[0].strip()
    if primary == "missense_variant":
        return "missense"
    if primary in {"5_prime_UTR_variant", "3_prime_UTR_variant"}:
        return "utr"
    if primary in {
        "stop_gained",
        "stop_lost",
        "frameshift_variant",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "splice_donor_5th_base_variant",
        "splice_region_variant",
    }:
        return "lof"
    return "other"


def _format_hgvs_consequence(hgvsc: str | None, hgvsp: str | None) -> str | None:
    if hgvsp:
        return hgvsp.split(":", 1)[-1] if ":" in hgvsp else hgvsp
    if hgvsc:
        return hgvsc.split(":", 1)[-1] if ":" in hgvsc else hgvsc
    return None


def _format_germline_label(
    classification: str | None,
    clinical_tier: str | None = None,
) -> str | None:
    if classification:
        normalized = classification.replace("_", " ").strip()
        lower = normalized.lower()
        if "pathogenic" in lower and "likely" in lower:
            return "Pathogenic/Likely pathogenic"
        if "uncertain" in lower:
            return "Uncertain significance"
        if "conflicting" in lower:
            return "Conflicting interpretations"
        if lower.startswith("pathogenic"):
            return "Pathogenic"
        if lower.startswith("likely benign"):
            return "Likely benign"
        if lower.startswith("benign"):
            return "Benign"
        return normalized
    if clinical_tier and clinical_tier in _GERMLINE_LABELS:
        return _GERMLINE_LABELS[clinical_tier]
    return None


def _gnomad_region_index(
    chrom: str,
    start: int,
    end: int,
) -> dict[tuple[int, str, str], dict[str, Any]]:
    path = gnomad_vcf_path(chrom)
    if not path.exists():
        return {}

    index: dict[tuple[int, str, str], dict[str, Any]] = {}
    with pysam.VariantFile(str(path)) as vcf:
        contig = _resolve_gnomad_contig(vcf, chrom)
        for rec in vcf.fetch(contig, start - 1, end):
            if rec.pos < start or rec.pos > end:
                continue
            info = dict(rec.info)
            for alt in rec.alts or []:
                allele_idx = _allele_index(rec, alt)
                ac = _info_allele_value(info, "AC", allele_idx)
                an = _info_allele_value(info, "AN", allele_idx)
                af = _info_allele_value(info, "AF", allele_idx)
                if af is None:
                    af = _gnomad_af_for_allele(rec, alt)
                vep_strings = list(info.get("vep") or [])
                annotation = pick_vep_for_alt(vep_strings, alt)
                index[(rec.pos, rec.ref, alt)] = {
                    "allele_count": int(ac) if ac is not None else None,
                    "allele_number": int(an) if an is not None else None,
                    "allele_frequency": float(af) if af is not None else None,
                    "homozygote_count": int(nhomalt) if (nhomalt := _info_allele_value(info, "nhomalt", allele_idx)) is not None else None,
                    "hgvsc": annotation.get("hgvsc") or None,
                    "hgvsp": annotation.get("hgvsp") or None,
                    "consequence": annotation.get("consequence") or None,
                }
    return index


def _gnomad_genomes_presence(
    chrom: str,
    start: int,
    end: int,
) -> set[tuple[int, str, str]]:
    path = gnomad_vcf_path(chrom, callset="genomes")
    if not path.exists():
        return set()

    found: set[tuple[int, str, str]] = set()
    with pysam.VariantFile(str(path)) as vcf:
        contig = _resolve_gnomad_contig(vcf, chrom)
        for rec in vcf.fetch(contig, start - 1, end):
            if rec.pos < start or rec.pos > end:
                continue
            for alt in rec.alts or []:
                found.add((rec.pos, rec.ref, alt))
    return found


def _germline_classification(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    clinvar_map: dict[tuple[int, str, str], dict[str, Any]],
    fallback_clinsig: str | None,
    fallback_tier: str | None,
) -> str | None:
    clin = clinvar_map.get((pos, ref, alt))
    if clin and clin.get("clinsig"):
        return str(clin["clinsig"])
    if fallback_clinsig:
        return str(fallback_clinsig)
    if fallback_tier and fallback_tier != "none":
        return str(fallback_tier)
    return None


def _build_variant_rows(
    rows_df: pd.DataFrame,
    chrom: str,
    root: Path,
) -> list[dict[str, Any]]:
    """Build detailed variant rows from a filtered DataFrame subset."""
    if rows_df.empty:
        return []

    pos_min = int(rows_df["pos"].min())
    pos_max = int(rows_df["pos"].max())
    gnomad_index = _gnomad_region_index(chrom, pos_min, pos_max)
    genomes_present = _gnomad_genomes_presence(chrom, pos_min, pos_max)
    clinvar_map = load_clinvar_map(chrom)

    variants: list[dict[str, Any]] = []
    for _, row in rows_df.iterrows():
        pos = int(row["pos"])
        ref = str(row["ref"])
        alt = str(row["alt"])
        key = (pos, ref, alt)
        gnomad = gnomad_index.get(key, {})
        hgvsc = gnomad.get("hgvsc")
        hgvsp = gnomad.get("hgvsp")
        consequence = gnomad.get("consequence")
        variant_id = _canonical_variant_id(chrom, pos, ref, alt)
        in_exome = key in gnomad_index
        in_genome = key in genomes_present
        raw_clinsig = _germline_classification(
            chrom,
            pos,
            ref,
            alt,
            clinvar_map,
            row.get("clinsig"),
            row.get("clinical_tier"),
        )
        gene_context = gene_context_for_position(chrom, pos, root)
        location_consequence = _location_consequence_from_gene_context(gene_context)
        variants.append(
            {
                "variant_id": variant_id,
                "variant_page": gnomad_variant_page_url(variant_id),
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "locus": f"{chrom}:{pos:,}",
                "is_novel": bool(row.get("is_novel")),
                "match_status": str(row.get("match_status") or ""),
                "source_exome": in_exome,
                "source_genome": in_genome,
                "hgvs_consequence": _format_hgvs_consequence(hgvsc, hgvsp),
                "vep_annotation": _format_vep_annotation(consequence),
                "vep_category": _vep_category(consequence),
                "location_consequence": location_consequence,
                "lof_curation": None,
                "germline_classification": _format_germline_label(
                    raw_clinsig,
                    str(row.get("clinical_tier") or "") or None,
                ),
                "flags": ["novel"] if bool(row.get("is_novel")) else [],
                "gene_context": gene_context,
                "allele_count": gnomad.get("allele_count"),
                "allele_number": gnomad.get("allele_number"),
                "allele_frequency": gnomad.get("allele_frequency")
                if gnomad.get("allele_frequency") is not None
                else (float(row["gnomad_af"]) if pd.notna(row.get("gnomad_af")) else None),
                "homozygote_count": gnomad.get("homozygote_count"),
            }
        )

    return variants


def _location_consequence_from_gene_context(gene_context: dict[str, Any] | None) -> str | None:
    """Return a simple consequence string from gene context when VEP is missing."""
    if not gene_context:
        return "intergenic"
    location = (gene_context.get("location") or "").lower()
    if "exon" in location:
        return "exonic"
    if "intron" in location:
        return "intronic"
    if "flanking" in location:
        return "flanking"
    return "genic"


def load_sample_bin_variant_details(
    sample_id: str,
    chrom: str,
    bin_start: int,
    bin_end: int,
    root: Path = ROOT,
    *,
    limit: int = 100,
) -> dict[str, Any]:
    """Return per-variant browser detail rows for one sample histogram bin."""
    frame = _load_sample_frame(sample_id, str(root))
    sub = frame[
        (frame["chrom"] == chrom)
        & (frame["pos"] >= int(bin_start))
        & (frame["pos"] <= int(bin_end))
    ].sort_values(["pos", "ref", "alt"])

    total_count = int(len(sub))
    if total_count == 0:
        return {
            "sample_id": sample_id,
            "chrom": chrom,
            "bin_start": int(bin_start),
            "bin_end": int(bin_end),
            "total_count": 0,
            "returned_count": 0,
            "truncated": False,
            "variants": [],
        }

    rows = sub.head(int(limit))
    variants = _build_variant_rows(rows, chrom, root)

    return {
        "sample_id": sample_id,
        "chrom": chrom,
        "bin_start": int(bin_start),
        "bin_end": int(bin_end),
        "total_count": total_count,
        "returned_count": len(variants),
        "truncated": total_count > len(variants),
        "variants": variants,
    }


def load_sample_viewport_variants(
    sample_id: str,
    chrom: str,
    start: int,
    end: int,
    root: Path = ROOT,
    *,
    limit: int = 500,
) -> dict[str, Any]:
    """Return per-variant browser detail rows across the full current viewport."""
    frame = _load_sample_frame(sample_id, str(root))
    sub = frame[
        (frame["chrom"] == chrom)
        & (frame["pos"] >= int(start))
        & (frame["pos"] <= int(end))
    ].sort_values(["pos", "ref", "alt"])

    total_count = int(len(sub))
    if total_count == 0:
        return {
            "sample_id": sample_id,
            "chrom": chrom,
            "start": int(start),
            "end": int(end),
            "total_count": 0,
            "returned_count": 0,
            "truncated": False,
            "variants": [],
        }

    rows = sub.head(int(limit))
    variants = _build_variant_rows(rows, chrom, root)

    return {
        "sample_id": sample_id,
        "chrom": chrom,
        "start": int(start),
        "end": int(end),
        "total_count": total_count,
        "returned_count": len(variants),
        "truncated": total_count > len(variants),
        "variants": variants,
    }
