"""Build single-variant JSON reports and render PDF from gnomAD + ClinVar + sample VCF."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from io import BytesIO
from pathlib import Path
from typing import Any

import pysam

from _bin_builder import clinical_tier
from _compare import classify_match, priority_score
from _config import (
    CLINVAR_VCF,
    CONSTRAINT_TSV,
    DATA_RELEASE,
    GNOMAD_EXOMES_TSV_DIR,
    GNOMAD_EXOMES_VCF_DIR,
    GNOMAD_GENOMES_TSV_DIR,
    GNOMAD_GENOMES_VCF_DIR,
    GNOMAD_POPULATIONS,
    REPORT_SCHEMA_VERSION,
    clinvar_contig,
    gnomad_vcf_path,
)
from _sample_loader import _resolve_contig, normalize_chrom

VEP_KEYS = (
    "allele",
    "consequence",
    "impact",
    "symbol",
    "gene_id",
    "feature_type",
    "feature",
    "biotype",
    "exon",
    "intron",
    "hgvsc",
    "hgvsp",
)

VERDICT_TEXT = {
    "known_gnomad_clinvar": "Known in gnomAD exomes and ClinVar",
    "known_gnomad": "Known in gnomAD exomes",
    "known_clinvar_only": "Known in ClinVar only",
    "novel_in_sample": "Novel relative to gnomAD exomes and ClinVar",
}

PDF_CLINICAL_COLORS = {
    "pathogenic": (198, 40, 40),
    "likely_pathogenic": (249, 168, 37),
    "vus": (105, 126, 21),
    "conflicting": (249, 168, 37),
    "benign": (124, 179, 66),
    "likely_benign": (124, 179, 66),
    "none": (117, 117, 117),
    "other": (117, 117, 117),
}


def first_value(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def as_str_list(value: Any) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, tuple):
        return [str(item) for item in value]
    return [str(value)]


def parse_vep_entry(vep_str: str) -> dict[str, str]:
    parts = vep_str.split("|")
    return {key: parts[i] if i < len(parts) else "" for i, key in enumerate(VEP_KEYS)}


def pick_canonical_vep(vep_entries: list[dict[str, str]], gene_hint: str | None = None) -> dict[str, str]:
    if not vep_entries:
        return {}
    if gene_hint:
        for entry in vep_entries:
            if entry.get("symbol") == gene_hint:
                return entry
    for entry in vep_entries:
        if entry.get("impact") in {"HIGH", "MODERATE"}:
            return entry
    return vep_entries[0]


def pick_vep_for_alt(vep_strings: list[str], alt: str) -> dict[str, str]:
    entries = [parse_vep_entry(item) for item in vep_strings]
    if not entries:
        return {}
    for entry in entries:
        if entry.get("allele") == alt:
            return entry
    # gnomAD may annotate insertions/deletions with surrogate alleles.
    surrogates = {alt}
    if len(alt) > len("A") and len(alt) > 1:
        surrogates.add(alt[-1])
        surrogates.add("-")
    for entry in entries:
        if entry.get("allele") in surrogates:
            return entry
    return pick_canonical_vep(entries)


def parse_gnomad_variant_id(variant_id: str) -> tuple[str, int, str, str]:
    """Parse gnomAD variant ID: 1-55051215-G-GA -> chr1, pos, ref, alt."""
    variant_id = variant_id.strip()
    parts = variant_id.split("-")
    if len(parts) < 4:
        raise ValueError(f"invalid gnomAD variant_id: {variant_id!r}")
    chrom_token, pos_s, ref = parts[0], parts[1], parts[2]
    alt = "-".join(parts[3:])
    chrom = normalize_chrom(chrom_token)
    return chrom, int(pos_s), ref, alt


def gnomad_variant_page_url(variant_id: str) -> str:
    return f"https://gnomad.broadinstitute.org/variant/{variant_id}?dataset=gnomad_r4"


def _allele_index(rec: pysam.VariantRecord, alt: str) -> int:
    alts = list(rec.alts or [])
    return alts.index(alt)


def _info_allele_value(info: dict[str, Any], key: str, allele_idx: int) -> Any:
    value = info.get(key)
    if value is None:
        return None
    if isinstance(value, tuple):
        return value[allele_idx] if allele_idx < len(value) else value[0]
    return value


def priority_tier(score: float) -> str:
    if score >= 80:
        return "high"
    if score >= 30:
        return "elevated"
    return "routine"


def report_stem(chrom: str, pos: int, ref: str, alt: str) -> str:
    chrom_slug = chrom.replace("chr", "")
    return f"{chrom_slug}-{pos}-{ref}-{alt}"


def should_render_pdf(comparison: dict[str, Any], clinical_tier: str | None) -> bool:
    if comparison.get("is_novel"):
        return True
    if clinical_tier in {"pathogenic", "likely_pathogenic", "vus"}:
        return True
    if float(comparison.get("priority_score") or 0) >= 30:
        return True
    return False


def _resolve_gnomad_contig(vcf: pysam.VariantFile, chrom: str) -> str:
    for name in (chrom, chrom.replace("chr", "")):
        if name in vcf.header.contigs:
            return name
    contigs = list(vcf.header.contigs)
    if not contigs:
        raise ValueError(f"no contigs in gnomAD VCF for {chrom}")
    return contigs[0]


def _gnomad_af_for_allele(rec: pysam.VariantRecord, alt: str) -> float | None:
    alts = rec.alts or []
    if alt not in alts:
        return None
    af = rec.info.get("AF")
    if af is None:
        return None
    if isinstance(af, tuple):
        idx = alts.index(alt)
        return float(af[idx]) if idx < len(af) else float(af[0])
    return float(af)


def _safe_float(value: str | None) -> float | None:
    if value is None or value == "" or value.upper() == "NA":
        return None
    return float(value)


def _load_constraint_by_gene() -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    if not CONSTRAINT_TSV.exists():
        return out
    with CONSTRAINT_TSV.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("canonical") != "true":
                continue
            gene = row.get("gene")
            if not gene or gene in out:
                continue
            out[gene] = {
                "gene": gene,
                "transcript": row.get("transcript"),
                "lof_oe": _safe_float(row.get("lof.oe")),
                "pli": _safe_float(row.get("lof.pLI")),
                "mis_z": _safe_float(row.get("mis.z_score")),
            }
    return out


_CONSTRAINT_CACHE: dict[str, dict[str, Any]] | None = None


def get_gene_constraint(gene_symbol: str | None) -> dict[str, Any] | None:
    global _CONSTRAINT_CACHE
    if not gene_symbol:
        return None
    if _CONSTRAINT_CACHE is None:
        _CONSTRAINT_CACHE = _load_constraint_by_gene()
    return _CONSTRAINT_CACHE.get(gene_symbol)


def fetch_gnomad_variant(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    callset: str = "exomes",
) -> dict[str, Any] | None:
    gnomad_path = gnomad_vcf_path(chrom, callset=callset)
    if not gnomad_path.exists():
        raise FileNotFoundError(f"missing gnomAD VCF: {gnomad_path}")

    with pysam.VariantFile(str(gnomad_path)) as vcf:
        contig = _resolve_gnomad_contig(vcf, chrom)
        for rec in vcf.fetch(contig, pos - 1, pos):
            if rec.pos != pos or rec.ref != ref or alt not in (rec.alts or []):
                continue
            info = dict(rec.info)
            allele_idx = _allele_index(rec, alt)
            populations: dict[str, dict[str, Any]] = {}
            for pop in GNOMAD_POPULATIONS:
                ac = _info_allele_value(info, f"AC_{pop}", allele_idx)
                af = _info_allele_value(info, f"AF_{pop}", allele_idx)
                an = _info_allele_value(info, f"AN_{pop}", allele_idx)
                if ac is None and af is None:
                    continue
                populations[pop.upper()] = {
                    "ac": int(ac) if ac is not None else None,
                    "an": int(an) if an is not None else None,
                    "af": float(af) if af is not None else None,
                }

            vep_strings = list(info.get("vep") or [])
            annotation = pick_vep_for_alt(vep_strings, alt)
            gene_symbol = annotation.get("symbol") or None
            transcript_consequences = []
            seen = set()
            for item in vep_strings:
                entry = parse_vep_entry(item)
                key = (entry.get("feature"), entry.get("consequence"))
                if key in seen:
                    continue
                seen.add(key)
                if entry.get("symbol") or entry.get("feature"):
                    transcript_consequences.append(
                        {
                            "allele": entry.get("allele"),
                            "gene_symbol": entry.get("symbol") or None,
                            "transcript": entry.get("feature") or None,
                            "consequence": entry.get("consequence") or None,
                            "impact": entry.get("impact") or None,
                            "hgvsc": entry.get("hgvsc") or None,
                            "hgvsp": entry.get("hgvsp") or None,
                        }
                    )

            ac = _info_allele_value(info, "AC", allele_idx)
            an = _info_allele_value(info, "AN", allele_idx)
            nhomalt = _info_allele_value(info, "nhomalt", allele_idx)

            return {
                "callset": f"{callset}_v4.1",
                "filter": list(rec.filter),
                "overall": {
                    "ac": int(ac) if ac is not None else 0,
                    "an": int(an) if an is not None else 0,
                    "af": float(_info_allele_value(info, "AF", allele_idx))
                    if _info_allele_value(info, "AF", allele_idx) is not None
                    else _gnomad_af_for_allele(rec, alt),
                    "nhomalt": int(nhomalt) if nhomalt is not None else 0,
                },
                "populations": populations,
                "non_ukb": {
                    "ac": int(v) if (v := _info_allele_value(info, "AC_non_ukb", allele_idx)) is not None else None,
                    "an": int(v) if (v := _info_allele_value(info, "AN_non_ukb", allele_idx)) is not None else None,
                    "af": float(v) if (v := _info_allele_value(info, "AF_non_ukb", allele_idx)) is not None else None,
                },
                "grpmax": str(first_value(info.get("grpmax"))) if info.get("grpmax") is not None else None,
                "predictions": {
                    "cadd_phred": float(first_value(info["cadd_phred"]))
                    if info.get("cadd_phred") is not None
                    else None,
                    "revel_max": float(first_value(info["revel_max"]))
                    if info.get("revel_max") is not None
                    else None,
                    "spliceai_ds_max": float(first_value(info["spliceai_ds_max"]))
                    if info.get("spliceai_ds_max") is not None
                    else None,
                    "phylop": float(first_value(info["phylop"])) if info.get("phylop") is not None else None,
                    "sift_max": float(first_value(info["sift_max"])) if info.get("sift_max") is not None else None,
                    "polyphen_max": float(first_value(info["polyphen_max"]))
                    if info.get("polyphen_max") is not None
                    else None,
                },
                "annotation": {
                    "canonical_transcript": annotation.get("feature") or None,
                    "gene_symbol": gene_symbol,
                    "consequence": annotation.get("consequence") or None,
                    "impact": annotation.get("impact") or None,
                    "hgvsc": annotation.get("hgvsc") or None,
                    "hgvsp": annotation.get("hgvsp") or None,
                    "exon": annotation.get("exon") or None,
                    "intron": annotation.get("intron") or None,
                },
                "transcript_consequences": transcript_consequences[:12],
                "gene_constraint": get_gene_constraint(gene_symbol),
                "source_vcf": str(gnomad_path),
            }
    return None


def fetch_clinvar_variant(chrom: str, pos: int, ref: str, alt: str) -> dict[str, Any] | None:
    if not CLINVAR_VCF.exists():
        return None
    contig = clinvar_contig(chrom)
    with pysam.VariantFile(str(CLINVAR_VCF)) as vcf:
        if contig not in vcf.header.contigs:
            return None
        for rec in vcf.fetch(contig, pos - 1, pos):
            if rec.pos != pos or rec.ref != ref or alt not in (rec.alts or []):
                continue
            info = dict(rec.info)
            clinsig_raw = info.get("CLNSIG")
            if isinstance(clinsig_raw, tuple):
                clinsig_raw = ",".join(str(x) for x in clinsig_raw)
            elif clinsig_raw is not None:
                clinsig_raw = str(clinsig_raw)
            return {
                "clinsig": as_str_list(info.get("CLNSIG")),
                "clinical_tier": clinical_tier(clinsig_raw),
                "review_status": as_str_list(info.get("CLNREVSTAT")),
                "condition": as_str_list(info.get("CLNDN")),
                "hgvs": first_value(info.get("CLNHGVS")),
                "gene_info": first_value(info.get("GENEINFO")),
                "molecular_consequence": as_str_list(info.get("MC")),
                "variant_type": first_value(info.get("CLNVC")),
                "submission": as_str_list(info.get("CLNSIGSCV")),
                "origin": as_str_list(info.get("ORIGIN")),
                "rsid": first_value(info.get("RS")),
                "source_vcf": str(CLINVAR_VCF),
            }
    return None


def lookup_sample_variant(
    sample_vcf: Path,
    sample_id: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> dict[str, Any]:
    chrom = normalize_chrom(chrom)
    with pysam.VariantFile(str(sample_vcf)) as vcf:
        contig = _resolve_contig(vcf, chrom)
        for rec in vcf.fetch(contig, pos - 1, pos):
            if rec.pos != pos or rec.ref != ref or alt not in (rec.alts or []):
                continue
            gt_str, zygosity = None, None
            dp = None
            ad_ref = ad_alt = None
            if sample_id in rec.samples:
                sample_data = rec.samples[sample_id]
                gt = sample_data.get("GT")
                if gt is not None and None not in gt:
                    gt_str = "/".join(str(a) for a in gt)
                    alleles = list(gt)
                    if len(set(alleles)) == 1 and alleles[0] == 0:
                        zygosity = "homozygous_ref"
                    elif len(set(alleles)) == 1:
                        zygosity = "homozygous_alt"
                    else:
                        zygosity = "heterozygous"
                dp_val = sample_data.get("DP")
                if dp_val is not None:
                    dp = int(dp_val)
                ad = sample_data.get("AD")
                if ad is not None and len(ad) >= 2:
                    ad_ref, ad_alt = int(ad[0]), int(ad[1])
            allele_depth = None
            if ad_ref is not None and ad_alt is not None:
                allele_depth = {"ref": ad_ref, "alt": ad_alt}
            filt = list(rec.filter) if rec.filter else None
            return {
                "sample_id": sample_id,
                "genotype": gt_str,
                "zygosity": zygosity,
                "depth": dp,
                "allele_depth": allele_depth,
                "quality": float(rec.qual) if rec.qual is not None else None,
                "filter": ",".join(filt) if filt else None,
                "source_vcf": str(sample_vcf),
            }
    return {
        "sample_id": sample_id,
        "genotype": None,
        "zygosity": None,
        "depth": None,
        "allele_depth": None,
        "quality": None,
        "source_vcf": str(sample_vcf),
    }


def build_comparison_block(
    in_gnomad_exomes: bool,
    in_clinvar: bool,
    clinical_tier_value: str,
    priority_score_value: float,
    *,
    in_gnomad_genomes: bool = False,
) -> dict[str, Any]:
    in_gnomad = in_gnomad_exomes
    match_status = classify_match(in_gnomad, in_clinvar)
    return {
        "match_status": match_status,
        "is_novel": match_status == "novel_in_sample",
        "in_gnomad": in_gnomad,
        "in_gnomad_exomes": in_gnomad_exomes,
        "in_gnomad_genomes": in_gnomad_genomes,
        "in_clinvar": in_clinvar,
        "verdict": VERDICT_TEXT[match_status],
        "priority_score": priority_score_value,
        "priority_tier": priority_tier(priority_score_value),
        "baseline_scope": "gnomad_exomes_and_clinvar",
    }


def _format_af(af: float | None) -> str:
    if af is None:
        return "not observed"
    if af == 0:
        return "0"
    if af < 0.001:
        return f"{af:.2e}"
    return f"{af:.6f}"


def _format_af_display(af: float | None) -> str:
    if af is None:
        return "n/a"
    if af < 0.001:
        return f"{_format_af(af)} ({af * 100:.4f}%)"
    return f"{_format_af(af)} ({af * 100:.2f}%)"


def _format_num(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if value == 0:
            return "0"
        if value.is_integer():
            return f"{int(value):,}"
        if abs(value) < 0.001:
            return f"{value:.2e}"
        return f"{value:.4f}"
    return str(value)


def _format_display(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return _format_num(value)
    return str(value)


def build_interpretation(
    sample: dict[str, Any],
    comparison: dict[str, Any],
    gnomad: dict[str, Any] | None,
    clinvar: dict[str, Any] | None,
) -> dict[str, Any]:
    flags: list[str] = []
    gene = (gnomad or {}).get("annotation", {}).get("gene_symbol")
    consequence = (gnomad or {}).get("annotation", {}).get("consequence") or ""
    af = (gnomad or {}).get("overall", {}).get("af")
    ac = (gnomad or {}).get("overall", {}).get("ac")
    an = (gnomad or {}).get("overall", {}).get("an")
    tier = (clinvar or {}).get("clinical_tier") or "none"

    if comparison.get("is_novel"):
        flags.append("novel_in_sample")
    if af is not None and af < 0.001:
        flags.append("low_frequency")
    if tier in {"pathogenic", "likely_pathogenic"}:
        flags.append("clinvar_pathogenic")
    elif tier == "vus":
        flags.append("clinvar_vus")
    elif tier == "conflicting":
        flags.append("clinvar_conflicting")
    if "splice" in consequence:
        flags.append("splice_region")
    if ((gnomad or {}).get("predictions", {}).get("spliceai_ds_max") or 0) >= 0.5:
        flags.append("spliceai_elevated")

    zyg = sample.get("zygosity") or "unknown zygosity"
    gene_txt = gene or "unknown gene"
    parts = [f"{zyg.replace('_', ' ')} variant"]
    if gene:
        parts.append(f"in {gene}")
    if consequence:
        parts.append(f"({consequence.replace('&', ', ')})")
    if comparison.get("is_novel"):
        parts.append("not observed in gnomAD exomes or ClinVar baseline.")
    elif in_gnomad := comparison.get("in_gnomad"):
        freq = _format_af(af)
        parts.append(f"observed in gnomAD exomes (AF={freq}, AC={_format_num(ac)}/{_format_num(an)}).")
        if comparison.get("in_clinvar"):
            clin_sig = ", ".join(clinvar.get("clinsig") or []) if clinvar else tier
            parts.append(f"ClinVar: {clin_sig}.")
    elif comparison.get("in_clinvar"):
        clin_sig = ", ".join(clinvar.get("clinsig") or []) if clinvar else tier
        parts.append(f"ClinVar only ({clin_sig}).")

    return {
        "summary": " ".join(parts),
        "flags": flags,
        "flag_details": build_flag_details(
            flags,
            gnomad=gnomad,
            clinvar=clinvar,
            sample=sample,
            comparison=comparison,
        ),
        "acmg_automation": "not_evaluated",
    }


def _explain_flag(
    flag: str,
    *,
    gnomad_exomes: dict[str, Any] | None,
    gnomad_genomes: dict[str, Any] | None,
    gnomad: dict[str, Any] | None,
    clinvar: dict[str, Any] | None,
    sample: dict[str, Any] | None,
    comparison: dict[str, Any] | None,
) -> str:
    gnomad = gnomad or gnomad_exomes or gnomad_genomes
    ann = (gnomad or {}).get("annotation") or {}
    consequence = ann.get("consequence") or ""
    pred = (gnomad or {}).get("predictions") or {}
    tier = (clinvar or {}).get("clinical_tier") or "none"
    clinsig = ", ".join(clinvar.get("clinsig") or []) if clinvar else "n/a"
    ex_af = (gnomad_exomes or {}).get("overall", {}).get("af")
    ge_af = (gnomad_genomes or {}).get("overall", {}).get("af")
    af_used = ex_af if ex_af is not None else ge_af
    spliceai = pred.get("spliceai_ds_max")

    if flag == "low_frequency":
        bits = ["rule: overall AF < 0.1% (0.001); exomes AF preferred when present"]
        if ex_af is not None:
            bits.append(f"this variant exomes AF = {_format_af_display(ex_af)}")
            if ex_af < 0.001:
                bits.append("below threshold -> flagged")
            if ge_af is not None:
                if ge_af >= 0.001:
                    bits.append(
                        f"genomes AF = {_format_af_display(ge_af)} is at or above 0.1% but flag follows exomes AF"
                    )
                else:
                    bits.append(f"genomes AF = {_format_af_display(ge_af)}")
        elif af_used is not None:
            bits.append(f"this variant AF = {_format_af_display(af_used)}")
            if af_used < 0.001:
                bits.append("below threshold -> flagged")
        return "; ".join(bits)

    if flag == "common":
        bits = ["rule: overall AF >= 1% (0.01); exomes AF preferred when present"]
        if af_used is not None:
            bits.append(f"this variant AF = {_format_af_display(af_used)} >= 0.01 -> flagged")
        return "; ".join(bits)

    if flag == "exome_genome_af_discordant":
        if ex_af is None or ge_af is None:
            return (
                "rule: |exomes AF - genomes AF| > 50% of max(exomes AF, genomes AF); "
                "both callsets required"
            )
        diff = abs(ex_af - ge_af)
        threshold = max(ex_af, ge_af) * 0.5
        bits = [
            "rule: |exomes AF - genomes AF| > 50% of max(exomes AF, genomes AF)",
            f"exomes AF = {_format_af_display(ex_af)}",
            f"genomes AF = {_format_af_display(ge_af)}",
            f"|diff| = {_format_af_display(diff)}",
            f"threshold = {_format_af_display(threshold)}",
        ]
        if diff > threshold:
            bits.append("difference exceeds threshold -> flagged")
        return "; ".join(bits)

    if flag == "clinvar_pathogenic":
        return (
            "rule: ClinVar clinical_tier is pathogenic or likely_pathogenic; "
            f"this variant tier = {tier}, significance = {clinsig or 'n/a'} -> flagged"
        )

    if flag == "clinvar_vus":
        return (
            "rule: ClinVar clinical_tier is VUS; "
            f"this variant tier = {tier}, significance = {clinsig or 'n/a'} -> flagged"
        )

    if flag == "clinvar_conflicting":
        return (
            "rule: ClinVar clinical_tier is conflicting; "
            f"this variant tier = {tier}, significance = {clinsig or 'n/a'} -> flagged"
        )

    if flag == "clinvar_benign":
        return (
            "rule: ClinVar clinical_tier is benign or likely_benign; "
            f"this variant tier = {tier}, significance = {clinsig or 'n/a'} -> flagged"
        )

    if flag == "splice_region":
        cons = consequence.replace("&", ", ") or "n/a"
        return (
            "rule: consequence annotation contains 'splice'; "
            f"this variant consequence = {cons} -> flagged"
        )

    if flag == "spliceai_elevated":
        score = _format_display(spliceai)
        return (
            "rule: SpliceAI max score >= 0.5; "
            f"this variant SpliceAI max = {score} >= 0.5 -> flagged"
        )

    if flag == "novel_in_sample":
        status = (comparison or {}).get("match_status") or "novel_in_sample"
        return (
            "rule: variant not observed in gnomAD exomes baseline; "
            f"this sample match_status = {status} -> flagged"
        )

    return "rule not documented for this flag"


def build_flag_details(
    flags: list[str],
    *,
    gnomad_exomes: dict[str, Any] | None = None,
    gnomad_genomes: dict[str, Any] | None = None,
    gnomad: dict[str, Any] | None = None,
    clinvar: dict[str, Any] | None = None,
    sample: dict[str, Any] | None = None,
    comparison: dict[str, Any] | None = None,
) -> list[dict[str, str]]:
    return [
        {
            "flag": flag,
            "explanation": _explain_flag(
                flag,
                gnomad_exomes=gnomad_exomes,
                gnomad_genomes=gnomad_genomes,
                gnomad=gnomad,
                clinvar=clinvar,
                sample=sample,
                comparison=comparison,
            ),
        }
        for flag in flags
    ]


def build_reference_interpretation(
    gnomad_exomes: dict[str, Any] | None,
    gnomad_genomes: dict[str, Any] | None,
    clinvar: dict[str, Any] | None,
) -> dict[str, Any]:
    flags: list[str] = []
    gnomad = gnomad_exomes or gnomad_genomes
    gene = (gnomad or {}).get("annotation", {}).get("gene_symbol")
    consequence = (gnomad or {}).get("annotation", {}).get("consequence") or ""
    tier = (clinvar or {}).get("clinical_tier") or "none"

    ex_af = (gnomad_exomes or {}).get("overall", {}).get("af")
    ge_af = (gnomad_genomes or {}).get("overall", {}).get("af")
    af = ex_af if ex_af is not None else ge_af

    if af is not None and af < 0.001:
        flags.append("low_frequency")
    elif af is not None and af >= 0.01:
        flags.append("common")
    if ex_af is not None and ge_af is not None and abs(ex_af - ge_af) > max(ex_af, ge_af) * 0.5:
        flags.append("exome_genome_af_discordant")
    if tier in {"pathogenic", "likely_pathogenic"}:
        flags.append("clinvar_pathogenic")
    elif tier == "vus":
        flags.append("clinvar_vus")
    elif tier == "conflicting":
        flags.append("clinvar_conflicting")
    elif tier in {"benign", "likely_benign"}:
        flags.append("clinvar_benign")
    if "splice" in consequence:
        flags.append("splice_region")
    if ((gnomad or {}).get("predictions", {}).get("spliceai_ds_max") or 0) >= 0.5:
        flags.append("spliceai_elevated")

    parts: list[str] = []
    if gene:
        parts.append(f"Variant in {gene}")
    if consequence:
        parts.append(f"({consequence.replace('&', ', ')})")
    if gnomad_exomes:
        o = gnomad_exomes["overall"]
        parts.append(
            f"gnomAD exomes v4.1: AF={_format_af(o.get('af'))}, AC={_format_num(o.get('ac'))}/{_format_num(o.get('an'))}, "
            f"filter={','.join(gnomad_exomes.get('filter') or ['-'])}."
        )
    if gnomad_genomes:
        o = gnomad_genomes["overall"]
        parts.append(
            f"gnomAD genomes v4.1: AF={_format_af(o.get('af'))}, AC={_format_num(o.get('ac'))}/{_format_num(o.get('an'))}, "
            f"filter={','.join(gnomad_genomes.get('filter') or ['-'])}."
        )
    if not gnomad_exomes and not gnomad_genomes:
        parts.append("not found in local gnomAD exomes or genomes sites VCF.")
    if clinvar:
        parts.append(f"ClinVar significance: {', '.join(clinvar.get('clinsig') or [])}.")
    elif gnomad:
        parts.append("No matching ClinVar record for this allele.")

    return {
        "summary": " ".join(parts),
        "flags": flags,
        "flag_details": build_flag_details(
            flags,
            gnomad_exomes=gnomad_exomes,
            gnomad_genomes=gnomad_genomes,
            gnomad=gnomad,
            clinvar=clinvar,
        ),
        "acmg_automation": "not_evaluated",
    }


def build_variant_brief_text(report: dict[str, Any]) -> str:
    """One-paragraph variant overview for the PDF header (position, gene, databases)."""
    variant = report["variant"]
    gnomad = report.get("gnomad") or {}
    gnomad_exomes = report.get("gnomad_exomes")
    gnomad_genomes = report.get("gnomad_genomes")
    clinvar = report.get("clinvar")
    comparison = report.get("comparison") or {}
    sample = report.get("sample") or {}
    pipeline = report.get("pipeline") or {}
    is_reference = report.get("report_type") == "reference_variant"
    ann = gnomad.get("annotation") or {}

    locus = variant.get("locus") or f"{variant.get('chrom')}:{variant.get('pos'):,}"
    change = variant.get("change") or f"{variant.get('ref')}>{variant.get('alt')}"
    genome = pipeline.get("reference_genome", "GRCh38")
    gene = ann.get("gene_symbol") or "unknown gene"
    consequence = (ann.get("consequence") or "unspecified consequence").replace("&", ", ")
    impact = ann.get("impact") or "unknown"

    variant_id = report.get("variant_id") or variant.get("variant_id")
    if not variant_id:
        chrom_key = str(variant.get("chrom", "")).replace("chr", "")
        variant_id = f"{chrom_key}-{variant.get('pos')}-{variant.get('ref')}-{variant.get('alt')}"

    position_line = (
        f"The position of variant {variant_id} is located in {locus} on {genome}, "
        f"the ref-alt is {change}, in {gene} "
        f"({consequence}; {impact} impact)."
    )

    def _brief_callset_freq(label: str, overall: dict[str, Any] | None) -> str:
        o = overall or {}
        return (
            f"{label} AF is {_format_af(o.get('af'))} "
            f"(AC {_format_num(o.get('ac'))}/{_format_num(o.get('an'))})"
        )

    if comparison.get("is_novel"):
        gnomad_sentence = "It is not observed in the gnomAD v4.1 baseline."
    else:
        freq_bits: list[str] = []
        if gnomad_exomes:
            freq_bits.append(_brief_callset_freq("exomes", gnomad_exomes.get("overall")))
        elif gnomad and str(gnomad.get("callset", "")).startswith("exomes"):
            freq_bits.append(_brief_callset_freq("exomes", gnomad.get("overall")))
        if gnomad_genomes:
            freq_bits.append(_brief_callset_freq("genomes", gnomad_genomes.get("overall")))
        if freq_bits:
            gnomad_sentence = "It is recorded in gnomAD v4.1: " + "; ".join(freq_bits) + "."
        elif comparison.get("in_gnomad") or gnomad:
            gnomad_sentence = "It is present in gnomAD v4.1."
        else:
            gnomad_sentence = "It is not recorded in gnomAD v4.1."

    if clinvar:
        sig = ", ".join(clinvar.get("clinsig") or []) or "reported"
        review = ", ".join(clinvar.get("review_status") or [])
        clinvar_sentence = f"ClinVar lists this allele as {sig}"
        if review:
            clinvar_sentence += f" ({review})"
        clinvar_sentence += "."
    elif comparison.get("in_clinvar") is False or (is_reference and not clinvar):
        clinvar_sentence = "No ClinVar entry matches this exact allele."
    else:
        clinvar_sentence = "ClinVar status not evaluated."

    database_line = f"{gnomad_sentence} {clinvar_sentence}"

    if not is_reference and sample.get("genotype"):
        zyg = (sample.get("zygosity") or "unknown").replace("_", " ")
        sample_line = (
            f"In sample {sample.get('sample_id', 'unknown')} the call is "
            f"{sample.get('genotype')} ({zyg}, DP {_format_display(sample.get('depth'))})."
        )
        return f"{position_line} {sample_line} {database_line}"

    return f"{position_line} {database_line}"


def _fetch_reference_layers(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> dict[str, Any]:
    chrom = normalize_chrom(chrom)
    gnomad_exomes = fetch_gnomad_variant(chrom, pos, ref, alt, callset="exomes")
    gnomad_genomes = fetch_gnomad_variant(chrom, pos, ref, alt, callset="genomes")
    clinvar = fetch_clinvar_variant(chrom, pos, ref, alt)
    primary = gnomad_exomes or gnomad_genomes
    rsid = (clinvar or {}).get("rsid")
    if rsid and not str(rsid).startswith("rs"):
        rsid = f"rs{rsid}"
    gene_symbol = (primary or {}).get("annotation", {}).get("gene_symbol")
    canonical_vid = f"{chrom.replace('chr', '')}-{pos}-{ref}-{alt}"
    return {
        "chrom": chrom,
        "canonical_vid": canonical_vid,
        "gnomad_exomes": gnomad_exomes,
        "gnomad_genomes": gnomad_genomes,
        "clinvar": clinvar,
        "gnomad": primary,
        "gene_symbol": gene_symbol,
        "rsid": rsid,
    }


def _data_sources_block(
    *,
    gnomad_exomes: dict[str, Any] | None,
    gnomad_genomes: dict[str, Any] | None,
    clinvar: dict[str, Any] | None,
    sample_vcf: Path | None = None,
) -> dict[str, Any]:
    return {
        "gnomad_exomes_vcf": gnomad_exomes is not None,
        "gnomad_genomes_vcf": gnomad_genomes is not None,
        "gnomad_exomes_vcf_path": str(GNOMAD_EXOMES_VCF_DIR),
        "gnomad_genomes_vcf_path": str(GNOMAD_GENOMES_VCF_DIR),
        "gnomad_exomes_tsv_path": str(GNOMAD_EXOMES_TSV_DIR),
        "gnomad_genomes_tsv_path": str(GNOMAD_GENOMES_TSV_DIR),
        "clinvar_vcf_path": str(CLINVAR_VCF),
        "constraint_tsv_path": str(CONSTRAINT_TSV),
        "clinvar_vcf": clinvar is not None,
        "constraint_tsv": CONSTRAINT_TSV.exists(),
        "sample_vcf": sample_vcf is not None,
        "sample_vcf_path": str(sample_vcf) if sample_vcf else None,
        "notes": (
            "Primary annotation and exome frequencies from exomes sites VCF; "
            "genome frequencies from genomes sites VCF. "
            "Supplementary TSV bundles under gnomAD/4.1/tsv/tsv/{exomes,genomes}."
        ),
    }


def _apply_compare_row_to_sample(sample: dict[str, Any], compare_row: dict[str, Any] | None) -> None:
    if not compare_row:
        return
    if compare_row.get("sample_gt"):
        sample["genotype"] = compare_row.get("sample_gt")
    if compare_row.get("sample_zygosity"):
        z = str(compare_row["sample_zygosity"])
        sample["zygosity"] = "heterozygous" if z == "het" else (
            "homozygous_alt" if z == "hom_alt" else z
        )
    if compare_row.get("sample_dp") is not None:
        sample["depth"] = int(float(compare_row["sample_dp"]))
    if compare_row.get("variant_id"):
        sample["compare_variant_id"] = compare_row.get("variant_id")


def build_reference_variant_report(variant_id: str) -> dict[str, Any]:
    chrom, pos, ref, alt = parse_gnomad_variant_id(variant_id)
    layers = _fetch_reference_layers(chrom, pos, ref, alt)
    if (
        layers["gnomad_exomes"] is None
        and layers["gnomad_genomes"] is None
        and layers["clinvar"] is None
    ):
        raise ValueError(f"variant not found in gnomAD or ClinVar: {variant_id}")

    from _gene_info import resolve_gene_info

    gene_info = resolve_gene_info(layers["gene_symbol"])
    canonical_vid = layers["canonical_vid"]
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "report_type": "reference_variant",
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "variant_id": canonical_vid,
        "gnomad_variant_page": gnomad_variant_page_url(canonical_vid),
        "pipeline": {
            "name": "exome_baseline",
            "mode": "reference_only",
            "baseline_callset": f"gnomad_{DATA_RELEASE['gnomad_callset']}_{DATA_RELEASE['gnomad']}",
            "reference_genome": DATA_RELEASE["reference_genome"],
        },
        "data_sources": _data_sources_block(
            gnomad_exomes=layers["gnomad_exomes"],
            gnomad_genomes=layers["gnomad_genomes"],
            clinvar=layers["clinvar"],
        ),
        "variant": {
            "variant_id": canonical_vid,
            "chrom": layers["chrom"],
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "change": f"{ref}>{alt}",
            "locus": f"{layers['chrom']}:{pos:,}",
            "rsid": layers["rsid"],
            "variant_type": (layers["clinvar"] or {}).get("variant_type"),
        },
        "gnomad": layers["gnomad"],
        "gnomad_exomes": layers["gnomad_exomes"],
        "gnomad_genomes": layers["gnomad_genomes"],
        "clinvar": layers["clinvar"],
        "gene_info": gene_info,
        "interpretation": build_reference_interpretation(
            layers["gnomad_exomes"],
            layers["gnomad_genomes"],
            layers["clinvar"],
        ),
    }


def build_variant_report(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    sample_id: str,
    sample_vcf: Path,
    compare_row: dict[str, Any] | None = None,
    *,
    compare_dir: Path | None = None,
) -> dict[str, Any]:
    layers = _fetch_reference_layers(chrom, pos, ref, alt)
    chrom = layers["chrom"]
    sample = lookup_sample_variant(sample_vcf, sample_id, chrom, pos, ref, alt)
    _apply_compare_row_to_sample(sample, compare_row)

    in_gnomad_exomes = layers["gnomad_exomes"] is not None
    in_gnomad_genomes = layers["gnomad_genomes"] is not None
    in_clinvar = layers["clinvar"] is not None
    tier = (layers["clinvar"] or {}).get("clinical_tier") or (
        compare_row.get("clinical_tier") if compare_row else "none"
    )
    row_for_score = {
        "is_novel": not (in_gnomad_exomes or in_clinvar),
        "clinical_tier": tier,
        "gnomad_af": (layers["gnomad_exomes"] or layers["gnomad"] or {}).get("overall", {}).get("af"),
        "sample_zygosity": sample.get("zygosity"),
    }
    if compare_row and compare_row.get("priority_score") is not None:
        score = float(compare_row["priority_score"])
    else:
        score = priority_score(row_for_score)

    comparison = build_comparison_block(
        in_gnomad_exomes,
        in_clinvar,
        tier,
        score,
        in_gnomad_genomes=in_gnomad_genomes,
    )
    canonical_vid = layers["canonical_vid"]
    from _gene_info import resolve_gene_info

    gene_info = resolve_gene_info(layers["gene_symbol"])
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "report_type": "sample_variant",
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "variant_id": canonical_vid,
        "gnomad_variant_page": gnomad_variant_page_url(canonical_vid),
        "pipeline": {
            "name": "exome_baseline",
            "mode": "sample_vs_baseline",
            "sample_id": sample_id,
            "sample_vcf": str(sample_vcf),
            "compare_dir": str(compare_dir) if compare_dir else None,
            "baseline_callset": f"gnomad_{DATA_RELEASE['gnomad_callset']}_{DATA_RELEASE['gnomad']}",
            "reference_genome": DATA_RELEASE["reference_genome"],
        },
        "data_sources": _data_sources_block(
            gnomad_exomes=layers["gnomad_exomes"],
            gnomad_genomes=layers["gnomad_genomes"],
            clinvar=layers["clinvar"],
            sample_vcf=sample_vcf,
        ),
        "variant": {
            "variant_id": canonical_vid,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "change": f"{ref}>{alt}",
            "locus": f"{chrom}:{pos:,}",
            "rsid": layers["rsid"],
            "variant_type": (layers["clinvar"] or {}).get("variant_type"),
        },
        "sample": sample,
        "comparison": comparison,
        "gene_info": gene_info,
        "gnomad": layers["gnomad"],
        "gnomad_exomes": layers["gnomad_exomes"],
        "gnomad_genomes": layers["gnomad_genomes"],
        "clinvar": layers["clinvar"],
        "interpretation": build_interpretation(
            sample,
            comparison,
            layers["gnomad"],
            layers["clinvar"],
        ),
    }


build_sample_variant_report = build_variant_report


def write_report_json(report: dict[str, Any], path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    path.write_text(text, encoding="utf-8")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def load_report_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _pdf_safe(text: Any) -> str:
    if text is None:
        return "-"
    s = str(text)
    s = s.replace("\u2013", "-").replace("\u2014", "-").replace("\u00d7", "x")
    return s.encode("latin-1", errors="replace").decode("latin-1")


def _pdf_safe_cjk(text: Any) -> str:
    if text is None:
        return "-"
    return str(text)


_PDF_RENDER_CTX: dict[str, Any] = {"locale": "en", "titles": None}


def _pdf_titles() -> Any:
    from _variant_report_i18n import get_pdf_titles

    titles = _PDF_RENDER_CTX.get("titles")
    if titles is not None:
        return titles
    return get_pdf_titles(_PDF_RENDER_CTX.get("locale", "en"))


def _pdf_locale() -> str:
    return _PDF_RENDER_CTX.get("locale", "en")


def _pdf_register_cjk(pdf: Any) -> None:
    from _variant_report_i18n import NOTO_CJK_BOLD, NOTO_CJK_REGULAR, NOTO_CJK_SC_INDEX

    if getattr(pdf, "_cjk_registered", False):
        return
    pdf.add_font("NotoCJK", "", NOTO_CJK_REGULAR, collection_font_number=NOTO_CJK_SC_INDEX)
    pdf.add_font("NotoCJK", "B", NOTO_CJK_BOLD, collection_font_number=NOTO_CJK_SC_INDEX)
    pdf._cjk_registered = True


def _pdf_bind_locale(pdf: Any, locale: str) -> None:
    from _variant_report_i18n import get_pdf_titles

    pdf._locale = locale
    pdf._titles = get_pdf_titles(locale)
    _PDF_RENDER_CTX["locale"] = locale
    _PDF_RENDER_CTX["titles"] = pdf._titles
    if locale == "zh":
        _pdf_register_cjk(pdf)


_PDF_HRULE = "B"  # horizontal rules only (no vertical cell borders)
_PDF_TABLE_HEAD_FILL = (64, 64, 64)
_PDF_CONTENT_SECTION_FILL = (255, 236, 210)  # light orange
_PDF_SUBSECTION_FILL = (242, 242, 242)  # light gray
_PDF_SEPARATOR_COLOR = (230, 81, 0)  # orange rules (distinct from table black borders)
_PDF_KV_LABEL_W = 28  # default label column width (mm)
_PDF_HALF_LINE = 2.5
_PDF_PAGE_FOOTER_H = 12  # reserved bottom area for "N of M" footer
_PDF_TABLE_FONT_SIZE = 8
_PDF_TABLE_LINE_H = 5
# Column widths (mm) — edit these to tune PDF layout
_PDF_CALLSET_COL_W = (28, 22, 28, 28)  # Callset, AC, AN, AF; Filter uses remaining width
_PDF_FREQ_COL_W = (12, 14, 18, 24, 18)  # Pop, AC, AN, AF, rel. AF (per side-by-side block)
_PDF_FREQ_SIDE_GAP = 4  # horizontal gap between exomes / genomes tables
_PDF_IDEOGRAM_WIDTH_FRAC = 0.80
_PDF_IDEOGRAM_LABEL_GAP = 2.5
_PDF_FREQ_CHART_PAD_R = 3  # right inset so max bar does not touch column edge (mm)
_PDF_DATA_SOURCES_FONT_SIZE = 6
_PDF_DATA_SOURCES_LINE_H = 2.1
_PDF_GLOSSARY_FONT_SIZE = 6
_PDF_GLOSSARY_LINE_H = 1.92
_PDF_INTERPRETATION_FLAG_LINE_H = 5
_PDF_INTERPRETATION_RULE_LINE_H = 2.5

# Term definitions for ANNOTATION & SCORES (rendered at end of report PDF).
_PDF_ANNOTATION_GLOSSARY: list[tuple[str, str]] = [
    ("Gene", "HGNC gene symbol at the variant locus; full name is shown in parentheses when available."),
    (
        "Function",
        "Brief gene function summary from curated gene databases (e.g. NCBI Gene via MyGene.info).",
    ),
    ("Transcript", "Canonical Ensembl transcript (MANE Select when available) used for primary HGVS notation."),
    ("Consequence", "Sequence Ontology term(s) describing the effect on the transcript (e.g. intron_variant)."),
    ("Impact", "VEP predicted consequence severity: HIGH, MODERATE, LOW, or MODIFIER."),
    ("Exon", "Affected exon as current/total exon number; n/a when the variant is not in an exon."),
    ("Intron", "Affected intron as current/total intron number for intronic variants."),
    ("HGVS c.", "HGVS coding-DNA description relative to the listed transcript."),
    ("HGVS p.", "HGVS protein change; n/a when no protein-level change is annotated."),
    (
        "Constr. tx",
        "Transcript used for gnomAD v4 gene constraint metrics when it differs from the canonical transcript.",
    ),
    (
        "pLI",
        "Probability of loss-of-function intolerance (gnomAD); values near 1 indicate the gene tolerates LoF poorly.",
    ),
    (
        "LOEUF",
        "Loss-of-function observed/expected upper-bound fraction (gnomAD); lower values indicate stronger LoF constraint.",
    ),
    (
        "mis_z",
        "Missense constraint Z-score (gnomAD); higher values indicate depletion of missense variation.",
    ),
    (
        "CADD Phred",
        "PHRED-scaled CADD score integrating multiple annotations; higher values suggest greater deleteriousness.",
    ),
    (
        "REVEL max",
        "Maximum REVEL missense pathogenicity score (0-1); higher values support pathogenic missense effects.",
    ),
    (
        "SpliceAI max",
        "Maximum SpliceAI delta score (0–1) for altered splicing; values >=0.5 are often considered notable.",
    ),
    (
        "SIFT max",
        "Maximum SIFT score across affected transcripts; lower values generally indicate tolerated substitutions.",
    ),
    (
        "PolyPhen max",
        "Maximum PolyPhen-2 score predicting damaging missense substitutions.",
    ),
    (
        "phyloP",
        "Maximum phyloP nucleotide conservation score; positive values indicate evolutionary conservation.",
    ),
]

# Sample finding and baseline comparison terms (appended for sample_variant reports).
_PDF_SAMPLE_GLOSSARY: list[tuple[str, str]] = [
    (
        "Sample call",
        "Genotype and zygosity observed in the individual sample VCF at this locus.",
    ),
    (
        "Genotype",
        "VCF GT field for the sample at this locus (e.g. 0/1 = one reference and one alternate allele).",
    ),
    (
        "Zygosity",
        "Call class derived from genotype: heterozygous, homozygous_alt, homozygous_ref, etc.",
    ),
    (
        "Depth",
        "Sequencing depth (FORMAT DP) at the variant site in the sample VCF.",
    ),
    (
        "Allele depth",
        "Read counts supporting reference and alternate alleles (FORMAT AD: ref / alt).",
    ),
    (
        "Quality",
        "Variant quality from the sample VCF record (VCF QUAL or sample GQ when available).",
    ),
    (
        "FILTER",
        "VCF FILTER status for the sample record (e.g. PASS or failing filter labels).",
    ),
    (
        "Match status",
        "Machine-readable classification vs gnomAD exomes + ClinVar baseline "
        "(e.g. known_gnomad_clinvar, novel_in_sample).",
    ),
    (
        "Verdict",
        "Short human-readable summary of how the sample allele relates to the baseline.",
    ),
    (
        "Novel",
        "Yes when the allele is absent from both gnomAD exomes and ClinVar baseline (novel_in_sample).",
    ),
    (
        "Priority",
        "Heuristic triage score and tier (routine / elevated / urgent); higher for novel or clinical variants.",
    ),
    (
        "In gnomAD exomes",
        "Whether this exact allele is present in the local gnomAD v4.1 exomes sites VCF (WES baseline).",
    ),
    (
        "In gnomAD genomes",
        "Whether this exact allele is present in the local gnomAD v4.1 genomes sites VCF (supplementary; "
        "not used for novel determination).",
    ),
    (
        "In ClinVar",
        "Whether this exact allele is present in the local ClinVar GRCh38 VCF.",
    ),
]


class _VariantReportPDF:
    """FPDF subclass with English page numbers in footer."""

    @staticmethod
    def create() -> Any:
        from fpdf import FPDF

        class VariantReportPDF(FPDF):
            def footer(self) -> None:
                self.set_y(-12)
                self.set_font("Helvetica", "", 8)
                self.set_text_color(120, 120, 120)
                self.cell(0, 8, f"{self.page_no()} of {{nb}}", align="C")

        return VariantReportPDF()


def _pdf_section_title(pdf: Any, title: str, content_w: float) -> None:
    pdf.ln(1)
    use_cjk = _pdf_locale() == "zh"
    if use_cjk:
        pdf.set_font("NotoCJK", "B", 10)
        safe_title = _pdf_safe_cjk(title)
    else:
        pdf.set_font("Helvetica", "B", 10)
        safe_title = _pdf_safe(title)
    pdf.set_fill_color(*_PDF_CONTENT_SECTION_FILL)
    pdf.set_text_color(0, 0, 0)
    pdf.cell(content_w, 6, safe_title, new_x="LMARGIN", new_y="NEXT", fill=True)


def _pdf_orange_section(pdf: Any, content_w: float, title: str) -> None:
    """Section heading above an orange horizontal rule."""
    pdf.ln(_PDF_TABLE_LINE_H)
    _pdf_section_title(pdf, title, content_w)
    _pdf_dashed_rule(pdf, content_w)


def _pdf_subsection_title(pdf: Any, title: str, content_w: float) -> None:
    pdf.ln(_PDF_HALF_LINE)
    use_cjk = _pdf_locale() == "zh"
    if use_cjk:
        pdf.set_font("NotoCJK", "B", 10)
        safe_title = _pdf_safe_cjk(title)
    else:
        pdf.set_font("Helvetica", "B", 10)
        safe_title = _pdf_safe(title)
    pdf.set_fill_color(*_PDF_SUBSECTION_FILL)
    pdf.set_text_color(0, 0, 0)
    pdf.cell(content_w, 6, safe_title, new_x="LMARGIN", new_y="NEXT", fill=True)


def _pdf_table_section_title(
    pdf: Any,
    number: int,
    title: str,
    content_w: float,
    *,
    start_x: float | None = None,
) -> None:
    pdf.ln(1)
    pdf.set_font("Helvetica", "B", 10)
    pdf.set_text_color(0, 0, 0)
    label = f"{number}. {title}"
    pdf.set_x(start_x if start_x is not None else pdf.l_margin)
    pdf.cell(content_w, 6, _pdf_safe(label), new_x="LMARGIN", new_y="NEXT", border="B")


def _pdf_table_header_row_at(
    pdf: Any,
    start_x: float,
    widths: list[float],
    values: list[Any],
    *,
    font_size: int = _PDF_TABLE_FONT_SIZE,
    border: str = _PDF_HRULE,
) -> None:
    pdf.set_fill_color(*_PDF_TABLE_HEAD_FILL)
    pdf.set_text_color(255, 255, 255)
    pdf.set_font("Helvetica", "B", font_size)
    pdf.set_x(start_x)
    for width, value in zip(widths, values):
        pdf.cell(width, 5, _pdf_safe(value), border=border, fill=True)
    pdf.ln()
    pdf.set_text_color(0, 0, 0)


def _pdf_table_header_row(
    pdf: Any,
    widths: list[float],
    values: list[Any],
    *,
    font_size: int = _PDF_TABLE_FONT_SIZE,
    border: str = _PDF_HRULE,
) -> None:
    pdf.set_fill_color(*_PDF_TABLE_HEAD_FILL)
    pdf.set_text_color(255, 255, 255)
    pdf.set_font("Helvetica", "B", font_size)
    for width, value in zip(widths, values):
        pdf.cell(width, 5, _pdf_safe(value), border=border, fill=True)
    pdf.ln()
    pdf.set_text_color(0, 0, 0)


def _pdf_header_id_row(
    pdf: Any,
    left_text: str,
    right_text: str,
    content_w: float,
    *,
    left_size: int = 18,
    right_size: int = 9,
    row_h: float = 8,
) -> None:
    pdf.set_text_color(0, 0, 0)
    pdf.set_font("Helvetica", "B", left_size)
    left = _pdf_safe(left_text)
    left_w = pdf.get_string_width(left)
    pdf.set_x(pdf.l_margin)
    pdf.cell(left_w + 1, row_h, left)
    pdf.set_font("Helvetica", "", right_size)
    right = _pdf_safe(right_text)
    right_w = pdf.get_string_width(right)
    gap = max(2, content_w - left_w - 1 - right_w)
    pdf.cell(gap, row_h, "")
    pdf.cell(right_w, row_h, right, new_x="LMARGIN", new_y="NEXT")


def _pdf_text_right(
    pdf: Any,
    text: str,
    content_w: float,
    *,
    height: float = 4,
    font_size: float = 6,
    color: tuple[int, int, int] = (90, 90, 90),
) -> None:
    pdf.set_text_color(*color)
    safe = _pdf_safe(text)
    size = font_size
    pdf.set_font("Helvetica", "", size)
    while size > 4 and pdf.get_string_width(safe) > content_w:
        size -= 0.5
        pdf.set_font("Helvetica", "", size)
    text_w = pdf.get_string_width(safe)
    pdf.set_x(pdf.l_margin + max(0, content_w - text_w))
    pdf.cell(min(text_w, content_w), height, safe, new_x="LMARGIN", new_y="NEXT")
    pdf.set_text_color(0, 0, 0)


def _pdf_dashed_rule(pdf: Any, content_w: float, *, half_line: float = 2.5) -> None:
    pdf.ln(half_line)
    y = pdf.get_y()
    pdf.set_draw_color(*_PDF_SEPARATOR_COLOR)
    pdf.set_line_width(0.2)
    pdf.line(pdf.l_margin, y, pdf.l_margin + content_w, y)
    pdf.set_draw_color(0, 0, 0)
    pdf.set_y(y + half_line)


def _pdf_page_break_limit(pdf: Any) -> float:
    return pdf.h - pdf.b_margin


def _pdf_ensure_vertical_space(pdf: Any, height: float) -> None:
    if pdf.get_y() + height > _pdf_page_break_limit(pdf):
        pdf.add_page()


def _pdf_callset_summary_height_estimate() -> float:
    title_h = 1 + 6
    header_h = 5
    rows_h = 2 * _PDF_TABLE_LINE_H
    tail_h = _PDF_HALF_LINE
    return title_h + header_h + rows_h + tail_h


def _pdf_orange_section_height_estimate() -> float:
    return _PDF_TABLE_LINE_H + _PDF_HALF_LINE * 2 + 1 + 6


def _pdf_gnomad_callset_freq_height_estimate(
    exomes: dict[str, Any] | None,
    genomes: dict[str, Any] | None,
) -> float:
    freq_h = max(
        _pdf_freq_block_height_estimate(exomes),
        _pdf_freq_block_height_estimate(genomes),
    ) + _PDF_HALF_LINE
    return _pdf_orange_section_height_estimate() + _pdf_callset_summary_height_estimate() + freq_h


def _pdf_freq_block_height_estimate(gnomad_block: dict[str, Any] | None) -> float:
    title_h = 1 + 6
    if not gnomad_block:
        return title_h + 5
    n_rows = 1 + len(gnomad_block.get("populations") or {})
    body_h = 5 + n_rows * _PDF_TABLE_LINE_H
    note_h = 2 * (_PDF_TABLE_LINE_H - 1) + 1
    if gnomad_block.get("source_vcf"):
        note_h += _PDF_TABLE_LINE_H - 1.5
    return title_h + body_h + note_h


def _pdf_register_dejavu(pdf: Any) -> None:
    if getattr(pdf, "_dejavu_registered", False):
        return
    regular = Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
    bold = Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf")
    if not regular.is_file():
        raise RuntimeError(f"DejaVu font not found: {regular}")
    pdf.add_font("DejaVu", "", str(regular))
    if bold.is_file():
        pdf.add_font("DejaVu", "B", str(bold))
    pdf._dejavu_registered = True


def _pdf_draw_presence_checkbox(pdf: Any, x: float, y: float, size: float, present: bool) -> None:
    if present:
        pdf.set_draw_color(34, 160, 80)
        pdf.set_text_color(34, 160, 80)
        symbol = "\u2713"
        font_style = "B"
    else:
        pdf.set_draw_color(180, 180, 180)
        pdf.set_text_color(140, 140, 140)
        symbol = "-"
        font_style = ""
    pdf.set_line_width(0.35)
    pdf.rect(x, y, size, size, style="D")
    pdf.set_font("DejaVu", font_style, size * 2.4)
    sym_w = pdf.get_string_width(symbol)
    pdf.set_xy(x + (size - sym_w) / 2, y + size * 0.08)
    pdf.cell(sym_w, size * 0.85, symbol, align="C")
    pdf.set_draw_color(0, 0, 0)
    pdf.set_line_width(0.2)


def _pdf_presence_status_line(
    pdf: Any,
    content_w: float,
    *,
    in_exomes: bool,
    in_genomes: bool,
    in_clinvar: bool,
) -> None:
    _pdf_register_dejavu(pdf)
    pdf.ln(_PDF_HALF_LINE)
    line_h = 6
    box_size = 3.2
    item_gap = 5
    y0 = pdf.get_y()
    x = pdf.l_margin
    pdf.set_font("Helvetica", "B", 10)
    pdf.set_text_color(12, 94, 196)
    prefix = "Present in:"
    pdf.set_xy(x, y0)
    pdf.cell(pdf.get_string_width(prefix) + 2, line_h, prefix)
    x = pdf.get_x() + 1
    for label, present in (
        ("gnomAD exomes", in_exomes),
        ("gnomAD genomes", in_genomes),
        ("ClinVar", in_clinvar),
    ):
        x += item_gap
        box_y = y0 + (line_h - box_size) / 2
        _pdf_draw_presence_checkbox(pdf, x, box_y, box_size, present)
        x += box_size + 1.5
        pdf.set_xy(x, y0)
        pdf.set_font("Helvetica", "B", 10)
        pdf.set_text_color(12, 94, 196)
        pdf.cell(pdf.get_string_width(label) + 0.5, line_h, label)
        x = pdf.get_x()
    pdf.set_y(y0 + line_h)
    pdf.set_text_color(0, 0, 0)


def _is_sample_report(report: dict[str, Any]) -> bool:
    return report.get("report_type") in ("sample_variant", "single_variant")


def _pdf_sample_call_line(
    pdf: Any,
    sample: dict[str, Any],
    comparison: dict[str, Any],
    content_w: float,
) -> None:
    pdf.ln(_PDF_HALF_LINE)
    gt = sample.get("genotype")
    zyg = (sample.get("zygosity") or "").replace("_", " ")
    if gt:
        call_text = f"{gt} ({zyg})" if zyg else str(gt)
    else:
        call_text = "not called"
    pdf.set_font("Helvetica", "B", 10)
    if comparison.get("is_novel"):
        pdf.set_text_color(230, 81, 0)
    else:
        pdf.set_text_color(12, 94, 196)
    pdf.cell(content_w, 6, _pdf_safe(f"Sample call: {call_text}"), new_x="LMARGIN", new_y="NEXT")
    pdf.set_text_color(0, 0, 0)


def _pdf_sample_baseline_comparison_section(
    pdf: Any,
    sample: dict[str, Any],
    comparison: dict[str, Any],
    content_w: float,
) -> None:
    titles = _pdf_titles()
    _pdf_orange_section(pdf, content_w, titles.section_sample_baseline_comparison)
    _pdf_subsection_title(pdf, titles.subsection_sample_finding, content_w)
    ad = sample.get("allele_depth") or {}
    ad_value = None
    if ad:
        ad_value = f"ref {_format_display(ad.get('ref'))} / alt {_format_display(ad.get('alt'))}"
    finding_fields: list[tuple[str, Any]] = [
        ("Genotype", sample.get("genotype")),
        ("Zygosity", sample.get("zygosity")),
        ("Depth", _format_display(sample.get("depth"))),
    ]
    if ad_value:
        finding_fields.append(("Allele depth", ad_value))
    finding_fields.extend(
        [
            ("Quality", _format_display(sample.get("quality"))),
            ("FILTER", sample.get("filter")),
        ]
    )
    _pdf_kv_block(pdf, finding_fields, content_w)

    in_exomes = comparison.get("in_gnomad_exomes", comparison.get("in_gnomad"))
    in_genomes = comparison.get("in_gnomad_genomes", False)
    _pdf_subsection_title(pdf, titles.subsection_baseline_verdict, content_w)
    _pdf_kv_block(
        pdf,
        [
            ("Match status", comparison.get("match_status")),
            ("Verdict", comparison.get("verdict")),
            ("Novel", "Yes" if comparison.get("is_novel") else "No"),
            (
                "Priority",
                f"{_format_num(comparison.get('priority_score'))} ({comparison.get('priority_tier')})",
            ),
            ("In gnomAD exomes", "Yes" if in_exomes else "No"),
            ("In gnomAD genomes", "Yes" if in_genomes else "No"),
            ("In ClinVar", "Yes" if comparison.get("in_clinvar") else "No"),
        ],
        content_w,
    )


def _has_gnomad_frequency_data(
    gnomad_exomes: dict[str, Any] | None,
    gnomad_genomes: dict[str, Any] | None,
) -> bool:
    return gnomad_exomes is not None or gnomad_genomes is not None


def _has_gnomad_annotation(gnomad: dict[str, Any]) -> bool:
    if not gnomad:
        return False
    ann = gnomad.get("annotation") or {}
    return bool(ann.get("gene_symbol") or ann.get("consequence") or ann.get("canonical_transcript"))


def _pdf_variant_brief(pdf: Any, text: str, content_w: float) -> None:
    pdf.ln(_PDF_HALF_LINE)
    pdf.set_x(pdf.l_margin)
    pdf.set_font("Helvetica", "", 9)
    pdf.set_text_color(50, 50, 50)
    pdf.multi_cell(content_w, 5, _pdf_safe(text), align="L")
    pdf.set_text_color(0, 0, 0)
    pdf.ln(_PDF_HALF_LINE)


def _pdf_kv_line(pdf: Any, label: str, value: Any, content_w: float, label_w: float = 42) -> None:
    del label_w  # stacked layout uses full width
    pdf.set_x(pdf.l_margin)
    pdf.set_font("Helvetica", "B", 8)
    pdf.cell(content_w, 5, _pdf_safe(label), new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 8)
    pdf.multi_cell(content_w, 5, _pdf_safe(value))


def _pdf_table_tail_spacer(pdf: Any) -> None:
    pdf.ln(_PDF_HALF_LINE)


def _pdf_kv_inline(
    pdf: Any,
    label: str,
    value: Any,
    content_w: float,
    label_w: float = _PDF_KV_LABEL_W,
    *,
    value_color: tuple[int, int, int] | None = None,
) -> None:
    pdf.set_x(pdf.l_margin)
    if value is None or value == "" or value == "-":
        value = "n/a"
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        value = _format_num(value)
    value_text = _pdf_safe(value)
    value_w = max(20, content_w - label_w)
    row_h = 5

    pdf.set_font("Helvetica", "B", 8)
    pdf.cell(label_w, row_h, _pdf_safe(f"{label}:"))
    pdf.set_font("Helvetica", "", 8)
    if value_color:
        pdf.set_text_color(*value_color)
    if pdf.get_string_width(value_text) <= value_w - 0.5:
        pdf.cell(value_w, row_h, value_text, new_x="LMARGIN", new_y="NEXT")
    else:
        pdf.multi_cell(value_w, row_h, value_text)
        pdf.set_x(pdf.l_margin)
    if value_color:
        pdf.set_text_color(0, 0, 0)

    y = pdf.get_y()
    pdf.set_draw_color(210, 210, 210)
    pdf.set_line_width(0.15)
    pdf.line(pdf.l_margin, y, pdf.l_margin + content_w, y)
    pdf.ln(0.4)


def _pdf_place_doc_footer(pdf: Any, text: str, content_w: float) -> None:
    """Place generated metadata; push down only when close to page bottom."""
    line_h = 4
    spacer = _PDF_TABLE_LINE_H
    target_y = pdf.h - _PDF_PAGE_FOOTER_H - spacer - line_h - 3
    gap = target_y - pdf.get_y()
    if 0 < gap <= 35:
        pdf.ln(gap)
    pdf.ln(spacer)
    pdf.set_x(pdf.l_margin)
    pdf.set_font("Helvetica", "", 6)
    pdf.set_text_color(120, 120, 120)
    safe = _pdf_safe(text)
    if pdf.get_string_width(safe) > content_w - 0.5:
        pdf.multi_cell(content_w, line_h, safe, align="L")
    else:
        pdf.cell(content_w, line_h, safe, new_x="LMARGIN", new_y="NEXT")
    pdf.set_text_color(0, 0, 0)


def _pdf_wrap_lines(
    pdf: Any,
    text: Any,
    width: float,
    *,
    font_size: int = 7,
    font_style: str = "",
) -> list[str]:
    pdf.set_font("Helvetica", font_style, font_size)
    safe = _pdf_safe(text)
    if not safe or safe == "-":
        return ["-"]
    if pdf.get_string_width(safe) <= width - 0.5:
        return [safe]
    words = safe.replace(",", " ,").split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = word if not current else f"{current} {word}"
        if pdf.get_string_width(candidate) <= width - 0.5:
            current = candidate
            continue
        if current:
            lines.append(current)
        if pdf.get_string_width(word) <= width - 0.5:
            current = word
        else:
            chunk = ""
            for ch in word:
                trial = f"{chunk}{ch}"
                if pdf.get_string_width(trial) <= width - 0.5:
                    chunk = trial
                else:
                    if chunk:
                        lines.append(chunk)
                    chunk = ch
            current = chunk
    if current:
        lines.append(current)
    return lines or ["-"]


def _pdf_table_row_wrapped_at(
    pdf: Any,
    start_x: float,
    widths: list[float],
    values: list[Any],
    *,
    bold: bool = False,
    font_size: int = 7,
    line_h: float = 4,
    border: str = _PDF_HRULE,
) -> float:
    style = "B" if bold else ""
    wrapped = [
        _pdf_wrap_lines(pdf, _format_display(value), width, font_size=font_size, font_style=style)
        for width, value in zip(widths, values)
    ]
    row_lines = max(len(lines) for lines in wrapped)
    for line_idx in range(row_lines):
        pdf.set_font("Helvetica", style, font_size)
        pdf.set_x(start_x)
        for width, lines in zip(widths, wrapped):
            chunk = lines[line_idx] if line_idx < len(lines) else ""
            cell_border = border if line_idx == row_lines - 1 else ""
            pdf.cell(width, line_h, chunk, border=cell_border)
        pdf.ln(line_h)
    return row_lines * line_h


def _pdf_draw_h_bar(
    pdf: Any,
    x: float,
    y: float,
    w: float,
    h: float,
    fraction: float,
    *,
    highlight: bool = False,
    total: bool = False,
) -> None:
    fraction = max(0.0, min(1.0, fraction))
    if fraction <= 0:
        return
    bar_w = max(0.4, w * fraction)
    if total:
        pdf.set_fill_color(100, 100, 100)
    elif highlight:
        pdf.set_fill_color(230, 81, 0)
    else:
        pdf.set_fill_color(12, 94, 196)
    pdf.rect(x, y, bar_w, h, style="F")
    pdf.set_draw_color(0, 0, 0)


def _pdf_table_row_wrapped(
    pdf: Any,
    widths: list[float],
    values: list[Any],
    *,
    bold: bool = False,
    font_size: int = 7,
    line_h: float = 4,
    border: str = _PDF_HRULE,
) -> None:
    style = "B" if bold else ""
    wrapped = [
        _pdf_wrap_lines(
            pdf,
            _format_display(value) if not (isinstance(value, str) and value.startswith("#")) else value,
            width,
            font_size=font_size,
            font_style=style,
        )
        for width, value in zip(widths, values)
    ]
    row_lines = max(len(lines) for lines in wrapped)
    for line_idx in range(row_lines):
        pdf.set_font("Helvetica", style, font_size)
        pdf.set_x(pdf.l_margin)
        for col_idx, (width, lines) in enumerate(zip(widths, wrapped)):
            chunk = lines[line_idx] if line_idx < len(lines) else ""
            cell_border = border if line_idx == row_lines - 1 else ""
            pdf.cell(width, line_h, chunk, border=cell_border)
        pdf.ln(line_h)


def _pdf_format_optional(value: Any) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return _format_num(value)
    return str(value)


def _annotation_score_fields(
    gnomad: dict[str, Any],
    gene_info: dict[str, Any] | None = None,
) -> list[tuple[str, Any]]:
    from _gene_info import format_gene_annotation_value, format_gene_function_value

    ann = gnomad.get("annotation") or {}
    constraint = gnomad.get("gene_constraint") or {}
    pred = gnomad.get("predictions") or {}
    gene_value = format_gene_annotation_value(gene_info) or ann.get("gene_symbol")
    fields: list[tuple[str, Any]] = [
        ("Gene", gene_value),
    ]
    function_value = format_gene_function_value(gene_info)
    if function_value:
        fields.append(("Function", function_value))
    fields.extend(
        [
        ("Transcript", ann.get("canonical_transcript")),
        ("Consequence", ann.get("consequence")),
        ("Impact", ann.get("impact")),
        ("Exon", ann.get("exon")),
        ("Intron", ann.get("intron")),
        ("HGVS c.", ann.get("hgvsc")),
        ("HGVS p.", ann.get("hgvsp")),
        ]
    )
    if constraint:
        c_tx = constraint.get("transcript")
        if c_tx and c_tx != ann.get("canonical_transcript"):
            fields.append(("Constr. tx", c_tx))
        fields.extend(
            [
                ("pLI", _format_num(constraint.get("pli"))),
                ("LOEUF", _format_num(constraint.get("lof_oe"))),
                ("mis_z", _format_num(constraint.get("mis_z"))),
            ]
        )
    fields.extend(
        [
            ("CADD Phred", _pdf_format_optional(pred.get("cadd_phred"))),
            ("REVEL max", _pdf_format_optional(pred.get("revel_max"))),
            ("SpliceAI max", _pdf_format_optional(pred.get("spliceai_ds_max"))),
            ("SIFT max", _pdf_format_optional(pred.get("sift_max"))),
            ("PolyPhen max", _pdf_format_optional(pred.get("polyphen_max"))),
            ("phyloP", _pdf_format_optional(pred.get("phylop"))),
        ]
    )
    return fields


def _pdf_kv_block(pdf: Any, fields: list[tuple[str, Any]], content_w: float) -> None:
    pdf.set_font("Helvetica", "B", 8)
    label_w = max(
        _PDF_KV_LABEL_W,
        max(pdf.get_string_width(_pdf_safe(f"{label}:")) + 1.5 for label, _ in fields),
    )
    label_w = min(label_w, 32)
    for label, value in fields:
        _pdf_kv_inline(pdf, label, value, content_w, label_w=label_w)


def _pdf_table_row(
    pdf: Any,
    widths: list[float],
    values: list[Any],
    *,
    bold: bool = False,
    font_size: int = _PDF_TABLE_FONT_SIZE,
    border: str = _PDF_HRULE,
) -> None:
    _pdf_table_row_wrapped(
        pdf,
        widths,
        values,
        bold=bold,
        font_size=font_size,
        line_h=_PDF_TABLE_LINE_H,
        border=border,
    )


def _pdf_freq_table(
    pdf: Any,
    gnomad_block: dict[str, Any] | None,
    content_w: float,
    *,
    title: str,
    section_no: int,
    start_x: float | None = None,
) -> float:
    x0 = pdf.l_margin if start_x is None else start_x
    _pdf_table_section_title(pdf, section_no, title, content_w, start_x=x0)
    if not gnomad_block:
        pdf.set_font("Helvetica", "", 9)
        pdf.set_x(x0)
        pdf.cell(content_w, 5, _pdf_safe("Not observed in this callset."), new_x="LMARGIN", new_y="NEXT")
        return pdf.get_y()

    col_w = list(_PDF_FREQ_COL_W)
    rel_af_w = col_w[4]
    chart_x = x0 + sum(col_w[:-1])

    overall = gnomad_block.get("overall") or {}
    populations = dict(gnomad_block.get("populations") or {})
    grpmax = (gnomad_block.get("grpmax") or "").upper()
    rows: list[tuple[str, dict[str, Any]]] = [("Total", overall)]
    for pop in sorted(populations.keys(), key=lambda p: populations[p].get("af") or 0, reverse=True):
        rows.append((pop, populations[pop]))

    header_h = 5
    pdf.set_fill_color(*_PDF_TABLE_HEAD_FILL)
    pdf.set_text_color(255, 255, 255)
    pdf.set_font("Helvetica", "B", _PDF_TABLE_FONT_SIZE)
    pdf.set_x(x0)
    for width, label in zip(col_w, ["Pop", "AC", "AN", "AF", "rel. AF"]):
        pdf.cell(width, header_h, _pdf_safe(label), border=_PDF_HRULE, fill=True, align="L")
    pdf.ln(header_h)
    pdf.set_text_color(0, 0, 0)

    bar_w = rel_af_w - _PDF_FREQ_CHART_PAD_R
    max_af = max((row[1].get("af") or 0) for row in rows) or 1e-12
    for pop_name, stats in rows:
        af = stats.get("af") or 0
        is_grpmax = pop_name.upper() == grpmax and pop_name != "Total"
        is_total = pop_name == "Total"
        y_row = pdf.get_y()
        row_h = _pdf_table_row_wrapped_at(
            pdf,
            x0,
            col_w,
            [
                pop_name,
                _format_num(stats.get("ac")),
                _format_num(stats.get("an")),
                _format_af(af),
                " ",
            ],
            bold=is_total,
            font_size=_PDF_TABLE_FONT_SIZE,
            line_h=_PDF_TABLE_LINE_H,
        )
        pad = min(1.0, row_h * 0.15)
        _pdf_draw_h_bar(
            pdf,
            chart_x,
            y_row + pad,
            bar_w,
            max(2.5, row_h - 2 * pad),
            af / max_af,
            highlight=is_grpmax,
            total=is_total,
        )

    filt = ", ".join(gnomad_block.get("filter") or []) or "-"
    meta_parts = [f"Filter: {filt}", f"Hom alt: {_format_num(overall.get('nhomalt', 0))}"]
    non_ukb = gnomad_block.get("non_ukb") or {}
    if any(non_ukb.get(k) is not None for k in ("ac", "an", "af")):
        meta_parts.append(
            f"non-UKB AC={_format_num(non_ukb.get('ac'))} AN={_format_num(non_ukb.get('an'))}"
            f" AF={_format_af(non_ukb.get('af'))}"
        )
    meta_line = " | ".join(meta_parts)
    src = gnomad_block.get("source_vcf")
    source_line = f"source: {Path(src).name}" if src else None

    note_h = _PDF_TABLE_LINE_H - 1
    pdf.set_x(x0)
    pdf.set_font("Helvetica", "", _PDF_TABLE_FONT_SIZE - 1)
    pdf.multi_cell(content_w, note_h, _pdf_safe(meta_line), align="L")
    if source_line:
        pdf.set_x(x0)
        pdf.set_font("Helvetica", "", _PDF_TABLE_FONT_SIZE - 2)
        pdf.multi_cell(content_w, note_h - 0.5, _pdf_safe(source_line), align="L")
    pdf.set_font("Helvetica", "", _PDF_TABLE_FONT_SIZE)
    return pdf.get_y()


def _pdf_freq_tables_side_by_side(
    pdf: Any,
    exomes: dict[str, Any] | None,
    genomes: dict[str, Any] | None,
    content_w: float,
    *,
    section_no_ex: int,
    section_no_gen: int,
) -> None:
    gap = _PDF_FREQ_SIDE_GAP
    half_w = (content_w - gap) / 2
    x_left = pdf.l_margin
    x_right = pdf.l_margin + half_w + gap
    y0 = pdf.get_y()
    auto_break = pdf.auto_page_break
    break_margin = pdf.b_margin
    pdf.set_auto_page_break(auto=False)
    try:
        y_left = _pdf_freq_table(
            pdf,
            exomes,
            half_w,
            title="gnomAD EXOMES - POPULATION FREQUENCY",
            section_no=section_no_ex,
            start_x=x_left,
        )
        pdf.set_y(y0)
        y_right = _pdf_freq_table(
            pdf,
            genomes,
            half_w,
            title="gnomAD GENOMES - POPULATION FREQUENCY",
            section_no=section_no_gen,
            start_x=x_right,
        )
        pdf.set_y(max(y_left, y_right))
    finally:
        pdf.set_auto_page_break(auto_break, margin=break_margin)
    _pdf_table_tail_spacer(pdf)


def _pdf_callset_summary(
    pdf: Any,
    exomes: dict[str, Any] | None,
    genomes: dict[str, Any] | None,
    content_w: float,
    *,
    section_no: int,
) -> None:
    _pdf_table_section_title(pdf, section_no, "gnomAD CALLSET SUMMARY (v4.1)", content_w)
    col_w = list(_PDF_CALLSET_COL_W) + [content_w - sum(_PDF_CALLSET_COL_W)]
    _pdf_table_header_row(pdf, col_w, ["Callset", "AC", "AN", "AF", "Filter"])
    for label, block in (("Exomes", exomes), ("Genomes", genomes)):
        if block:
            o = block.get("overall") or {}
            _pdf_table_row(
                pdf,
                col_w,
                [label, _format_num(o.get("ac")), _format_num(o.get("an")), _format_num(o.get("af")), ",".join(block.get("filter") or [])],
            )
        else:
            _pdf_table_row(pdf, col_w, [label, "-", "-", "-", "not observed"])
    _pdf_table_tail_spacer(pdf)


def _pdf_transcript_table(
    pdf: Any,
    transcripts: list[dict[str, Any]],
    content_w: float,
    *,
    section_no: int,
) -> None:
    _pdf_table_section_title(pdf, section_no, "TRANSCRIPT CONSEQUENCES (VEP)", content_w)
    if not transcripts:
        pdf.set_font("Helvetica", "", 9)
        pdf.cell(content_w, 5, _pdf_safe("No VEP transcript entries."), new_x="LMARGIN", new_y="NEXT")
        _pdf_table_tail_spacer(pdf)
        return
    col_w = [18, 26, 46, 16, content_w - 106]
    _pdf_table_header_row(
        pdf,
        col_w,
        ["Gene", "Transcript", "Consequence", "Impact", "HGVS c."],
    )
    for tx in transcripts:
        cons = (tx.get("consequence") or "-").replace("&", ", ")
        hgvsc = tx.get("hgvsc") or "n/a"
        _pdf_table_row_wrapped(
            pdf,
            col_w,
            [
                tx.get("gene_symbol") or "n/a",
                tx.get("transcript") or "n/a",
                cons,
                tx.get("impact") or "n/a",
                hgvsc,
            ],
            font_size=_PDF_TABLE_FONT_SIZE,
            line_h=_PDF_TABLE_LINE_H,
        )
    _pdf_table_tail_spacer(pdf)


def _pdf_display_path(path: str | None) -> str | None:
    if not path:
        return None
    text = str(path)
    prefix = "/mnt/data2/0_database/"
    if text.startswith(prefix):
        return text[len(prefix) :]
    return text


def _clinvar_fields(clinvar: dict[str, Any] | None) -> list[tuple[str, Any]]:
    if not clinvar:
        return [("Status", "No ClinVar record for this allele")]
    return [
        ("Significance", ", ".join(clinvar.get("clinsig") or []) or "n/a"),
        ("Review status", ", ".join(clinvar.get("review_status") or []) or "n/a"),
        ("Condition", " | ".join(clinvar.get("condition") or []) or "n/a"),
        ("HGVS g.", clinvar.get("hgvs")),
        ("Submission", ", ".join(clinvar.get("submission") or []) or "n/a"),
    ]


def _pdf_chrom_ideogram(pdf: Any, chrom: str, pos: int, content_w: float) -> None:
    try:
        from _pdf_ideogram import chrom_ideogram_label, render_chrom_ideogram_png

        png = render_chrom_ideogram_png(chrom, pos)
        label = chrom_ideogram_label(chrom, pos)
    except (KeyError, ValueError, OSError):
        return
    pdf.ln(_PDF_HALF_LINE)
    block_w = content_w * _PDF_IDEOGRAM_WIDTH_FRAC
    img_h = 4.0
    pdf.set_font("Helvetica", "", 10)
    pdf.set_text_color(0, 0, 0)
    label_text = _pdf_safe(label)
    label_w = pdf.get_string_width(label_text)
    gap = _PDF_IDEOGRAM_LABEL_GAP
    img_w = max(20.0, block_w - label_w - gap)
    y0 = pdf.get_y()
    pdf.set_xy(pdf.l_margin, y0)
    pdf.cell(label_w, img_h, label_text, align="L")
    pdf.image(BytesIO(png), x=pdf.l_margin + label_w + gap, y=y0, w=img_w, h=img_h)
    pdf.set_y(y0 + img_h + _PDF_HALF_LINE * 2)


def _pdf_clinvar_section(pdf: Any, clinvar: dict[str, Any] | None, content_w: float) -> None:
    _pdf_subsection_title(pdf, _pdf_titles().subsection_clinvar, content_w)
    fields = _clinvar_fields(clinvar)
    pdf.set_font("Helvetica", "B", 8)
    label_w = max(
        _PDF_KV_LABEL_W,
        max(pdf.get_string_width(_pdf_safe(f"{label}:")) + 1.5 for label, _ in fields),
    )
    label_w = min(label_w, 32)
    for label, value in fields:
        value_color = None
        if label == "Significance" and clinvar:
            tier = clinvar.get("clinical_tier") or "none"
            value_color = PDF_CLINICAL_COLORS.get(tier, (80, 80, 80))
        _pdf_kv_inline(pdf, label, value, content_w, label_w=label_w, value_color=value_color)


def _interpretation_for_pdf(
    interpretation: dict[str, Any],
    *,
    gnomad_exomes: dict[str, Any] | None,
    gnomad_genomes: dict[str, Any] | None,
    gnomad: dict[str, Any] | None,
    clinvar: dict[str, Any] | None,
    sample: dict[str, Any] | None,
    comparison: dict[str, Any] | None,
) -> dict[str, Any]:
    if interpretation.get("flag_details"):
        return interpretation
    flags = interpretation.get("flags") or []
    if not flags:
        return interpretation
    return {
        **interpretation,
        "flag_details": build_flag_details(
            flags,
            gnomad_exomes=gnomad_exomes,
            gnomad_genomes=gnomad_genomes,
            gnomad=gnomad,
            clinvar=clinvar,
            sample=sample,
            comparison=comparison,
        ),
    }


def _pdf_interpretation_flag_line(
    pdf: Any,
    flag: str,
    explanation: str,
    content_w: float,
    *,
    indent: float,
    flag_line_h: float = _PDF_INTERPRETATION_FLAG_LINE_H,
    rule_line_h: float = _PDF_INTERPRETATION_RULE_LINE_H,
) -> None:
    x0 = pdf.l_margin + indent
    body_w = content_w - indent
    flag_text = _pdf_safe(flag)
    expl_text = _pdf_safe(explanation) if explanation else ""

    pdf.set_x(x0)
    pdf.set_text_color(0, 0, 0)
    pdf.set_font("Helvetica", "B", 8)
    pdf.cell(body_w, flag_line_h, flag_text, new_x="LMARGIN", new_y="NEXT")

    if not expl_text:
        return

    pdf.set_x(x0)
    pdf.set_font("Helvetica", "I", 8)
    pdf.set_text_color(90, 90, 90)
    pdf.multi_cell(body_w, rule_line_h, expl_text)
    pdf.set_text_color(0, 0, 0)


def _pdf_interpretation_section(pdf: Any, interpretation: dict[str, Any], content_w: float) -> None:
    _pdf_subsection_title(pdf, _pdf_titles().subsection_interpretation, content_w)
    flag_details = interpretation.get("flag_details") or []
    if not flag_details:
        flags = interpretation.get("flags") or []
        if not flags:
            pdf.set_x(pdf.l_margin)
            pdf.set_font("Helvetica", "", 8)
            pdf.cell(content_w, 5, _pdf_safe("No interpretation flags."), new_x="LMARGIN", new_y="NEXT")
            return
        flag_details = [{"flag": flag, "explanation": ""} for flag in flags]

    pdf.set_x(pdf.l_margin)
    pdf.set_font("Helvetica", "B", 8)
    pdf.cell(content_w, 5, _pdf_safe("Flags:"), new_x="LMARGIN", new_y="NEXT")
    indent = 6
    for idx, item in enumerate(flag_details):
        if idx > 0:
            pdf.ln(_PDF_HALF_LINE / 2)
        _pdf_interpretation_flag_line(
            pdf,
            item.get("flag") or "",
            item.get("explanation") or "",
            content_w,
            indent=indent,
        )


def _pdf_data_sources(pdf: Any, report: dict[str, Any], content_w: float) -> None:
    _pdf_orange_section(pdf, content_w, _pdf_titles().section_data_sources)
    ds = report.get("data_sources") or {}
    gnomad = report.get("gnomad") or {}
    gnomad_exomes = report.get("gnomad_exomes") or {}
    gnomad_genomes = report.get("gnomad_genomes") or {}
    paths = [
        ("Sample VCF", ds.get("sample_vcf_path"), None),
        ("Exomes VCF", ds.get("gnomad_exomes_vcf_path"), gnomad_exomes.get("source_vcf") or gnomad.get("source_vcf")),
        ("Genomes VCF", ds.get("gnomad_genomes_vcf_path"), gnomad_genomes.get("source_vcf")),
        ("Exomes TSV", ds.get("gnomad_exomes_tsv_path"), None),
        ("Genomes TSV", ds.get("gnomad_genomes_tsv_path"), None),
        ("ClinVar VCF", ds.get("clinvar_vcf_path"), None),
        ("Constraint TSV", ds.get("constraint_tsv_path"), None),
    ]
    line_h = _PDF_DATA_SOURCES_LINE_H
    pdf.set_font("Helvetica", "", _PDF_DATA_SOURCES_FONT_SIZE)
    pdf.set_text_color(90, 90, 90)
    for label, path, source_file in paths:
        if not path:
            continue
        if label == "Sample VCF" and not ds.get("sample_vcf"):
            continue
        display = _pdf_display_path(path)
        line = f"{label}: {display}"
        if source_file:
            line += f" ({Path(source_file).name})"
        pdf.cell(content_w, line_h, _pdf_safe(line), new_x="LMARGIN", new_y="NEXT")
    notes = ds.get("notes")
    if notes:
        pdf.multi_cell(content_w, line_h, _pdf_safe(notes))
    pdf.set_text_color(0, 0, 0)


def _pdf_glossary_section(pdf: Any, content_w: float, report: dict[str, Any] | None = None) -> None:
    from _variant_report_i18n import glossary_definition

    _pdf_orange_section(pdf, content_w, _pdf_titles().section_glossary)
    line_h = _PDF_GLOSSARY_LINE_H
    locale = _pdf_locale()
    glossary_entries = list(_PDF_ANNOTATION_GLOSSARY)
    if report and _is_sample_report(report):
        glossary_entries.extend(_PDF_SAMPLE_GLOSSARY)
    for term, definition in glossary_entries:
        if term == "Genotype":
            pdf.ln(_PDF_HALF_LINE)
        definition_text = glossary_definition(term, definition, locale)
        pdf.set_x(pdf.l_margin)
        pdf.set_font("Helvetica", "B", _PDF_GLOSSARY_FONT_SIZE)
        term_prefix = _pdf_safe(f"{term} — ")
        prefix_w = pdf.get_string_width(term_prefix)
        pdf.cell(prefix_w, line_h, term_prefix, new_x="RIGHT", new_y="TOP")
        if locale == "zh":
            pdf.set_font("NotoCJK", "", _PDF_GLOSSARY_FONT_SIZE)
            body_text = _pdf_safe_cjk(definition_text)
        else:
            pdf.set_font("Helvetica", "", _PDF_GLOSSARY_FONT_SIZE)
            body_text = _pdf_safe(definition_text)
        pdf.set_text_color(60, 60, 60)
        pdf.multi_cell(content_w - prefix_w, line_h, body_text, align="L")
        pdf.ln(0.36)
    pdf.set_text_color(0, 0, 0)


def render_variant_report_pdf(
    report: dict[str, Any],
    pdf_path: Path,
    json_hash: str | None = None,
    *,
    locale: str = "en",
) -> None:
    pdf = _VariantReportPDF.create()
    _pdf_bind_locale(pdf, locale)
    from _gene_info import format_hero_gene, resolve_gene_info

    titles = _pdf_titles()
    pdf.set_auto_page_break(auto=True, margin=_PDF_PAGE_FOOTER_H + 2)
    pdf.alias_nb_pages()
    pdf.add_page()
    margin = 12
    pdf.set_left_margin(margin)
    pdf.set_right_margin(margin)
    content_w = pdf.w - 2 * margin

    variant = report["variant"]
    gnomad = report.get("gnomad") or {}
    gnomad_exomes = report.get("gnomad_exomes") or (gnomad if gnomad.get("callset", "").startswith("exomes") else None)
    gnomad_genomes = report.get("gnomad_genomes")
    clinvar = report.get("clinvar") or {}
    pipeline = report.get("pipeline") or {}
    is_reference = report.get("report_type") == "reference_variant"
    is_sample = _is_sample_report(report)
    sample = report.get("sample") or {}
    comparison = report.get("comparison") or {}
    if is_reference:
        interpretation = build_reference_interpretation(gnomad_exomes, gnomad_genomes, clinvar)
    elif is_sample:
        interpretation = build_interpretation(sample, comparison, gnomad, clinvar)
    else:
        interpretation = report.get("interpretation") or {}
    interpretation = _interpretation_for_pdf(
        interpretation,
        gnomad_exomes=gnomad_exomes,
        gnomad_genomes=gnomad_genomes,
        gnomad=gnomad,
        clinvar=clinvar or None,
        sample=sample or None,
        comparison=comparison or None,
    )

    variant_id = report.get("variant_id") or variant.get("variant_id")
    meta = f"{pipeline.get('reference_genome', 'GRCh38')} | gnomAD v4.1"

    # Header
    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(content_w, 7, _pdf_safe("HUMAN VARIANT REPORT"), new_x="LMARGIN", new_y="NEXT")
    pdf.ln(2.5)
    if is_sample:
        _pdf_header_id_row(pdf, f"Sample {sample.get('sample_id')}", meta, content_w, left_size=24)
        if variant_id:
            pdf.set_font("Helvetica", "B", 11)
            pdf.cell(content_w, 6, _pdf_safe(str(variant_id)), new_x="LMARGIN", new_y="NEXT")
    else:
        _pdf_header_id_row(pdf, str(variant_id or variant.get("variant_id")), meta, content_w, left_size=24)
    if report.get("gnomad_variant_page"):
        _pdf_text_right(pdf, report["gnomad_variant_page"], content_w, height=4, font_size=6)
        pdf.ln(5)

    _pdf_orange_section(pdf, content_w, titles.section_variant_overview)

    pdf.set_text_color(0, 0, 0)

    pdf.set_font("Helvetica", "B", 12)
    rsid = variant.get("rsid")
    gene_symbol = (gnomad.get("annotation") or {}).get("gene_symbol")
    gene_info = resolve_gene_info(gene_symbol, report.get("gene_info"))
    hero_gene = format_hero_gene(gene_info) or gene_symbol
    hero_parts = [p for p in (variant.get("locus"), variant.get("change"), rsid, hero_gene) if p]
    hero = "|  " + "  |  ".join(hero_parts) + "  |" if hero_parts else ""
    pdf.cell(content_w, 7, _pdf_safe(hero), new_x="LMARGIN", new_y="NEXT")

    chrom_key = variant.get("chrom")
    variant_pos = variant.get("pos")
    if chrom_key and variant_pos is not None:
        _pdf_chrom_ideogram(pdf, str(chrom_key), int(variant_pos), content_w)

    pdf.set_font("Helvetica", "B", 10)
    if is_sample:
        _pdf_sample_call_line(pdf, sample, comparison, content_w)
        in_exomes = bool(comparison.get("in_gnomad_exomes", comparison.get("in_gnomad")))
        in_genomes = bool(comparison.get("in_gnomad_genomes"))
        in_clinvar = bool(comparison.get("in_clinvar"))
    else:
        in_exomes = gnomad_exomes is not None
        in_genomes = gnomad_genomes is not None
        in_clinvar = bool(clinvar)
    if is_reference or is_sample:
        _pdf_presence_status_line(
            pdf,
            content_w,
            in_exomes=in_exomes,
            in_genomes=in_genomes,
            in_clinvar=in_clinvar,
        )
    pdf.set_text_color(0, 0, 0)
    _pdf_variant_brief(pdf, build_variant_brief_text(report), content_w)

    table_section_no = 0

    if is_sample:
        _pdf_sample_baseline_comparison_section(pdf, sample, comparison, content_w)

    has_annotation = _has_gnomad_annotation(gnomad) or bool(gene_info)
    if has_annotation or clinvar or interpretation:
        _pdf_orange_section(pdf, content_w, titles.section_annotation_interpretation)
        if has_annotation:
            _pdf_subsection_title(pdf, titles.subsection_annotation_scores, content_w)
            _pdf_kv_block(pdf, _annotation_score_fields(gnomad, gene_info), content_w)
        if clinvar:
            _pdf_clinvar_section(pdf, clinvar, content_w)
        if interpretation:
            _pdf_interpretation_section(pdf, interpretation, content_w)

    if _has_gnomad_frequency_data(gnomad_exomes, gnomad_genomes):
        table_section_no += 1
        _pdf_ensure_vertical_space(
            pdf,
            _pdf_gnomad_callset_freq_height_estimate(gnomad_exomes, gnomad_genomes),
        )
        _pdf_orange_section(pdf, content_w, titles.section_frequency_transcripts)
        _pdf_callset_summary(pdf, gnomad_exomes, gnomad_genomes, content_w, section_no=table_section_no)
        _pdf_freq_tables_side_by_side(
            pdf,
            gnomad_exomes,
            gnomad_genomes,
            content_w,
            section_no_ex=table_section_no + 1,
            section_no_gen=table_section_no + 2,
        )
        table_section_no += 2

    transcripts = gnomad.get("transcript_consequences") or []
    if transcripts:
        table_section_no += 1
        _pdf_transcript_table(
            pdf,
            transcripts,
            content_w,
            section_no=table_section_no,
        )

    _pdf_data_sources(pdf, report, content_w)
    _pdf_glossary_section(pdf, content_w, report)

    digest = (json_hash or "")[:8]
    doc_footer = f"Generated {report.get('generated_at')} | {pipeline.get('name')} | JSON sha256 {digest}"
    _pdf_place_doc_footer(pdf, doc_footer, content_w)

    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    pdf.output(str(pdf_path))
    _PDF_RENDER_CTX["locale"] = "en"
    _PDF_RENDER_CTX["titles"] = None


def parse_variant_spec(spec: str) -> tuple[str, int, str, str]:
    """Parse chr21:25881776:AAAAC>A, chr21:25881776:AAAAC:A, or chr21:25881776 AAAAC>A."""
    spec = spec.strip().strip("'\"")
    if ">" in spec:
        left, alt = spec.rsplit(">", 1)
        if ":" in left:
            ref = left.rsplit(":", 1)[-1]
            chrom_pos = left[: -len(ref) - 1]
            chrom, pos_s = chrom_pos.rsplit(":", 1)
            return normalize_chrom(chrom), int(pos_s), ref, alt
    if spec.count(":") >= 3:
        chrom, pos_s, ref, alt = spec.split(":", 3)
        return normalize_chrom(chrom), int(pos_s), ref, alt
    match = re.match(r"^(chr[\w]+):(\d+)\s+(\S+)>(\S+)$", spec)
    if match:
        return normalize_chrom(match.group(1)), int(match.group(2)), match.group(3), match.group(4)
    raise ValueError(f"invalid variant spec: {spec!r} (quote shell args, e.g. 'chr21:25881776:AAAAC>A')")


def append_report_index(index_path: Path, rows: list[dict[str, Any]]) -> None:
    import pandas as pd

    index_path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    if index_path.exists():
        existing = pd.read_csv(index_path)
        frame = pd.concat([existing, frame], ignore_index=True)
        frame = frame.drop_duplicates(subset=["report_stem"], keep="last")
    frame.to_csv(index_path, index=False)


def _read_compare_variants(compare_dir: Path) -> "pd.DataFrame":
    import pandas as pd

    parquet_path = compare_dir / "compare_variants.parquet"
    csv_path = compare_dir / "compare_variants.csv"
    if parquet_path.exists():
        return pd.read_parquet(parquet_path)
    if csv_path.exists():
        return pd.read_csv(csv_path)
    raise FileNotFoundError(f"missing compare variants table under {compare_dir}")


def _compare_bool(series: "pd.Series") -> "pd.Series":
    import pandas as pd

    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin({"true", "yes", "1"})


def select_variants_from_compare(
    compare_dir: Path,
    flagged_only: bool = True,
    limit: int | None = None,
) -> list[dict[str, Any]]:
    df = _read_compare_variants(compare_dir)
    if flagged_only:
        is_novel = _compare_bool(df["is_novel"])
        mask = (
            is_novel
            | df["clinical_tier"].isin(["pathogenic", "likely_pathogenic", "vus"])
            | (df["priority_score"].astype(float) >= 30)
        )
        df = df[mask]
    df = df.sort_values(["priority_score", "chrom", "pos"], ascending=[False, True, True])
    if limit is not None:
        df = df.head(limit)
    return df.to_dict(orient="records")
