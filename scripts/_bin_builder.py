"""Stream gnomAD exome VCF per chromosome and build multi-resolution bin tables."""

from __future__ import annotations

import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pysam

from _config import (
    ALL_CHROMSOMES,
    BIN_RESOLUTIONS,
    BY_CHROM,
    BINS,
    CLINVAR_VCF,
    DATA_RELEASE,
    GNOMAD_DIR,
    ROOT,
    clinvar_contig,
    gnomad_vcf_path,
    resolution_label,
)


def clinical_tier(clinsig: str | None) -> str:
    if not clinsig:
        return "none"
    s = clinsig.lower()
    if "conflict" in s:
        return "conflicting"
    if "pathogenic" in s and "likely" not in s:
        return "pathogenic"
    if "likely_pathogenic" in s.replace(" ", "_") or "likely pathogenic" in s:
        return "likely_pathogenic"
    if "uncertain" in s or "vus" in s:
        return "vus"
    if "likely_benign" in s.replace(" ", "_") or "likely benign" in s:
        return "likely_benign"
    if "benign" in s:
        return "benign"
    return "other"


def bin_key(pos: int, resolution: int) -> int:
    return ((pos - 1) // resolution) * resolution + 1


def chrom_length_from_vcf(vcf_path: Path, contig: str) -> int:
    with pysam.VariantFile(str(vcf_path)) as vcf:
        if contig not in vcf.header.contigs:
            raise KeyError(f"contig {contig!r} not in {vcf_path}")
        return int(vcf.header.contigs[contig].length)


def load_clinvar_map(
    chrom: str,
    clinvar_vcf: Path = CLINVAR_VCF,
) -> dict[tuple[int, str, str], dict[str, Any]]:
    out: dict[tuple[int, str, str], dict[str, Any]] = {}
    if not clinvar_vcf.exists():
        return out

    contig = clinvar_contig(chrom)
    with pysam.VariantFile(str(clinvar_vcf)) as vcf:
        if contig not in vcf.header.contigs:
            return out
        for rec in vcf.fetch(contig):
            clnsig = rec.info.get("CLNSIG")
            if isinstance(clnsig, (tuple, list)):
                clnsig = ",".join(str(x) for x in clnsig)
            elif clnsig is not None:
                clnsig = str(clnsig)
            for alt in rec.alts or []:
                out[(rec.pos, rec.ref, alt)] = {"clinsig": clnsig}
    return out


def stream_chrom_bins(
    chrom: str,
    vcf_path: Path,
    clinvar_map: dict[tuple[int, str, str], dict[str, Any]],
    chrom_length: int,
    resolutions: tuple[int, ...] | None = None,
) -> tuple[dict[int, list[dict[str, Any]]], int]:
    resolutions = resolutions or BIN_RESOLUTIONS
    counters: dict[int, dict[int, dict[str, int]]] = {
        res: defaultdict(lambda: defaultdict(int)) for res in resolutions
    }
    variant_alleles = 0

    with pysam.VariantFile(str(vcf_path)) as vcf:
        contig = chrom if chrom in vcf.header.contigs else chrom.replace("chr", "")
        fetch_iter = vcf.fetch(contig) if contig in vcf.header.contigs else vcf

        for rec in fetch_iter:
            pos = rec.pos
            ref = rec.ref
            af = rec.info.get("AF")
            if isinstance(af, tuple):
                af = af[0]
            elif af is not None:
                af = float(af)

            for alt in rec.alts or []:
                variant_alleles += 1
                clin = clinvar_map.get((pos, ref, alt), {})
                tier = clinical_tier(clin.get("clinsig"))

                for res in resolutions:
                    bk = bin_key(pos, res)
                    bucket = counters[res][bk]
                    bucket["variant_count"] += 1
                    if len(ref) == 1 and len(alt) == 1:
                        bucket["snv_count"] += 1
                    else:
                        bucket["indel_count"] += 1
                    if af is not None:
                        if af < 0.001:
                            bucket["rare_count"] += 1
                        else:
                            bucket["common_count"] += 1
                    if tier in {"pathogenic", "likely_pathogenic"}:
                        bucket["pathogenic_count"] += 1
                    elif tier == "vus":
                        bucket["vus_count"] += 1

    results: dict[int, list[dict[str, Any]]] = {}
    for res in resolutions:
        rows: list[dict[str, Any]] = []
        for bk in sorted(counters[res].keys()):
            data = counters[res][bk]
            rows.append(
                {
                    "chrom": chrom,
                    "bin_start": bk,
                    "bin_end": min(bk + res - 1, chrom_length),
                    "resolution": res,
                    "variant_count": int(data.get("variant_count", 0)),
                    "snv_count": int(data.get("snv_count", 0)),
                    "indel_count": int(data.get("indel_count", 0)),
                    "rare_count": int(data.get("rare_count", 0)),
                    "common_count": int(data.get("common_count", 0)),
                    "pathogenic_count": int(data.get("pathogenic_count", 0)),
                    "vus_count": int(data.get("vus_count", 0)),
                }
            )
        results[res] = rows
    return results, variant_alleles


def write_chrom_bins(chrom: str, bin_sets: dict[int, list[dict[str, Any]]]) -> dict[str, str]:
    BY_CHROM.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    for res, rows in bin_sets.items():
        label = resolution_label(res)
        out_path = BY_CHROM / f"{chrom}_bins_{label}.parquet"
        pd.DataFrame(rows).to_parquet(out_path, index=False)
        paths[label] = str(out_path.relative_to(ROOT))
    return paths


def build_chrom_bins(chrom: str) -> dict[str, Any]:
    vcf_path = gnomad_vcf_path(chrom)
    if not vcf_path.exists():
        raise FileNotFoundError(f"missing gnomAD VCF: {vcf_path}")

    started = datetime.now(timezone.utc)
    chrom_length = chrom_length_from_vcf(vcf_path, chrom)
    clinvar_map = load_clinvar_map(chrom)
    bin_sets, variant_alleles = stream_chrom_bins(
        chrom, vcf_path, clinvar_map, chrom_length
    )
    output_paths = write_chrom_bins(chrom, bin_sets)
    finished = datetime.now(timezone.utc)

    return {
        "chrom": chrom,
        "chrom_length": chrom_length,
        "variant_allele_count": variant_alleles,
        "clinvar_keys": len(clinvar_map),
        "gnomad_vcf": str(vcf_path),
        "output_paths": output_paths,
        "started_at": started.isoformat(),
        "finished_at": finished.isoformat(),
        "duration_sec": round((finished - started).total_seconds(), 1),
    }


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    if chrom.startswith("chr"):
        body = chrom[3:]
    else:
        body = chrom
    if body.isdigit():
        return (0, f"{int(body):02d}")
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, f"{order.get(body, 99):02d}{body}")


def merge_genome_bins(chromosomes: tuple[str, ...] | None = None) -> dict[str, str]:
    chromosomes = chromosomes or ALL_CHROMSOMES
    BINS.mkdir(parents=True, exist_ok=True)
    merged_paths: dict[str, str] = {}

    for res in BIN_RESOLUTIONS:
        label = resolution_label(res)
        frames: list[pd.DataFrame] = []
        for chrom in sorted(chromosomes, key=chrom_sort_key):
            path = BY_CHROM / f"{chrom}_bins_{label}.parquet"
            if not path.exists():
                raise FileNotFoundError(f"missing per-chrom bin file: {path}")
            frames.append(pd.read_parquet(path))
        genome_df = pd.concat(frames, ignore_index=True)
        genome_df["_chrom_order"] = genome_df["chrom"].map(chrom_sort_key)
        genome_df = genome_df.sort_values(["_chrom_order", "bin_start"]).drop(
            columns="_chrom_order"
        )
        out_path = BINS / f"genome_bins_{label}.parquet"
        genome_df.to_parquet(out_path, index=False)
        merged_paths[label] = str(out_path.relative_to(ROOT))

    return merged_paths


def write_manifest(
    chrom_stats: list[dict[str, Any]],
    genome_paths: dict[str, str] | None = None,
) -> Path:
    genome_paths = genome_paths or {}
    total_variants = sum(item["variant_allele_count"] for item in chrom_stats)
    manifest = {
        **DATA_RELEASE,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "paths": {
            "gnomad_dir": str(GNOMAD_DIR),
            "clinvar_vcf": str(CLINVAR_VCF),
            "bins_dir": str(BINS),
        },
        "resolutions": {
            resolution_label(res): res for res in BIN_RESOLUTIONS
        },
        "totals": {
            "chromosomes": len(chrom_stats),
            "variant_allele_count": total_variants,
        },
        "genome_bins": genome_paths,
        "chromosomes": sorted(chrom_stats, key=lambda x: chrom_sort_key(x["chrom"])),
    }
    manifest_path = ROOT / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest_path
