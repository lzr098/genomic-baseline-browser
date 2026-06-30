"""Compare individual sample variants against gnomAD exome + ClinVar baseline."""

from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pysam

from _bin_builder import bin_key, clinical_tier, load_clinvar_map
from _config import (
    ALL_CHROMSOMES,
    BINS,
    CLINVAR_VCF,
    DATA_RELEASE,
    gnomad_vcf_path,
)


def classify_match(in_gnomad: bool, in_clinvar: bool) -> str:
    if in_gnomad and in_clinvar:
        return "known_gnomad_clinvar"
    if in_gnomad:
        return "known_gnomad"
    if in_clinvar:
        return "known_clinvar_only"
    return "novel_in_sample"


def priority_score(row: dict[str, Any]) -> float:
    if not row.get("is_novel"):
        return 0.0
    score = 20.0
    tier = row.get("clinical_tier") or "none"
    if tier in {"pathogenic", "likely_pathogenic"}:
        score += 80.0
    elif tier == "vus":
        score += 30.0
    elif tier == "conflicting":
        score += 15.0
    af = row.get("gnomad_af")
    if af is not None and af < 0.001:
        score += 10.0
    if row.get("sample_zygosity") == "hom_alt":
        score += 5.0
    return score


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


def _clinvar_meta(clinvar_map: dict[tuple[int, str, str], dict[str, Any]], row: dict[str, Any]) -> dict[str, Any]:
    meta = clinvar_map.get((row["pos"], row["ref"], row["alt"]), {})
    clinsig = meta.get("clinsig")
    return {
        "clinsig": clinsig,
        "clinical_tier": clinical_tier(clinsig),
        "rsid": meta.get("rsid"),
    }


def compare_chrom_variants(
    chrom: str,
    sample_rows: list[dict[str, Any]],
    sample_id: str,
    clinvar_map: dict[tuple[int, str, str], dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    if not sample_rows:
        return []

    gnomad_path = gnomad_vcf_path(chrom)
    if not gnomad_path.exists():
        raise FileNotFoundError(f"missing gnomAD VCF: {gnomad_path}")

    clinvar_map = clinvar_map if clinvar_map is not None else load_clinvar_map(chrom)
    compared: list[dict[str, Any]] = []

    with pysam.VariantFile(str(gnomad_path)) as gnomad_vcf:
        contig = _resolve_gnomad_contig(gnomad_vcf, chrom)
        for row in sample_rows:
            gnomad_af: float | None = None
            in_gnomad = False
            for rec in gnomad_vcf.fetch(contig, row["pos"] - 1, row["pos"]):
                if rec.pos != row["pos"] or rec.ref != row["ref"]:
                    continue
                if row["alt"] not in (rec.alts or []):
                    continue
                in_gnomad = True
                gnomad_af = _gnomad_af_for_allele(rec, row["alt"])
                break

            clin = _clinvar_meta(clinvar_map, row)
            in_clinvar = (row["pos"], row["ref"], row["alt"]) in clinvar_map
            match_status = classify_match(in_gnomad, in_clinvar)
            in_baseline = in_gnomad or in_clinvar
            is_novel = match_status == "novel_in_sample"

            compared.append(
                {
                    "sample_id": sample_id,
                    "chrom": chrom,
                    "pos": row["pos"],
                    "ref": row["ref"],
                    "alt": row["alt"],
                    "variant_id": row["variant_id"],
                    "sample_gt": row.get("sample_gt"),
                    "sample_zygosity": row.get("sample_zygosity"),
                    "sample_dp": row.get("sample_dp"),
                    "match_status": match_status,
                    "is_novel": is_novel,
                    "in_baseline": in_baseline,
                    "gnomad_af": gnomad_af,
                    "clinsig": clin.get("clinsig"),
                    "clinical_tier": clin.get("clinical_tier"),
                    "rsid": clin.get("rsid"),
                    "gene_symbol": None,
                    "consequence": None,
                    "impact": None,
                }
            )

    for row in compared:
        row["priority_score"] = priority_score(row)
    return compared


def build_compare_bins(
    variants_df: pd.DataFrame,
    resolution: int = 100_000,
) -> pd.DataFrame:
    from _config import resolution_label

    baseline_path = BINS / f"genome_bins_{resolution_label(resolution)}.parquet"
    if not baseline_path.exists():
        raise FileNotFoundError(f"missing baseline bins: {baseline_path}")

    baseline = pd.read_parquet(baseline_path)
    baseline = baseline[baseline["resolution"] == resolution].copy()

    counters: dict[tuple[str, int], dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    for _, row in variants_df.iterrows():
        bk = bin_key(int(row["pos"]), resolution)
        key = (str(row["chrom"]), bk)
        bucket = counters[key]
        bucket["sample_variant_count"] += 1
        if row["is_novel"]:
            bucket["novel_count"] += 1
        else:
            bucket["known_count"] += 1
        tier = row.get("clinical_tier") or "none"
        if tier in {"pathogenic", "likely_pathogenic"}:
            bucket["pathogenic_count"] += 1
        elif tier == "vus":
            bucket["vus_count"] += 1

    rows: list[dict[str, Any]] = []
    for (chrom, bin_start), bucket in sorted(counters.items()):
        bin_end = bin_start + resolution - 1
        baseline_row = baseline[
            (baseline["chrom"] == chrom) & (baseline["bin_start"] == bin_start)
        ]
        baseline_count = (
            int(baseline_row["variant_count"].iloc[0]) if not baseline_row.empty else 0
        )
        novel_count = int(bucket["novel_count"])
        rows.append(
            {
                "chrom": chrom,
                "bin_start": bin_start,
                "bin_end": bin_end,
                "resolution": resolution,
                "sample_variant_count": int(bucket["sample_variant_count"]),
                "novel_count": novel_count,
                "known_count": int(bucket["known_count"]),
                "pathogenic_count": int(bucket["pathogenic_count"]),
                "vus_count": int(bucket["vus_count"]),
                "baseline_variant_count": baseline_count,
                "novel_enrichment": novel_count / max(baseline_count, 1),
            }
        )
    return pd.DataFrame(rows)


def build_compare_by_gene(variants_df: pd.DataFrame) -> pd.DataFrame:
    """Placeholder gene rollup until GFF3 annotation is wired in Phase 4."""
    if variants_df.empty:
        return pd.DataFrame(
            columns=[
                "gene_symbol",
                "chrom",
                "sample_variant_count",
                "novel_count",
                "known_count",
                "pathogenic_count",
                "lof_count",
                "missense_count",
                "has_novel_in_exon",
                "max_priority_score",
            ]
        )

    grouped = variants_df.groupby("chrom", dropna=False)
    rows: list[dict[str, Any]] = []
    for chrom, frame in grouped:
        rows.append(
            {
                "gene_symbol": None,
                "chrom": chrom,
                "sample_variant_count": int(len(frame)),
                "novel_count": int(frame["is_novel"].sum()),
                "known_count": int((~frame["is_novel"]).sum()),
                "pathogenic_count": int(
                    frame["clinical_tier"].isin({"pathogenic", "likely_pathogenic"}).sum()
                ),
                "lof_count": 0,
                "missense_count": 0,
                "has_novel_in_exon": False,
                "max_priority_score": float(frame["priority_score"].max()),
            }
        )
    return pd.DataFrame(rows)


def build_compare_summary(
    sample_id: str,
    sample_vcf: Path,
    variants_df: pd.DataFrame,
    started_at: datetime,
    finished_at: datetime,
) -> dict[str, Any]:
    by_status = variants_df["match_status"].value_counts().to_dict() if not variants_df.empty else {}
    by_chrom = variants_df.groupby("chrom").size().to_dict() if not variants_df.empty else {}
    by_tier = (
        variants_df["clinical_tier"].fillna("none").value_counts().to_dict()
        if not variants_df.empty
        else {}
    )
    return {
        **DATA_RELEASE,
        "sample_id": sample_id,
        "sample_vcf": str(sample_vcf),
        "started_at": started_at.isoformat(),
        "finished_at": finished_at.isoformat(),
        "duration_sec": round((finished_at - started_at).total_seconds(), 1),
        "counts": {
            "total_variants": int(len(variants_df)),
            "novel_variants": int(variants_df["is_novel"].sum()) if not variants_df.empty else 0,
            "known_variants": int((~variants_df["is_novel"]).sum()) if not variants_df.empty else 0,
            "by_match_status": {str(k): int(v) for k, v in by_status.items()},
            "by_chrom": {str(k): int(v) for k, v in by_chrom.items()},
            "by_clinical_tier": {str(k): int(v) for k, v in by_tier.items()},
        },
        "reference": {
            "gnomad_dir": str(gnomad_vcf_path("chr1").parent),
            "clinvar_vcf": str(CLINVAR_VCF),
        },
    }


def _format_variant_site(row: pd.Series) -> str:
    return f"{int(row['pos'])}:{row['ref']}>{row['alt']}"


def _format_variant_locus(row: pd.Series) -> str:
    return f"{row['chrom']}:{int(row['pos'])} {row['ref']}>{row['alt']}"


def _variants_for_csv(df: pd.DataFrame) -> pd.DataFrame:
    """Select and order columns for human-readable CSV export."""
    out = df.copy()
    out["variant_change"] = out["ref"].astype(str) + ">" + out["alt"].astype(str)
    out["genomic_locus"] = (
        out["chrom"].astype(str) + ":" + out["pos"].astype(int).astype(str)
    )
    columns = [
        "sample_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "variant_change",
        "genomic_locus",
        "match_status",
        "is_novel",
        "in_baseline",
        "sample_gt",
        "sample_zygosity",
        "sample_dp",
        "gnomad_af",
        "clinsig",
        "clinical_tier",
        "rsid",
        "priority_score",
        "variant_id",
    ]
    for col in columns:
        if col not in out.columns:
            out[col] = None
    out = out[columns]
    if "is_novel" in out.columns:
        out["is_novel"] = out["is_novel"].map({True: "yes", False: "no"})
    if "in_baseline" in out.columns:
        out["in_baseline"] = out["in_baseline"].map({True: "yes", False: "no"})
    return out.sort_values(["chrom", "pos", "ref", "alt"]).reset_index(drop=True)


def _join_sites(frame: pd.DataFrame) -> str:
    if frame.empty:
        return ""
    return "; ".join(_format_variant_site(row) for _, row in frame.sort_values("pos").iterrows())


def _chrom_summary_for_csv(variants_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for chrom, frame in variants_df.groupby("chrom", sort=False):
        known = frame[~frame["is_novel"]]
        novel = frame[frame["is_novel"]]
        rows.append(
            {
                "sample_id": sample_id,
                "chrom": chrom,
                "total_variants": int(len(frame)),
                "known_variants": int(len(known)),
                "novel_variants": int(len(novel)),
                "pathogenic_or_likely": int(
                    frame["clinical_tier"].isin({"pathogenic", "likely_pathogenic"}).sum()
                ),
                "vus": int((frame["clinical_tier"] == "vus").sum()),
                "max_priority_score": float(frame["priority_score"].max()),
                "known_variant_sites": _join_sites(known),
                "novel_variant_sites": _join_sites(novel),
            }
        )
    from _bin_builder import chrom_sort_key

    return pd.DataFrame(rows).sort_values("chrom", key=lambda s: s.map(chrom_sort_key))


def _bins_with_variant_sites(
    bins_df: pd.DataFrame,
    variants_df: pd.DataFrame,
    resolution: int = 100_000,
) -> pd.DataFrame:
    """Attach pos:ref>alt site lists to each occupied 100 kb bin."""
    if bins_df.empty:
        return bins_df

    site_map: dict[tuple[str, int], dict[str, list[pd.Series]]] = {}
    for _, row in variants_df.iterrows():
        bk = bin_key(int(row["pos"]), resolution)
        key = (str(row["chrom"]), bk)
        site_map.setdefault(key, {"all": [], "known": [], "novel": []})
        site_map[key]["all"].append(row)
        bucket = "novel" if row["is_novel"] else "known"
        site_map[key][bucket].append(row)

    out = bins_df.copy()
    all_sites: list[str] = []
    known_sites: list[str] = []
    novel_sites: list[str] = []
    for _, brow in out.iterrows():
        key = (str(brow["chrom"]), int(brow["bin_start"]))
        buckets = site_map.get(key)
        if not buckets:
            all_sites.append("")
            known_sites.append("")
            novel_sites.append("")
            continue
        all_df = pd.DataFrame(buckets["all"])
        known_df = pd.DataFrame(buckets["known"]) if buckets["known"] else pd.DataFrame()
        novel_df = pd.DataFrame(buckets["novel"]) if buckets["novel"] else pd.DataFrame()
        all_sites.append(_join_sites(all_df))
        known_sites.append(_join_sites(known_df))
        novel_sites.append(_join_sites(novel_df))

    out["variant_sites"] = all_sites
    out["known_variant_sites"] = known_sites
    out["novel_variant_sites"] = novel_sites
    return out


def _variants_long_with_bins(
    variants_df: pd.DataFrame,
    bins_df: pd.DataFrame,
    sample_id: str,
    resolution: int = 100_000,
) -> pd.DataFrame:
    """One row per variant with bin context — primary diff detail export."""
    rows: list[dict[str, Any]] = []
    bin_lookup = {
        (str(r["chrom"]), int(r["bin_start"])): r for _, r in bins_df.iterrows()
    }
    for _, row in _variants_for_csv(variants_df).iterrows():
        bk = bin_key(int(row["pos"]), resolution)
        bin_row = bin_lookup.get((str(row["chrom"]), bk), {})
        rows.append(
            {
                "sample_id": sample_id,
                "chrom": row["chrom"],
                "pos": int(row["pos"]),
                "ref": row["ref"],
                "alt": row["alt"],
                "variant_change": row["variant_change"],
                "genomic_locus": row["genomic_locus"],
                "match_status": row["match_status"],
                "is_novel": row["is_novel"],
                "in_baseline": row["in_baseline"],
                "sample_gt": row["sample_gt"],
                "sample_zygosity": row["sample_zygosity"],
                "sample_dp": row["sample_dp"],
                "gnomad_af": row["gnomad_af"],
                "clinsig": row["clinsig"],
                "clinical_tier": row["clinical_tier"],
                "priority_score": row["priority_score"],
                "bin_start": bk,
                "bin_end": int(bin_row.get("bin_end", bk + resolution - 1))
                if len(bin_row)
                else bk + resolution - 1,
                "baseline_variant_count": bin_row.get("baseline_variant_count", ""),
                "bin_novel_count": bin_row.get("novel_count", ""),
            }
        )
    return pd.DataFrame(rows)


def export_compare_csv(
    output_dir: Path,
    sample_id: str,
    variants_df: pd.DataFrame,
    novel_df: pd.DataFrame,
    bins_df: pd.DataFrame,
    gene_df: pd.DataFrame,
) -> dict[str, str]:
    """Write spreadsheet-friendly CSV exports alongside parquet outputs."""
    paths: dict[str, str] = {}

    diff_csv = output_dir / "compare_diff.csv"
    _variants_long_with_bins(variants_df, bins_df, sample_id).to_csv(
        diff_csv, index=False, encoding="utf-8-sig"
    )
    paths["compare_diff_csv"] = str(diff_csv)

    variants_csv = output_dir / "compare_variants.csv"
    _variants_for_csv(variants_df).to_csv(variants_csv, index=False, encoding="utf-8-sig")
    paths["compare_variants_csv"] = str(variants_csv)

    novel_csv = output_dir / "novel_candidates.csv"
    _variants_for_csv(novel_df).to_csv(novel_csv, index=False, encoding="utf-8-sig")
    paths["novel_candidates_csv"] = str(novel_csv)

    chrom_csv = output_dir / "compare_by_chrom.csv"
    _chrom_summary_for_csv(variants_df, sample_id).to_csv(
        chrom_csv, index=False, encoding="utf-8-sig"
    )
    paths["compare_by_chrom_csv"] = str(chrom_csv)

    if not bins_df.empty:
        bins_csv = output_dir / "compare_bins_100kb.csv"
        bins_out = _bins_with_variant_sites(bins_df, variants_df).sort_values(
            ["novel_count", "novel_enrichment"], ascending=False
        ).reset_index(drop=True)
        bins_out.to_csv(bins_csv, index=False, encoding="utf-8-sig")
        paths["compare_bins_100kb_csv"] = str(bins_csv)

    if not gene_df.empty:
        gene_csv = output_dir / "compare_by_gene.csv"
        gene_df.to_csv(gene_csv, index=False, encoding="utf-8-sig")
        paths["compare_by_gene_csv"] = str(gene_csv)

    return paths


def export_compare_csv_from_dir(compare_dir: Path, sample_id: str | None = None) -> dict[str, str]:
    """Regenerate CSV exports from existing compare parquet outputs."""
    import json

    summary_path = compare_dir / "compare_summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"missing {summary_path}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    sample_id = sample_id or str(summary.get("sample_id", compare_dir.name))

    variants_df = pd.read_parquet(compare_dir / "compare_variants.parquet")
    novel_df = pd.read_parquet(compare_dir / "novel_candidates.parquet")
    bins_df = pd.read_parquet(compare_dir / "compare_bins_100kb.parquet")
    gene_df = pd.read_parquet(compare_dir / "compare_by_gene.parquet")
    return export_compare_csv(compare_dir, sample_id, variants_df, novel_df, bins_df, gene_df)


def compare_sample(
    sample_vcf: Path,
    sample_id: str,
    output_dir: Path,
    chromosomes: tuple[str, ...] | None = None,
) -> dict[str, Any]:
    from _bin_builder import chrom_sort_key
    from _sample_loader import load_sample_variants_by_chrom
    import json

    chromosomes = chromosomes or ALL_CHROMSOMES
    started = datetime.now(timezone.utc)
    output_dir.mkdir(parents=True, exist_ok=True)

    by_chrom = load_sample_variants_by_chrom(sample_vcf, sample_name=sample_id)
    all_rows: list[dict[str, Any]] = []
    for chrom in sorted(chromosomes, key=chrom_sort_key):
        rows = by_chrom.get(chrom, [])
        if not rows:
            continue
        clinvar_map = load_clinvar_map(chrom)
        all_rows.extend(compare_chrom_variants(chrom, rows, sample_id, clinvar_map))

    variants_df = pd.DataFrame(all_rows)
    if variants_df.empty:
        raise ValueError(f"no variants found in sample VCF: {sample_vcf}")

    novel_df = (
        variants_df[variants_df["is_novel"]]
        .sort_values("priority_score", ascending=False)
        .reset_index(drop=True)
    )
    bins_df = build_compare_bins(variants_df)
    gene_df = build_compare_by_gene(variants_df)
    finished = datetime.now(timezone.utc)
    summary = build_compare_summary(sample_id, sample_vcf, variants_df, started, finished)

    variants_path = output_dir / "compare_variants.parquet"
    novel_path = output_dir / "novel_candidates.parquet"
    bins_path = output_dir / "compare_bins_100kb.parquet"
    gene_path = output_dir / "compare_by_gene.parquet"
    summary_path = output_dir / "compare_summary.json"

    variants_df.to_parquet(variants_path, index=False)
    novel_df.to_parquet(novel_path, index=False)
    bins_df.to_parquet(bins_path, index=False)
    gene_df.to_parquet(gene_path, index=False)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    csv_paths = export_compare_csv(
        output_dir, sample_id, variants_df, novel_df, bins_df, gene_df
    )

    return {
        "summary": summary,
        "paths": {
            "compare_variants": str(variants_path),
            "novel_candidates": str(novel_path),
            "compare_bins_100kb": str(bins_path),
            "compare_by_gene": str(gene_path),
            "compare_summary": str(summary_path),
            **csv_paths,
        },
    }
