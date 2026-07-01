#!/usr/bin/env python3
"""Extract per-chromosome gene models from Ensembl GFF3 for the interactive browser."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from _config import ALL_CHROMSOMES, gff_contig

PREFER_TRANSCRIPT_TAGS: tuple[str, ...] = (
    "MANE_Select",
    "Ensembl_canonical",
    "gencode_primary",
)


def _parse_attrs(raw: str) -> dict[str, str]:
    return dict(item.split("=", 1) for item in raw.split(";") if "=" in item)


def _pick_transcript(
    gene_id: str,
    transcripts: dict[str, dict[str, Any]],
    exons_by_tx: dict[str, list[dict[str, int]]],
) -> tuple[dict[str, Any] | None, list[dict[str, int]]]:
    candidates = [tx for tx in transcripts.values() if tx.get("parent_gene_id") == gene_id]
    if not candidates:
        return None, []

    chosen: dict[str, Any] | None = None
    for tag in PREFER_TRANSCRIPT_TAGS:
        for tx in candidates:
            if tag in tx.get("tags", set()):
                chosen = tx
                break
        if chosen:
            break
    if not chosen:
        chosen = max(candidates, key=lambda t: t["end"] - t["start"])

    exons = sorted(exons_by_tx.get(chosen["transcript_id"], []), key=lambda x: x["start"])
    merged: list[dict[str, int]] = []
    for ex in exons:
        if merged and ex["start"] <= merged[-1]["end"] + 1:
            merged[-1]["end"] = max(merged[-1]["end"], ex["end"])
        else:
            merged.append({"start": ex["start"], "end": ex["end"]})
    return chosen, merged


def extract_genes_from_gff3(
    gff3: Path,
    chromosomes: tuple[str, ...] = ALL_CHROMSOMES,
) -> dict[str, list[dict[str, Any]]]:
    """One pass over GFF3; return gene rows keyed by chrom (chr21 style)."""
    target_contigs = {gff_contig(chrom) for chrom in chromosomes}
    contig_to_chrom = {gff_contig(chrom): chrom for chrom in chromosomes}

    genes: dict[str, dict[str, Any]] = {}
    transcripts: dict[str, dict[str, Any]] = {}
    exons_by_tx: dict[str, list[dict[str, int]]] = {}
    gene_ids_by_contig: dict[str, set[str]] = {c: set() for c in target_contigs}

    with gff3.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            contig = parts[0]
            if contig not in target_contigs:
                continue

            feature = parts[2]
            start, end = int(parts[3]), int(parts[4])
            strand = parts[6]
            attr = _parse_attrs(parts[8])

            if feature == "gene":
                gene_id = attr.get("ID", attr.get("gene_id", ""))
                if not gene_id:
                    continue
                gene_ids_by_contig[contig].add(gene_id)
                genes[gene_id] = {
                    "gene_id": gene_id,
                    "gene_name": attr.get("Name", attr.get("gene_id", gene_id)),
                    "chrom": contig_to_chrom[contig],
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "biotype": attr.get("biotype", "unknown"),
                }
                continue

            if feature in {"mRNA", "transcript"}:
                tx_id = attr.get("ID", attr.get("transcript_id", ""))
                parent = attr.get("Parent", "")
                if not tx_id or not parent:
                    continue
                tags = {t for t in attr.get("tag", "").replace(",", ";").split(";") if t}
                transcripts[tx_id] = {
                    "transcript_id": tx_id,
                    "parent_gene_id": parent,
                    "transcript_name": attr.get("Name", tx_id),
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "biotype": attr.get("biotype", "unknown"),
                    "tags": tags,
                }
                continue

            if feature == "exon":
                parent = attr.get("Parent", "")
                if parent in transcripts:
                    exons_by_tx.setdefault(parent, []).append({"start": start, "end": end})

    out: dict[str, list[dict[str, Any]]] = {chrom: [] for chrom in chromosomes}
    for contig, gene_ids in gene_ids_by_contig.items():
        chrom = contig_to_chrom[contig]
        for gene_id in sorted(gene_ids):
            gene = genes.get(gene_id)
            if not gene:
                continue
            tx, exons = _pick_transcript(gene_id, transcripts, exons_by_tx)
            row = {
                **gene,
                "transcript_id": tx["transcript_id"] if tx else None,
                "transcript_name": tx["transcript_name"] if tx else None,
                "exon_count": len(exons),
                "exons_json": json.dumps(exons, separators=(",", ":")),
            }
            out[chrom].append(row)
    return out


def _location_in_transcript(
    pos: int,
    gene: dict[str, Any],
    exons: list[dict[str, int]],
) -> tuple[str, int | None]:
    """Return (location_label, distance_to_nearest_exon_bp) for a position within a gene.

    location_label examples: "exon", "intron", "5' flanking", "3' flanking".
    """
    strand = gene.get("strand") or "+"
    exons_sorted = sorted(exons, key=lambda e: e["start"])
    if not exons_sorted:
        return "genic", None

    first_exon = exons_sorted[0]
    last_exon = exons_sorted[-1]

    # Within an exon?
    for i, ex in enumerate(exons_sorted):
        if ex["start"] <= pos <= ex["end"]:
            return f"exon {i + 1}", 0

    # Upstream / downstream of transcript (using strand-aware labels)
    if strand == "+":
        if pos < first_exon["start"]:
            return "5' flanking", first_exon["start"] - pos
        if pos > last_exon["end"]:
            return "3' flanking", pos - last_exon["end"]
    else:
        if pos < first_exon["start"]:
            return "3' flanking", first_exon["start"] - pos
        if pos > last_exon["end"]:
            return "5' flanking", pos - last_exon["end"]

    # In intron: find nearest exons and report distance
    nearest_dist = None
    for ex in exons_sorted:
        if pos < ex["start"]:
            nearest_dist = ex["start"] - pos
            break
    if nearest_dist is None:
        nearest_dist = pos - last_exon["end"]

    return "intron", nearest_dist


@lru_cache(maxsize=24)
def load_genes_frame(chrom: str, root_s: str) -> pd.DataFrame:
    path = Path(root_s) / "processed" / "browser" / "genes" / f"{chrom}_genes.parquet"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_parquet(path)


def gene_context_for_position(
    chrom: str,
    pos: int,
    root: Path = Path("."),
) -> dict[str, Any] | None:
    """Return the gene/exon context for a single variant position.

    If multiple genes overlap, prefer protein_coding and the smallest gene span.
    """
    frame = load_genes_frame(chrom, str(root))
    if frame.empty:
        return None

    overlapping = frame[(frame["start"] <= pos) & (frame["end"] >= pos)]
    if overlapping.empty:
        return None

    # Prefer protein_coding, then smallest span
    def _score(row: pd.Series) -> tuple[int, int]:
        is_coding = 1 if (row.get("biotype") or "") == "protein_coding" else 0
        span = int(row["end"]) - int(row["start"])
        return (is_coding, -span)

    best = overlapping.iloc[0]
    best_score = _score(best)
    for _, row in overlapping.iterrows():
        s = _score(row)
        if s > best_score:
            best = row
            best_score = s

    exons = json.loads(best["exons_json"]) if best.get("exons_json") else []
    location, distance = _location_in_transcript(pos, best.to_dict(), exons)

    return {
        "gene_name": best.get("gene_name"),
        "gene_id": best.get("gene_id"),
        "transcript_id": best.get("transcript_id"),
        "transcript_name": best.get("transcript_name"),
        "biotype": best.get("biotype"),
        "strand": best.get("strand"),
        "gene_start": int(best["start"]),
        "gene_end": int(best["end"]),
        "location": location,
        "distance_to_exon": distance,
    }


def genes_in_window(
    genes: list[dict[str, Any]],
    start: int,
    end: int,
) -> list[dict[str, Any]]:
    """Filter genes overlapping [start, end] and attach parsed exons."""
    rows: list[dict[str, Any]] = []
    for gene in genes:
        if gene["end"] < start or gene["start"] > end:
            continue
        item = dict(gene)
        item["exons"] = json.loads(gene["exons_json"])
        rows.append(item)
    return rows
