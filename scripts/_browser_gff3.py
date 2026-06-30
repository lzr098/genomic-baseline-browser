#!/usr/bin/env python3
"""Extract per-chromosome gene models from Ensembl GFF3 for the interactive browser."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

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
