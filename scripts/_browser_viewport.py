#!/usr/bin/env python3
"""Viewport payload builder for the interactive baseline genome browser."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from _browser_clinvar import load_clinvar_track
from _browser_gff3 import genes_in_window
from _browser_sample import list_browser_samples, load_sample_track_bins
from _config import CLICK_JUMP_BP, DEFAULT_BROWSER_SAMPLE, ROOT, resolution_label


def pick_gnomad_resolution(window_bp: int) -> int:
    if window_bp > 50_000_000:
        return 1_000_000
    if window_bp > 500_000:
        return 100_000
    return 10_000


def pick_clinvar_resolution(window_bp: int) -> int:
    if window_bp > 50_000:
        return 100_000
    return 10_000


def _auto_ruler_step(span: int, track_width_px: float = 1048.0, min_label_px: float = 88.0) -> int:
    nice_steps = [
        1_000,
        2_000,
        5_000,
        10_000,
        20_000,
        50_000,
        100_000,
        200_000,
        500_000,
        1_000_000,
        2_000_000,
        5_000_000,
        10_000_000,
    ]
    for step in nice_steps:
        if track_width_px * step / max(1, span) >= min_label_px:
            return step
    return nice_steps[-1]


def _format_ruler_label(pos: int, span: int) -> str:
    if span >= 10_000_000:
        if pos % 1_000_000 == 0:
            return f"{pos // 1_000_000} Mb"
        return f"{pos / 1_000_000:.1f} Mb"
    if pos >= 1_000_000:
        return f"{pos / 1_000_000:.2f} Mb"
    if span >= 50_000:
        return f"{pos / 1_000:.0f} kb"
    return f"{pos:,}"


def build_ruler_ticks(start: int, end: int, track_width_px: float = 1048.0) -> dict[str, Any]:
    span = end - start + 1
    step = _auto_ruler_step(span, track_width_px=track_width_px)
    tick_start = ((start + step - 1) // step) * step
    positions = list(range(tick_start, end + 1, step))
    if not positions or positions[0] != start:
        positions.insert(0, start)
    if positions[-1] != end:
        positions.append(end)

    min_label_gap = 72.0
    last_label_x = -1e9
    ticks: list[dict[str, Any]] = []
    for i, pos in enumerate(positions):
        x = (pos - start) / max(1, span) * track_width_px
        show_label = i == 0 or i == len(positions) - 1 or (x - last_label_x >= min_label_gap)
        tick: dict[str, Any] = {"pos": pos}
        if show_label:
            tick["label"] = _format_ruler_label(pos, span)
            last_label_x = x
        ticks.append(tick)
    return {"step_bp": step, "ticks": ticks}


def slice_histogram_bins(
    bins_df: pd.DataFrame,
    start: int,
    end: int,
    value_col: str,
) -> list[dict[str, int]]:
    if bins_df.empty:
        return []
    sub = bins_df[(bins_df["bin_start"] < end) & (bins_df["bin_end"] > start)]
    extra_cols = ("snv_count", "indel_count", "rare_count", "common_count")
    rows: list[dict[str, int]] = []
    for _, row in sub.iterrows():
        value = int(row[value_col])
        if value <= 0:
            continue
        item: dict[str, int] = {
            "start": int(max(row["bin_start"], start)),
            "end": int(min(row["bin_end"], end)),
            "value": value,
        }
        for col in extra_cols:
            if col in row.index:
                item[col] = int(row[col])
        rows.append(item)
    return rows


def filter_cytobands(
    bands: list[dict[str, Any]],
    start: int,
    end: int,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for band in bands:
        if int(band["end"]) < start or int(band["start"]) > end:
            continue
        out.append(
            {
                "start": int(band["start"]),
                "end": int(band["end"]),
                "name": str(band["name"]),
                "stain": str(band.get("stain", "")),
                "fill": str(band.get("fill", "#EEEEEE")),
            }
        )
    return out


class BrowserDataStore:
    """Load and cache browser manifest, cytobands, genes, and bin tables."""

    def __init__(self, root: Path = ROOT) -> None:
        self.root = root
        self.browser_dir = root / "processed" / "browser"
        self._manifest: dict[str, Any] | None = None
        self._cytobands: dict[str, list[dict[str, Any]]] | None = None
        self._chrom_index: dict[str, dict[str, Any]] | None = None

    def manifest(self) -> dict[str, Any]:
        if self._manifest is None:
            path = self.browser_dir / "manifest.json"
            if not path.exists():
                raise FileNotFoundError(f"missing browser manifest: {path}")
            self._manifest = json.loads(path.read_text(encoding="utf-8"))
            self._chrom_index = {item["chrom"]: item for item in self._manifest["chromosomes"]}
        return self._manifest

    def chrom_index(self) -> dict[str, dict[str, Any]]:
        self.manifest()
        assert self._chrom_index is not None
        return self._chrom_index

    def chrom_meta(self, chrom: str) -> dict[str, Any]:
        index = self.chrom_index()
        if chrom not in index:
            raise KeyError(f"unknown chrom: {chrom}")
        return index[chrom]

    def cytobands_for(self, chrom: str) -> list[dict[str, Any]]:
        if self._cytobands is None:
            path = self.browser_dir / "cytobands.json"
            self._cytobands = json.loads(path.read_text(encoding="utf-8"))
        return self._cytobands.get(chrom, [])

    @lru_cache(maxsize=48)
    def load_genes(self, chrom: str) -> tuple[dict[str, Any], ...]:
        meta = self.chrom_meta(chrom)
        path = self.root / meta["genes_parquet"]
        if not path.exists():
            return tuple()
        frame = pd.read_parquet(path)
        return tuple(frame.to_dict(orient="records"))

    @lru_cache(maxsize=96)
    def load_bins(self, chrom: str, resolution: int) -> pd.DataFrame:
        meta = self.chrom_meta(chrom)
        label = resolution_label(resolution)
        rel = meta["bins"][label]
        path = self.root / rel
        if not path.exists():
            return pd.DataFrame()
        return pd.read_parquet(path)

    def validate_viewport(self, chrom: str, start: int, end: int) -> tuple[int, int, int]:
        meta = self.chrom_meta(chrom)
        chrom_length = int(meta["length"])
        start_i = max(1, int(start))
        end_i = min(chrom_length, int(end))
        if start_i >= end_i:
            raise ValueError(f"invalid viewport: start={start_i} end={end_i}")
        return start_i, end_i, chrom_length

    def build_viewport(
        self,
        chrom: str,
        start: int,
        end: int,
        sample_id: str | None = None,
    ) -> dict[str, Any]:
        start_i, end_i, chrom_length = self.validate_viewport(chrom, start, end)
        manifest = self.manifest()
        window_bp = end_i - start_i + 1

        gnomad_res = pick_gnomad_resolution(window_bp)
        clinvar_res = pick_clinvar_resolution(window_bp)
        gnomad_bins = self.load_bins(chrom, gnomad_res)
        clinvar_data = load_clinvar_track(chrom, start_i, end_i, clinvar_res)

        gene_rows = genes_in_window(list(self.load_genes(chrom)), start_i, end_i)
        gene_features = [
            {
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": row["strand"],
                "biotype": row.get("biotype", "unknown"),
                "transcript_id": row.get("transcript_id"),
                "transcript_name": row.get("transcript_name"),
                "exon_count": int(row.get("exon_count", 0)),
                "exons": row.get("exons", []),
            }
            for row in gene_rows
        ]

        track_defs = manifest["tracks"]
        track_order = manifest["track_order"]
        tracks: list[dict[str, Any]] = []

        for track_id in track_order:
            tdef = track_defs[track_id]
            if track_id == "genes":
                tracks.append(
                    {
                        "id": "genes",
                        "label": tdef["label"],
                        "sublabel": tdef["sublabel"],
                        "type": "gene",
                        "feature_count": len(gene_features),
                        "features": gene_features,
                    }
                )
            elif track_id == "gnomad":
                value_col = tdef["value_col"]
                tracks.append(
                    {
                        "id": "gnomad",
                        "label": tdef["label"],
                        "sublabel": tdef["sublabel"],
                        "type": "histogram",
                        "resolution_bp": gnomad_res,
                        "fill": tdef["fill"],
                        "value_col": value_col,
                        "bin_count": len(gnomad_bins),
                        "bins": slice_histogram_bins(gnomad_bins, start_i, end_i, value_col),
                    }
                )
            elif track_id == "clinvar":
                tracks.append(
                    {
                        "id": "clinvar",
                        "label": tdef["label"],
                        "sublabel": tdef["sublabel"],
                        "type": "clinvar_stacked",
                        "resolution_bp": clinvar_res,
                        "bin_count": len(clinvar_data["bins"]),
                        "bins": clinvar_data["bins"],
                    }
                )
            elif track_id == "sample":
                if not sample_id:
                    continue
                sample_data = load_sample_track_bins(
                    sample_id, chrom, start_i, end_i, clinvar_res, self.root
                )
                label = tdef.get("label", "Sample").replace("{sample_id}", sample_id)
                sublabel = tdef.get("sublabel", "known · novel").replace("{sample_id}", sample_id)
                tracks.append(
                    {
                        "id": "sample",
                        "label": label,
                        "sublabel": sublabel,
                        "type": "sample_stacked",
                        "sample_id": sample_id,
                        "resolution_bp": clinvar_res,
                        "variant_count": sample_data["variant_count"],
                        "novel_count": sample_data["novel_count"],
                        "bin_count": len(sample_data["bins"]),
                        "bins": sample_data["bins"],
                    }
                )

        payload_sample = sample_id
        return {
            "chrom": chrom,
            "start": start_i,
            "end": end_i,
            "window_bp": window_bp,
            "chrom_length": chrom_length,
            "assembly": manifest.get("assembly", "GRCh38"),
            "click_jump_bp": int(manifest.get("click_jump_bp", CLICK_JUMP_BP)),
            "zoom_steps_bp": manifest.get("zoom_steps_bp", []),
            "resolutions": {
                "gnomad_bp": gnomad_res,
                "clinvar_bp": clinvar_res,
            },
            "cytobands": filter_cytobands(self.cytobands_for(chrom), start_i, end_i),
            "ideogram_cytobands": [
                {
                    "start": int(b["start"]),
                    "end": int(b["end"]),
                    "name": str(b["name"]),
                    "stain": str(b.get("stain", "")),
                    "fill": str(b.get("fill", "#EEEEEE")),
                }
                for b in self.cytobands_for(chrom)
            ],
            "ruler": build_ruler_ticks(start_i, end_i),
            "tracks": tracks,
            "sample_id": payload_sample,
        }

    def sample_bin_variants(
        self,
        chrom: str,
        sample_id: str,
        bin_start: int,
        bin_end: int,
        *,
        limit: int = 100,
    ) -> dict[str, Any]:
        from _browser_sample import load_sample_bin_variant_details

        self.chrom_meta(chrom)
        bin_start_i = int(bin_start)
        bin_end_i = int(bin_end)
        if bin_start_i > bin_end_i:
            raise ValueError("bin_start must be <= bin_end")
        return load_sample_bin_variant_details(
            sample_id,
            chrom,
            bin_start_i,
            bin_end_i,
            self.root,
            limit=limit,
        )

    def chrom_info(self, chrom: str) -> dict[str, Any]:
        meta = self.chrom_meta(chrom)
        manifest = self.manifest()
        chrom_length = int(meta["length"])
        return {
            "chrom": chrom,
            "length": chrom_length,
            "variant_count": int(meta["variant_count"]),
            "default_viewport": {"start": 1, "end": chrom_length},
            "click_jump_bp": int(manifest.get("click_jump_bp", CLICK_JUMP_BP)),
            "zoom_steps_bp": [chrom_length, *manifest.get("zoom_steps_bp", [])],
            "track_order": manifest["track_order"],
            "tracks": manifest["tracks"],
            "cytoband_count": len(self.cytobands_for(chrom)),
            "gene_count": len(self.load_genes(chrom)),
        }

    def api_manifest(self) -> dict[str, Any]:
        manifest = self.manifest()
        samples = manifest.get("samples") or list_browser_samples(self.root)
        default_sample = manifest.get("default_sample", DEFAULT_BROWSER_SAMPLE)
        return {
            "assembly": manifest.get("assembly", "GRCh38"),
            "gnomad": manifest.get("gnomad"),
            "gnomad_callset": manifest.get("gnomad_callset"),
            "clinvar": manifest.get("clinvar"),
            "default_chrom": manifest.get("default_chrom", "chr21"),
            "default_sample": default_sample,
            "samples": samples,
            "click_jump_bp": manifest.get("click_jump_bp", CLICK_JUMP_BP),
            "zoom_steps_bp": manifest.get("zoom_steps_bp", []),
            "track_order": manifest["track_order"],
            "tracks": manifest["tracks"],
            "chromosomes": [
                {
                    "chrom": item["chrom"],
                    "length": int(item["length"]),
                    "variant_count": int(item["variant_count"]),
                }
                for item in manifest["chromosomes"]
            ],
        }
