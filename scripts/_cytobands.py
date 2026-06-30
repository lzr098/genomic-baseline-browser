"""Load hg38 cytoband coordinates for any chromosome."""

from __future__ import annotations

from pathlib import Path

from _chr21_cytobands import STAIN_COLORS, acen_band_fill  # noqa: F401
from _config import CYTOBAND


def load_cytobands(chrom: str, path: Path = CYTOBAND) -> list[dict[str, object]]:
    if not path.exists():
        raise FileNotFoundError(f"missing cytoband file: {path}")
    bands: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            row_chrom, start, end, name, stain = line.rstrip("\n").split("\t")
            if row_chrom != chrom:
                continue
            bands.append(
                {
                    "chrom": row_chrom,
                    "start": int(start),
                    "end": int(end),
                    "name": name,
                    "stain": stain,
                    "fill": STAIN_COLORS.get(stain, "#EEEEEE"),
                }
            )
    return bands


def trim_bands_to_length(bands: list[dict[str, object]], chrom_length: int) -> list[dict[str, object]]:
    if not bands:
        return bands
    out = [dict(b) for b in bands]
    if out[-1]["end"] < chrom_length:
        out[-1]["end"] = chrom_length
    return out
