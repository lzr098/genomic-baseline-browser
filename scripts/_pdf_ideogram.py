"""UCSC ideogram track PNG for variant report PDFs (label rendered separately in PDF)."""

from __future__ import annotations

import html
import json
from io import BytesIO
from typing import Any

import cairosvg

from _config import BROWSER

IDEOGRAM_VIEW_WIDTH = 1000
IDEOGRAM_HEIGHT = 14

ACEN_P_FILL = "#3A3A3A"
ACEN_Q_FILL = "#707070"
MARKER_STROKE = "#CC0000"
BAND_STROKE = "#000000"
BAND_STROKE_W = 1.0
MARKER_STROKE_W = 3.2
MARKER_EXTEND = 4.5

_manifest: dict[str, Any] | None = None
_cytobands: dict[str, list[dict[str, Any]]] | None = None


def normalize_chrom(chrom: str) -> str:
    text = str(chrom).strip()
    return text if text.startswith("chr") else f"chr{text}"


def _load_manifest() -> dict[str, Any]:
    global _manifest
    if _manifest is None:
        _manifest = json.loads((BROWSER / "manifest.json").read_text(encoding="utf-8"))
    return _manifest


def _load_cytobands() -> dict[str, list[dict[str, Any]]]:
    global _cytobands
    if _cytobands is None:
        _cytobands = json.loads((BROWSER / "cytobands.json").read_text(encoding="utf-8"))
    return _cytobands


def chrom_length(chrom: str) -> int:
    chrom = normalize_chrom(chrom)
    for item in _load_manifest().get("chromosomes", []):
        if item.get("chrom") == chrom:
            return int(item["length"])
    bands = bands_for_chrom(chrom)
    if not bands:
        raise KeyError(f"unknown chromosome: {chrom}")
    return int(bands[-1]["end"])


def bands_for_chrom(chrom: str) -> list[dict[str, Any]]:
    chrom = normalize_chrom(chrom)
    return list(_load_cytobands().get(chrom, []))


def cytoband_at(bands: list[dict[str, Any]], pos: int) -> str:
    for band in bands:
        if pos >= int(band["start"]) and pos < int(band["end"]):
            return str(band.get("name") or "")
    return str(bands[-1].get("name") or "") if bands else ""


def chrom_ideogram_label(chrom: str, pos: int) -> str:
    chrom = normalize_chrom(chrom)
    bands = bands_for_chrom(chrom)
    length = chrom_length(chrom)
    pos = max(1, min(int(pos), length))
    band_name = cytoband_at(bands, pos)
    return f"{chrom} ({band_name})" if band_name else chrom


def ideogram_x(pos: int | float, chrom_length: int, track_width: int = IDEOGRAM_VIEW_WIDTH) -> float:
    length = max(1, int(chrom_length))
    return (float(pos) / length) * track_width


def _track_x(pos: int | float, chrom_length: int) -> float:
    return ideogram_x(pos, chrom_length, IDEOGRAM_VIEW_WIDTH)


def _append_ideogram_outlines(
    parts: list[str],
    bands: list[dict[str, Any]],
    acen_set: set[int],
    acen_indices: list[int],
    junction_bp: int | None,
    chrom_len: int,
    height: float,
    bar_y: float,
) -> None:
    y_top = bar_y
    y_bot = bar_y + height
    mid_y = bar_y + height / 2.0
    line = (
        f'stroke="{BAND_STROKE}" stroke-width="{BAND_STROKE_W}" fill="none" '
        f'shape-rendering="crispEdges" vector-effect="non-scaling-stroke"'
    )

    for i, band in enumerate(bands):
        if i in acen_set:
            continue

        x0 = _track_x(int(band["start"]), chrom_len)
        x1 = _track_x(int(band["end"]), chrom_len)
        before_acen = len(acen_indices) >= 2 and i == acen_indices[0] - 1
        after_acen = len(acen_indices) >= 2 and i == acen_indices[-1] + 1

        if before_acen and junction_bp is not None:
            left = bands[acen_indices[0]]
            right = bands[acen_indices[-1]]
            x0p = _track_x(int(left["start"]), chrom_len)
            x1q = _track_x(int(right["end"]), chrom_len)
            cx = _track_x(junction_bp, chrom_len)

            parts.append(f'<line x1="{x0:.2f}" y1="{y_top}" x2="{x0p:.2f}" y2="{y_top}" {line}/>')
            parts.append(f'<line x1="{x0:.2f}" y1="{y_bot}" x2="{x0p:.2f}" y2="{y_bot}" {line}/>')
            if x0 > 0.5:
                parts.append(f'<line x1="{x0:.2f}" y1="{y_top}" x2="{x0:.2f}" y2="{y_bot}" {line}/>')
            parts.append(
                f'<polyline points="{x0p:.2f},{y_top} {cx:.2f},{mid_y:.2f} {x0p:.2f},{y_bot}" {line}/>'
            )
            parts.append(
                f'<polyline points="{cx:.2f},{mid_y:.2f} {x1q:.2f},{y_top} {x1q:.2f},{y_bot}" {line}/>'
            )
            continue

        if after_acen:
            parts.append(f'<line x1="{x0:.2f}" y1="{y_top}" x2="{x1:.2f}" y2="{y_top}" {line}/>')
            parts.append(f'<line x1="{x0:.2f}" y1="{y_bot}" x2="{x1:.2f}" y2="{y_bot}" {line}/>')
            parts.append(f'<line x1="{x0:.2f}" y1="{y_top}" x2="{x0:.2f}" y2="{y_bot}" {line}/>')
            continue

        parts.append(f'<line x1="{x0:.2f}" y1="{y_top}" x2="{x1:.2f}" y2="{y_top}" {line}/>')
        parts.append(f'<line x1="{x0:.2f}" y1="{y_bot}" x2="{x1:.2f}" y2="{y_bot}" {line}/>')
        if x0 > 0.5:
            parts.append(f'<line x1="{x0:.2f}" y1="{y_top}" x2="{x0:.2f}" y2="{y_bot}" {line}/>')

    x_left = _track_x(0, chrom_len)
    x_right = _track_x(chrom_len, chrom_len)
    parts.append(f'<line x1="{x_left:.2f}" y1="{y_top}" x2="{x_left:.2f}" y2="{y_bot}" {line}/>')
    parts.append(f'<line x1="{x_right:.2f}" y1="{y_top}" x2="{x_right:.2f}" y2="{y_bot}" {line}/>')


def render_ideogram_svg(
    bands: list[dict[str, Any]],
    chrom_len: int,
    marker_pos: int,
) -> str:
    bar_h = IDEOGRAM_HEIGHT
    bar_y = MARKER_EXTEND
    svg_h = bar_h + 2 * MARKER_EXTEND
    length = max(1, int(chrom_len))
    total_w = IDEOGRAM_VIEW_WIDTH

    parts: list[str] = []
    acen_indices = [i for i, band in enumerate(bands) if band.get("stain") == "acen"]
    acen_set = set(acen_indices)
    junction_bp = int(bands[acen_indices[0]]["end"]) if len(acen_indices) >= 2 else None

    stroke_attrs = (
        f'stroke="{BAND_STROKE}" stroke-width="{BAND_STROKE_W}" '
        f'vector-effect="non-scaling-stroke" shape-rendering="crispEdges"'
    )

    for band in bands:
        if band.get("stain") == "acen":
            continue
        x0 = _track_x(int(band["start"]), length)
        x1 = _track_x(int(band["end"]), length)
        bw = max(0.25, x1 - x0)
        fill = html.escape(str(band.get("fill") or "#FFFFFF"))
        parts.append(
            f'<rect x="{x0:.2f}" y="{bar_y:.2f}" width="{bw:.2f}" height="{bar_h:.2f}" '
            f'fill="{fill}" {stroke_attrs}/>'
        )

    if len(acen_indices) >= 2:
        left = bands[acen_indices[0]]
        right = bands[acen_indices[-1]]
        x0p = _track_x(int(left["start"]), length)
        x1q = _track_x(int(right["end"]), length)
        cx = _track_x(junction_bp, length)
        mid_y = bar_y + bar_h / 2.0
        parts.append(
            f'<polygon points="{x0p:.2f},{bar_y} {cx:.2f},{mid_y:.2f} {x0p:.2f},{bar_y + bar_h}" '
            f'fill="{ACEN_P_FILL}" {stroke_attrs}/>'
        )
        parts.append(
            f'<polygon points="{cx:.2f},{mid_y:.2f} {x1q:.2f},{bar_y} {x1q:.2f},{bar_y + bar_h}" '
            f'fill="{ACEN_Q_FILL}" {stroke_attrs}/>'
        )

    _append_ideogram_outlines(parts, bands, acen_set, acen_indices, junction_bp, length, bar_h, bar_y)

    marker_x = _track_x(marker_pos, length)
    parts.append(
        f'<line x1="{marker_x:.2f}" y1="0" x2="{marker_x:.2f}" y2="{svg_h:.2f}" '
        f'stroke="{MARKER_STROKE}" stroke-width="{MARKER_STROKE_W}" '
        f'shape-rendering="crispEdges" vector-effect="non-scaling-stroke"/>'
    )

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total_w} {svg_h}" '
        f'preserveAspectRatio="none">{"".join(parts)}</svg>'
    )


def render_chrom_ideogram_png(
    chrom: str,
    pos: int,
    *,
    output_width: int = 1680,
) -> bytes:
    chrom = normalize_chrom(chrom)
    bands = bands_for_chrom(chrom)
    if not bands:
        raise ValueError(f"no cytobands for {chrom}")

    length = chrom_length(chrom)
    pos = max(1, min(int(pos), length))
    svg = render_ideogram_svg(bands, length, pos)
    svg_h = IDEOGRAM_HEIGHT + 2 * MARKER_EXTEND
    output_height = max(32, round(output_width * svg_h / IDEOGRAM_VIEW_WIDTH))
    return cairosvg.svg2png(
        bytestring=svg.encode("utf-8"),
        output_width=output_width,
        output_height=output_height,
    )


def render_chrom_ideogram_png_to(path: str | BytesIO, chrom: str, pos: int, **kwargs: Any) -> None:
    data = render_chrom_ideogram_png(chrom, pos, **kwargs)
    if isinstance(path, BytesIO):
        path.write(data)
        path.seek(0)
    else:
        with open(path, "wb") as fh:
            fh.write(data)
