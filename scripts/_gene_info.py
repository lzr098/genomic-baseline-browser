"""Gene symbol metadata for variant reports (full name + brief function)."""

from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
from functools import lru_cache
from typing import Any

from _config import GFF3

_SUMMARY_CACHE: dict[str, str] = {}


def _brief_summary(text: str, *, max_len: int = 320) -> str:
    text = " ".join(text.split())
    if len(text) <= max_len:
        return text
    clipped = text[:max_len]
    for sep in (". ", "; "):
        idx = clipped.rfind(sep)
        if idx >= 120:
            return clipped[: idx + 1].strip()
    return clipped.rstrip(" ,;") + "…"


@lru_cache(maxsize=1)
def _ensembl_gene_index() -> dict[str, str]:
    index: dict[str, str] = {}
    if not GFF3.exists():
        return index
    with GFF3.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or "\tgene\t" not in line:
                continue
            attrs = line.rstrip("\n").split("\t", 8)[-1]
            name_match = re.search(r"Name=([^;]+)", attrs)
            if not name_match:
                continue
            symbol = name_match.group(1)
            desc_match = re.search(r"description=([^;]+)", attrs)
            full_name = ""
            if desc_match:
                full_name = desc_match.group(1).replace("%3B", ";")
                full_name = re.sub(r"\s*\[Source:.*\]$", "", full_name).strip()
            if full_name:
                index[symbol] = full_name
    return index


def fetch_gene_summary(symbol: str, *, timeout: float = 8.0) -> str | None:
    symbol = symbol.strip()
    if not symbol:
        return None
    if symbol in _SUMMARY_CACHE:
        return _SUMMARY_CACHE[symbol]
    summary: str | None = None
    try:
        query = urllib.parse.urlencode(
            {
                "q": f"symbol:{symbol}",
                "species": "human",
                "fields": "summary",
                "size": "1",
            }
        )
        url = f"https://mygene.info/v3/query?{query}"
        req = urllib.request.Request(url, headers={"User-Agent": "genos-variant-report/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            payload = json.load(resp)
        hits = payload.get("hits") or []
        if hits:
            raw = hits[0].get("summary")
            if isinstance(raw, str) and raw.strip():
                summary = raw.strip()
    except Exception:
        summary = None
    if summary:
        _SUMMARY_CACHE[symbol] = summary
    return summary


def resolve_gene_info(symbol: str | None, cached: dict[str, Any] | None = None) -> dict[str, Any] | None:
    if not symbol:
        return None
    info = dict(cached or {})
    info.setdefault("symbol", symbol)
    if not info.get("full_name"):
        info["full_name"] = _ensembl_gene_index().get(symbol)
    if not info.get("summary"):
        info["summary"] = fetch_gene_summary(symbol)
    if not info.get("full_name") and not info.get("summary"):
        return {"symbol": symbol, "full_name": None, "summary": None}
    return {
        "symbol": symbol,
        "full_name": info.get("full_name") or None,
        "summary": info.get("summary") or None,
    }


def format_hero_gene(gene_info: dict[str, Any] | None) -> str | None:
    if not gene_info:
        return None
    symbol = gene_info.get("symbol")
    if not symbol:
        return None
    full_name = gene_info.get("full_name")
    if full_name:
        return f"{symbol} ({full_name})"
    return str(symbol)


def format_gene_overview_line(gene_info: dict[str, Any] | None) -> str | None:
    if not gene_info or not gene_info.get("summary"):
        return None
    return _brief_summary(str(gene_info["summary"]))


def format_gene_function_value(gene_info: dict[str, Any] | None) -> str | None:
    return format_gene_overview_line(gene_info)


def format_gene_annotation_value(gene_info: dict[str, Any] | None) -> str | None:
    return format_hero_gene(gene_info)
