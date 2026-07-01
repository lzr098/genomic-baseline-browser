"""Query gnomAD (and MyVariant.info fallback) for variant allele frequencies."""

from __future__ import annotations

import re
from typing import Any

import requests

from _variant_report import parse_gnomad_variant_id


def _myvariant_id(variant_id: str) -> str:
    """Convert gnomAD variant ID to MyVariant.info HGVS genomic notation.

    Example: 1-55051215-G-GA -> chr1:g.55051215G>A
             21-26034358-C-T -> chr21:g.26034358C>T
    """
    try:
        chrom, pos, ref, alt = parse_gnomad_variant_id(variant_id)
    except ValueError:
        # Fall back to manual parsing if normalize_chrom behaves unexpectedly.
        parts = variant_id.strip().split("-")
        if len(parts) < 4:
            raise ValueError(f"invalid gnomAD variant_id: {variant_id!r}") from None
        chrom = parts[0]
        pos = parts[1]
        ref = parts[2]
        alt = "-".join(parts[3:])
    chrom = chrom.lstrip("chr") if chrom.startswith("chr") else chrom
    return f"chr{chrom}:g.{pos}{ref}>{alt}"


def _extract_af(data: dict[str, Any]) -> dict[str, Any] | None:
    """Extract gnomAD exome/genome AF from MyVariant.info response."""
    exome = data.get("gnomad_exome") or {}
    genome = data.get("gnomad_genome") or {}

    def _get(d: dict[str, Any], *keys: str) -> Any:
        for key in keys:
            if key not in d or d[key] is None:
                continue
            # MyVariant stores values under a nested dict with the same key.
            if isinstance(d[key], dict):
                inner = d[key].get(key)
                if inner is not None:
                    return inner
            else:
                return d[key]
        return None

    exome_af = _get(exome, "af")
    exome_ac = _get(exome, "ac")
    exome_an = _get(exome, "an")
    genome_af = _get(genome, "af")
    genome_ac = _get(genome, "ac")
    genome_an = _get(genome, "an")

    if exome_af is None and genome_af is None:
        return None

    # Prefer exome AF; if only genome is available, use that.
    af = exome_af if exome_af is not None else genome_af

    return {
        "af": af,
        "exome": {
            "af": exome_af,
            "ac": exome_ac,
            "an": exome_an,
        },
        "genome": {
            "af": genome_af,
            "ac": genome_ac,
            "an": genome_an,
        },
        "source": "myvariant.info",
    }


def query_gnomad_af(variant_id: str, timeout: float = 15.0) -> dict[str, Any] | None:
    """Return allele-frequency data for a gnomAD-style variant ID.

    Uses MyVariant.info as the primary source so we avoid gnomAD GraphQL
    rate limits and schema churn. Returns None when the variant is not
    found or the request fails.
    """
    mv_id = _myvariant_id(variant_id)
    url = f"https://myvariant.info/v1/variant/{mv_id}"
    try:
        # Disable proxies to avoid local interceptors (Clash etc.).
        resp = requests.get(url, timeout=timeout, proxies={"http": None, "https": None})
        resp.raise_for_status()
        data = resp.json()
        if not isinstance(data, dict) or data.get("not_found"):
            return None
        return _extract_af(data)
    except Exception:
        return None
