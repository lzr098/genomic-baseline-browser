"""ClinVar stacked-bin data for the interactive browser."""

from __future__ import annotations

from collections import defaultdict
from typing import Any

import pysam

from _bin_builder import bin_key, clinical_tier
from _config import CLINVAR_VCF, clinvar_contig

DISPLAY_GROUP: dict[str, str] = {
    "pathogenic": "plp",
    "likely_pathogenic": "plp",
    "vus": "vus",
    "conflicting": "vus",
    "benign": "benign",
    "likely_benign": "benign",
}


def display_group(clinical_tier_name: str) -> str | None:
    return DISPLAY_GROUP.get(clinical_tier_name)


def _format_clnsig(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, (tuple, list)):
        return ",".join(str(x) for x in value)
    return str(value)


def load_clinvar_track(chrom: str, start: int, end: int, resolution: int) -> dict[str, Any]:
    """Return stacked ClinVar bins (P/LP, VUS/conflicting, benign) for a viewport."""
    if not CLINVAR_VCF.exists():
        return {"bins": []}

    contig = clinvar_contig(chrom)
    bin_counters: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    with pysam.VariantFile(str(CLINVAR_VCF)) as vcf:
        if contig not in vcf.header.contigs:
            return {"bins": []}

        for rec in vcf.fetch(contig, start - 1, end):
            pos = int(rec.pos)
            if pos < start or pos > end:
                continue
            clnsig = _format_clnsig(rec.info.get("CLNSIG"))
            tier = clinical_tier(clnsig)
            group = display_group(tier)
            if group is None:
                continue

            for _alt in rec.alts or []:
                bk = bin_key(pos, resolution)
                bin_counters[bk][group] += 1

    bins: list[dict[str, int]] = []
    for bk in sorted(bin_counters.keys()):
        counts = bin_counters[bk]
        plp = int(counts.get("plp", 0))
        vus = int(counts.get("vus", 0))
        benign = int(counts.get("benign", 0))
        total = plp + vus + benign
        if total <= 0:
            continue
        bins.append(
            {
                "start": int(max(bk, start)),
                "end": int(min(bk + resolution - 1, end)),
                "plp": plp,
                "vus": vus,
                "benign": benign,
                "total": total,
            }
        )

    return {"bins": bins}
