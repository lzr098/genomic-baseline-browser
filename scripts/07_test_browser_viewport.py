#!/usr/bin/env python3
"""Smoke-test browser viewport API without starting uvicorn."""

from __future__ import annotations

import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from _browser_reports import generate_reference_report_pdf  # noqa: E402
from _browser_viewport import BrowserDataStore  # noqa: E402


def main() -> None:
    store = BrowserDataStore()
    manifest = store.api_manifest()
    print("manifest chromosomes:", len(manifest["chromosomes"]))
    print("track_order:", manifest["track_order"])

    vp_full = store.build_viewport("chr21", 1, 46_709_983)
    print("chr21 full genes:", vp_full["tracks"][0]["feature_count"])
    print("chr21 full gnomAD bins:", len(vp_full["tracks"][1]["bins"]))

    vp_app = store.build_viewport("chr21", 25_800_000, 26_200_000)
    genes = vp_app["tracks"][0]["features"]
    app_gene = next((g for g in genes if g["gene_name"] == "APP"), None)
    print("APP window:", app_gene["gene_name"] if app_gene else "missing", "exons:", len(app_gene["exons"]) if app_gene else 0)

    vp_100k = store.build_viewport("chr21", 25_850_000, 25_950_000)
    print("100kb window resolutions:", vp_100k["resolutions"])
    print("100kb gnomAD bins:", len(vp_100k["tracks"][1]["bins"]))

    vp_clcn = store.build_viewport("chr1", 16_050_000, 16_053_000, sample_id="HG002")
    sample_track = next(t for t in vp_clcn["tracks"] if t["type"] == "sample_stacked")
    print("HG002 CLCNKB window bins:", sample_track["bin_count"], "variants:", sample_track["variant_count"])
    assert sample_track["variant_count"] >= 1

    bin_detail = store.sample_bin_variants("chr1", "HG002", 16_050_001, 16_060_000, limit=5)
    print("CLCNKB bin variants:", bin_detail["total_count"], "first:", bin_detail["variants"][0]["variant_id"])
    assert bin_detail["total_count"] >= 1
    assert bin_detail["variants"][0]["variant_id"] == "1-16051552-G-T"

    ref_pdf = generate_reference_report_pdf("1-16051552-G-T")
    print("reference PDF:", ref_pdf.name, ref_pdf.stat().st_size)
    assert ref_pdf.exists() and ref_pdf.stat().st_size > 1000


if __name__ == "__main__":
    main()
