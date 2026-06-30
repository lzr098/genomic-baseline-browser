#!/usr/bin/env python3
"""Prepare static data for the interactive baseline genome browser (Step A0)."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from _bin_builder import chrom_sort_key  # noqa: E402
from _browser_gff3 import extract_genes_from_gff3  # noqa: E402
from _config import (  # noqa: E402
    ALL_CHROMSOMES,
    BROWSER,
    BROWSER_GENES,
    CLICK_JUMP_BP,
    DEFAULT_BROWSER_CHROM,
    DEFAULT_BROWSER_SAMPLE,
    GFF3,
    ROOT,
    ZOOM_STEPS_BP,
)
from _browser_sample import list_browser_samples  # noqa: E402
from _cytobands import load_cytobands  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--chrom",
        default="all",
        help="One chromosome (chr21) or 'all' (default: all)",
    )
    parser.add_argument(
        "--skip-gff3",
        action="store_true",
        help="Skip GFF3 gene extraction (cytobands + manifest only)",
    )
    return parser.parse_args()


def load_build_manifest() -> dict:
    """Load the baseline build manifest; falls back to browser manifest if needed."""
    manifest_path = ROOT / "build_manifest.json"
    if not manifest_path.exists():
        raise SystemExit(f"missing {manifest_path}")
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def write_cytobands_json(out_path: Path) -> dict[str, list[dict[str, object]]]:
    all_bands: dict[str, list[dict[str, object]]] = {}
    for chrom in ALL_CHROMSOMES:
        all_bands[chrom] = load_cytobands(chrom)
    out_path.write_text(json.dumps(all_bands, ensure_ascii=False, indent=2), encoding="utf-8")
    return all_bands


def write_gene_parquets(
    genes_by_chrom: dict[str, list[dict[str, object]]],
    chroms: tuple[str, ...],
) -> dict[str, str]:
    BROWSER_GENES.mkdir(parents=True, exist_ok=True)
    rel_paths: dict[str, str] = {}
    for chrom in chroms:
        rows = genes_by_chrom.get(chrom, [])
        out = BROWSER_GENES / f"{chrom}_genes.parquet"
        frame = pd.DataFrame(rows)
        if frame.empty:
            frame = pd.DataFrame(
                columns=[
                    "gene_id",
                    "gene_name",
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "biotype",
                    "transcript_id",
                    "transcript_name",
                    "exon_count",
                    "exons_json",
                ]
            )
        frame.to_parquet(out, index=False)
        rel_paths[chrom] = str(out.relative_to(ROOT))
        print(f"  genes {chrom}: {len(frame):,} → {out}")
    return rel_paths


def write_browser_manifest(
    build_manifest: dict,
    gene_paths: dict[str, str],
    out_path: Path,
) -> dict:
    chromosomes: list[dict] = []
    for item in sorted(build_manifest["chromosomes"], key=lambda x: chrom_sort_key(x["chrom"])):
        chrom = item["chrom"]
        chromosomes.append(
            {
                "chrom": chrom,
                "length": int(item["chrom_length"]),
                "variant_count": int(item["variant_allele_count"]),
                "bins": item["output_paths"],
                "genes_parquet": gene_paths.get(chrom, f"processed/browser/genes/{chrom}_genes.parquet"),
            }
        )

    browser_manifest = {
        "assembly": build_manifest.get("reference_genome", "GRCh38"),
        "gnomad": build_manifest.get("gnomad", "v4.1"),
        "gnomad_callset": build_manifest.get("gnomad_callset", "exomes"),
        "clinvar": build_manifest.get("clinvar", "GRCh38 VCF"),
        "gff3": str(GFF3),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "default_chrom": DEFAULT_BROWSER_CHROM,
        "click_jump_bp": CLICK_JUMP_BP,
        "zoom_steps_bp": list(ZOOM_STEPS_BP),
        "track_order": [
            "genes",
            "gnomad",
            "clinvar",
            "sample",
        ],
        "tracks": {
            "genes": {
                "label": "Genes",
                "sublabel": "Ensembl 115 GFF3 · exons",
                "type": "gene",
            },
            "gnomad": {
                "label": "gnomAD_exome",
                "sublabel": "v4.1 exomes",
                "type": "histogram",
                "value_col": "variant_count",
                "fill": "#0C5EC4",
            },
            "clinvar": {
                "label": "ClinVar",
                "sublabel": "pathogenic · VUS/conflicting · benign",
                "type": "clinvar_stacked",
            },
            "sample": {
                "label": "{sample_id}",
                "sublabel": "known · novel stacked",
                "type": "sample_stacked",
            },
        },
        "default_sample": DEFAULT_BROWSER_SAMPLE,
        "samples": list_browser_samples(ROOT),
        "chromosomes": chromosomes,
        "paths": {
            "cytobands": "processed/browser/cytobands.json",
            "genes_dir": "processed/browser/genes",
        },
    }
    out_path.write_text(json.dumps(browser_manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return browser_manifest


def main() -> None:
    args = parse_args()
    chroms = ALL_CHROMSOMES if args.chrom == "all" else (args.chrom,)
    if args.chrom != "all" and args.chrom not in ALL_CHROMSOMES:
        raise SystemExit(f"unknown chrom: {args.chrom}")

    BROWSER.mkdir(parents=True, exist_ok=True)
    build_manifest = load_build_manifest()

    print("writing cytobands.json …")
    cytobands_path = BROWSER / "cytobands.json"
    all_bands = write_cytobands_json(cytobands_path)
    print(f"  {cytobands_path} ({sum(len(v) for v in all_bands.values())} bands)")

    gene_paths: dict[str, str] = {}
    if not args.skip_gff3:
        if not GFF3.exists():
            raise SystemExit(f"missing GFF3: {GFF3}")
        print(f"extracting genes from {GFF3} …")
        genes_by_chrom = extract_genes_from_gff3(GFF3, ALL_CHROMSOMES)
        gene_paths = write_gene_parquets(genes_by_chrom, chroms)
    else:
        for chrom in chroms:
            gene_paths[chrom] = f"processed/browser/genes/{chrom}_genes.parquet"

    manifest_path = BROWSER / "manifest.json"
    browser_manifest = write_browser_manifest(build_manifest, gene_paths, manifest_path)
    print(f"  {manifest_path}")
    print(
        f"done. chromosomes={len(browser_manifest['chromosomes'])}, "
        f"default={browser_manifest['default_chrom']}, "
        f"track_order={browser_manifest['track_order']}"
    )


if __name__ == "__main__":
    main()
