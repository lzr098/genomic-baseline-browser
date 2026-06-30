#!/usr/bin/env python3
"""Compare an individual sample VCF against gnomAD exome + ClinVar baseline.

If gnomAD VCFs are not configured (GNOMAD_EXOMES_VCF_DIR env var), sample
comparison is skipped with a clear message and no compare parquet is written.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from _compare import compare_sample
from _config import ALL_CHROMSOMES, COMPARE, GNOMAD_EXOMES_VCF_DIR, ROOT


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-vcf", type=Path, required=True, help="Individual VCF (.vcf or .vcf.gz)")
    parser.add_argument("--sample-id", required=True, help="Sample ID / VCF column name")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: compare/{sample_id})",
    )
    parser.add_argument(
        "--chrom",
        default="all",
        help="Limit to one chromosome or 'all' (default: all)",
    )
    return parser.parse_args()


def _gnomad_available_for(chromosomes: tuple[str, ...]) -> bool:
    """Return True when the configured gnomAD exome VCFs are actually present."""
    if not GNOMAD_EXOMES_VCF_DIR.exists():
        return False
    return all(gnomad_vcf.exists() for gnomad_vcf in (GNOMAD_EXOMES_VCF_DIR / f"gnomad.exomes.v4.1.sites.{c}.vcf.bgz" for c in chromosomes))


def main() -> None:
    args = parse_args()
    if not args.sample_vcf.exists():
        raise SystemExit(f"missing sample VCF: {args.sample_vcf}")

    build_manifest = ROOT / "build_manifest.json"
    if not build_manifest.exists():
        raise SystemExit(f"missing baseline manifest: {build_manifest}")

    output_dir = args.output_dir or (COMPARE / args.sample_id)
    chromosomes = ALL_CHROMSOMES if args.chrom == "all" else (args.chrom,)

    if not _gnomad_available_for(chromosomes):
        print("gnomAD exome VCFs not available; skipping sample comparison.")
        print(f"  set GNOMAD_EXOMES_VCF_DIR to enable comparison.")
        summary_path = output_dir / "compare_summary.json"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(
            json.dumps(
                {
                    "sample_id": args.sample_id,
                    "sample_vcf": str(args.sample_vcf.resolve()),
                    "skipped": True,
                    "reason": "gnomAD VCFs not configured",
                    "counts": {
                        "total_variants": 0,
                        "known_variants": 0,
                        "novel_variants": 0,
                        "by_match_status": {},
                    },
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"  wrote: {summary_path}")
        return

    result = compare_sample(
        sample_vcf=args.sample_vcf.resolve(),
        sample_id=args.sample_id,
        output_dir=output_dir.resolve(),
        chromosomes=chromosomes,
    )
    summary = result["summary"]
    print(f"sample {args.sample_id}: {summary['counts']['total_variants']:,} variants")
    print(f"  known: {summary['counts']['known_variants']:,}")
    print(f"  novel: {summary['counts']['novel_variants']:,}")
    print(f"  by_match_status: {summary['counts']['by_match_status']}")
    for name, path in result["paths"].items():
        size = Path(path).stat().st_size
        print(f"  {name}: {path} ({size:,} bytes)")
    print("done.")


if __name__ == "__main__":
    main()
