"""On-demand variant report PDF generation for the interactive browser."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from _browser_sample import _load_sample_frame
from _config import COMPARE, ROOT, SAMPLES
from _variant_report import (
    build_reference_variant_report,
    build_variant_report,
    parse_gnomad_variant_id,
    render_variant_report_pdf,
    report_stem,
    write_report_json,
)


def _report_cache_dir(root: Path) -> Path:
    path = root / "processed" / "browser" / "report_cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_sample_vcf(sample_id: str, root: Path = ROOT) -> Path:
    compare_dir = root / COMPARE / sample_id
    summary_path = compare_dir / "compare_summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        sample_vcf = Path(summary["sample_vcf"])
        if sample_vcf.exists():
            return sample_vcf.resolve()
    default = root / SAMPLES / f"{sample_id}.vcf.gz"
    if default.exists():
        return default.resolve()
    raise FileNotFoundError(f"missing sample VCF for {sample_id}")


def load_compare_row(
    sample_id: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    root: Path = ROOT,
) -> dict[str, Any] | None:
    frame = _load_sample_frame(sample_id, str(root))
    sub = frame[
        (frame["chrom"] == chrom)
        & (frame["pos"] == int(pos))
        & (frame["ref"] == ref)
        & (frame["alt"] == alt)
    ]
    if sub.empty:
        return None
    return sub.iloc[0].to_dict()


def reference_report_pdf_path(variant_id: str, root: Path = ROOT) -> Path:
    cache = _report_cache_dir(root) / "reference"
    cache.mkdir(parents=True, exist_ok=True)
    stem = variant_id.strip()
    return cache / f"{stem}.report.pdf"


def sample_report_pdf_path(variant_id: str, sample_id: str, root: Path = ROOT) -> Path:
    cache = _report_cache_dir(root) / "sample" / sample_id
    cache.mkdir(parents=True, exist_ok=True)
    stem = variant_id.strip()
    return cache / f"{stem}.report.pdf"


def generate_reference_report_pdf(variant_id: str, root: Path = ROOT) -> Path:
    """Build English reference (baseline) variant report PDF."""
    pdf_path = reference_report_pdf_path(variant_id, root)
    report = build_reference_variant_report(variant_id)
    json_path = pdf_path.with_suffix(".report.json")
    digest = write_report_json(report, json_path)
    render_variant_report_pdf(report, pdf_path, json_hash=digest, locale="en")
    return pdf_path


def generate_sample_report_pdf(
    variant_id: str,
    sample_id: str,
    root: Path = ROOT,
) -> Path:
    """Build English sample-vs-baseline variant report PDF."""
    chrom, pos, ref, alt = parse_gnomad_variant_id(variant_id)
    sample_vcf = resolve_sample_vcf(sample_id, root)
    compare_dir = root / COMPARE / sample_id
    compare_row = load_compare_row(sample_id, chrom, pos, ref, alt, root)

    report = build_variant_report(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        sample_id=sample_id,
        sample_vcf=sample_vcf,
        compare_row=compare_row,
        compare_dir=compare_dir,
    )
    pdf_path = sample_report_pdf_path(variant_id, sample_id, root)
    stem = report_stem(chrom, pos, ref, alt)
    json_path = pdf_path.parent / f"{stem}.report.json"
    digest = write_report_json(report, json_path)
    render_variant_report_pdf(report, pdf_path, json_hash=digest, locale="en")
    return pdf_path
