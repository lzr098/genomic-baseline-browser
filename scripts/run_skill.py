#!/usr/bin/env python3
"""One-click launcher for the exome baseline browser.

Usage:
    python scripts/run_skill.py /path/to/sample.vcf.gz [sample_id]

Steps:
1. Detect sample name and GRCh38 reference.
2. Compare sample against gnomAD baseline if gnomAD VCFs are configured.
3. Regenerate browser manifest.
4. Start FastAPI server and print the URL.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path

from pysam import VariantFile

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

from _config import GNOMAD_EXOMES_VCF_DIR  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_vcf", type=Path, help="GRCh38 sample VCF (.vcf or .vcf.gz)")
    parser.add_argument(
        "--sample-id",
        type=str,
        default=None,
        help="Sample ID (default: inferred from VCF filename or SAMPLE header)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=int(os.environ.get("PORT", "8765")),
        help="Browser API port (default: 8765)",
    )
    parser.add_argument(
        "--skip-compare",
        action="store_true",
        help="Skip gnomAD comparison even if VCFs are available",
    )
    return parser.parse_args()


def infer_sample_id(vcf_path: Path) -> str:
    stem = vcf_path.name.removesuffix(".gz").removesuffix(".vcf")
    if stem:
        return stem
    return "sample"


def sample_id_from_vcf(vcf_path: Path) -> str | None:
    try:
        with VariantFile(str(vcf_path)) as vcf:
            samples = list(vcf.header.samples)
            if samples:
                return str(samples[0])
    except Exception:
        return None
    return None


def detect_reference(vcf_path: Path) -> str:
    """Best-effort GRCh38 detection from VCF reference header."""
    ref = "unknown"
    try:
        with VariantFile(str(vcf_path)) as vcf:
            for record in vcf.header.records:
                if record.type == "GENERIC" and str(record.key) == "reference":
                    ref = str(record.value)
                    break
    except Exception:
        pass

    if not re.search(r"GRCh38|hg38|Homo_sapiens_assembly38", ref, re.IGNORECASE):
        print(f"WARNING: reference does not look like GRCh38: {ref}")
    return ref


def gnomad_available() -> bool:
    if not GNOMAD_EXOMES_VCF_DIR.exists():
        return False
    for chrom in (f"chr{i}" for i in range(1, 23)):
        vcf = GNOMAD_EXOMES_VCF_DIR / f"gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz"
        if not vcf.exists():
            return False
    return True


def run_compare(sample_vcf: Path, sample_id: str) -> None:
    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "02_compare_sample.py"),
        "--sample-vcf",
        str(sample_vcf),
        "--sample-id",
        sample_id,
        "--chrom",
        "all",
    ]
    print("$ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_prepare_browser() -> None:
    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "06_prepare_browser_data.py"),
        "--chrom",
        "all",
        "--skip-gff3",
    ]
    print("$ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def start_server(port: int) -> subprocess.Popen:
    env = os.environ.copy()
    env["PORT"] = str(port)
    cmd = [
        sys.executable,
        "-m",
        "uvicorn",
        "browser.api.main:app",
        "--host",
        "127.0.0.1",
        "--port",
        str(port),
    ]
    print("$ " + " ".join(cmd))
    return subprocess.Popen(
        cmd,
        cwd=ROOT,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


def wait_for_server(port: int, timeout: int = 30) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            import urllib.request
            with urllib.request.urlopen(f"http://127.0.0.1:{port}/api/health", timeout=2) as resp:
                if resp.status == 200:
                    return True
        except Exception:
            pass
        time.sleep(0.5)
    return False


def main() -> None:
    args = parse_args()
    sample_vcf = args.sample_vcf
    if not sample_vcf.exists():
        raise SystemExit(f"missing VCF: {sample_vcf}")

    sample_id = args.sample_id or sample_id_from_vcf(sample_vcf) or infer_sample_id(sample_vcf)
    print(f"sample VCF: {sample_vcf}")
    print(f"sample ID: {sample_id}")

    ref = detect_reference(sample_vcf)
    print(f"reference: {ref}")

    if not args.skip_compare and gnomad_available():
        run_compare(sample_vcf, sample_id)
    else:
        print("skipping gnomAD comparison (not configured or --skip-compare)")

    run_prepare_browser()

    port = args.port
    proc = start_server(port)
    print(f"\nStarting server on http://127.0.0.1:{port}/")
    if wait_for_server(port):
        print(f"Browser ready: http://127.0.0.1:{port}/")
        print(f"Manifest: http://127.0.0.1:{port}/api/manifest")
    else:
        print("ERROR: server did not start within timeout")
        proc.terminate()
        sys.exit(1)

    print("\nPress Ctrl+C to stop the server.")
    try:
        while True:
            line = proc.stdout.readline()
            if not line:
                break
            print(line, end="")
    except KeyboardInterrupt:
        print("\nStopping server...")
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()


if __name__ == "__main__":
    main()
