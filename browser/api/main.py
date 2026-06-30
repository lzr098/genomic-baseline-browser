"""FastAPI service for the baseline exome mutation landscape browser."""

from __future__ import annotations

import sys
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
BROWSER_ROOT = Path(__file__).resolve().parents[1]
STATIC_DIR = BROWSER_ROOT / "static"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from _browser_reports import (  # noqa: E402
    generate_reference_report_pdf,
    generate_sample_report_pdf,
)
from _browser_viewport import BrowserDataStore  # noqa: E402

app = FastAPI(
    title="Exome Baseline Genome Browser API",
    version="0.1.0",
    description="UCSC-style viewport API for gnomAD + ClinVar + GFF3 baseline tracks.",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET"],
    allow_headers=["*"],
)

store = BrowserDataStore()


@app.get("/api/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/api/manifest")
def get_manifest() -> dict:
    return store.api_manifest()


@app.get("/api/chrom/{chrom}/meta")
def get_chrom_meta(chrom: str) -> dict:
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    try:
        return store.chrom_info(chrom)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/api/chrom/{chrom}/viewport")
def get_viewport(
    chrom: str,
    start: int = Query(1, ge=1, description="1-based inclusive start"),
    end: int | None = Query(None, ge=1, description="1-based inclusive end"),
    sample: str | None = Query(None, description="Sample ID (e.g. HG002) for individual variant track"),
) -> dict:
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    try:
        meta = store.chrom_meta(chrom)
        end_val = int(end) if end is not None else int(meta["length"])
        sample_id = sample or store.manifest().get("default_sample")
        return store.build_viewport(chrom, start, end_val, sample_id=sample_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/api/chrom/{chrom}/sample/{sample_id}/bin-variants")
def get_sample_bin_variants(
    chrom: str,
    sample_id: str,
    bin_start: int = Query(..., ge=1, description="1-based inclusive bin start"),
    bin_end: int = Query(..., ge=1, description="1-based inclusive bin end"),
    limit: int = Query(100, ge=1, le=500, description="Max variants returned"),
) -> dict:
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    try:
        return store.sample_bin_variants(chrom, sample_id, bin_start, bin_end, limit=limit)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/api/variants/{variant_id}/report/reference.pdf")
def get_reference_variant_report_pdf(variant_id: str) -> FileResponse:
    try:
        pdf_path = generate_reference_report_pdf(variant_id, store.root)
        return FileResponse(
            pdf_path,
            media_type="application/pdf",
            filename=pdf_path.name,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc


@app.get("/api/variants/{variant_id}/report/sample.pdf")
def get_sample_variant_report_pdf(
    variant_id: str,
    sample_id: str = Query(..., description="Sample ID for compare report"),
) -> FileResponse:
    try:
        pdf_path = generate_sample_report_pdf(variant_id, sample_id, store.root)
        return FileResponse(
            pdf_path,
            media_type="application/pdf",
            filename=pdf_path.name,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc


if STATIC_DIR.is_dir():
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


@app.get("/")
def index() -> FileResponse:
    index_path = STATIC_DIR / "index.html"
    if not index_path.exists():
        raise HTTPException(status_code=404, detail="frontend not built")
    return FileResponse(index_path)
