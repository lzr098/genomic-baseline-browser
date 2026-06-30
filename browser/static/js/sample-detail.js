import { fetchSampleBinVariants } from "./api.js?v=12";

let selectedHit = null;
let context = null;
let tableState = {
  variants: [],
  sampleId: "",
  sortKey: "pos",
  sortDir: "asc",
};

const TABLE_COLUMNS = [
  { key: "variant_id", label: "Variant ID", className: "col-variant-id", sortable: true },
  { key: "source", label: "Source", className: "col-source", sortable: true },
  { key: "hgvs_consequence", label: "HGVS Consequence", className: "col-hgvs", sortable: true },
  { key: "vep_annotation", label: "VEP Annotation", className: "col-vep", sortable: true },
  { key: "lof_curation", label: "LoF Curation", className: "col-lof", sortable: true },
  { key: "germline_classification", label: "Germline classification", className: "col-germline", sortable: true },
  { key: "flags", label: "Flags", className: "col-flags", sortable: true },
  { key: "allele_count", label: "Allele Count", className: "col-num", sortable: true },
  { key: "allele_number", label: "Allele Number", className: "col-num", sortable: true },
  { key: "allele_frequency", label: "Allele Frequency", className: "col-num", sortable: true },
  { key: "homozygote_count", label: "Number of Homozygotes", className: "col-num", sortable: true },
  { key: "report", label: "report", className: "col-report", sortable: false },
];

function esc(text) {
  return String(text ?? "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

function formatAf(value) {
  if (value == null || Number.isNaN(Number(value))) return "—";
  const n = Number(value);
  if (n === 0) return "0";
  return n.toExponential(2);
}

function formatInt(value) {
  if (value == null || Number.isNaN(Number(value))) return "—";
  return Number(value).toLocaleString();
}

function reportPdfUrl(variantId, kind, sampleId) {
  const vid = encodeURIComponent(variantId);
  if (kind === "reference") {
    return `/api/variants/${vid}/report/reference.pdf`;
  }
  return `/api/variants/${vid}/report/sample.pdf?sample_id=${encodeURIComponent(sampleId)}`;
}

async function downloadReportPdf(button, variantId, kind, sampleId) {
  if (button.disabled) return;
  const label = button.textContent;
  button.disabled = true;
  button.textContent = "…";
  try {
    const res = await fetch(reportPdfUrl(variantId, kind, sampleId));
    if (!res.ok) {
      const err = await res.json().catch(() => ({}));
      throw new Error(err.detail || `PDF ${res.status}`);
    }
    const blob = await res.blob();
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${variantId}.${kind}.report.pdf`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  } catch (err) {
    window.alert(err.message || "Failed to generate PDF");
  } finally {
    button.disabled = false;
    button.textContent = label;
  }
}

function sortValue(variant, key) {
  switch (key) {
    case "variant_id":
      return variant.variant_id || "";
    case "source":
      return (variant.source_exome ? 2 : 0) + (variant.source_genome ? 1 : 0);
    case "flags":
      return (variant.flags || []).join(",");
    case "allele_count":
    case "allele_number":
    case "homozygote_count":
      return variant[key] == null ? Number.NEGATIVE_INFINITY : Number(variant[key]);
    case "allele_frequency":
      return variant.allele_frequency == null ? Number.NEGATIVE_INFINITY : Number(variant.allele_frequency);
    case "pos":
      return Number(variant.pos) || 0;
    default:
      return variant[key] ?? "";
  }
}

function compareSortValues(a, b) {
  const emptyA = a === "" || a === Number.NEGATIVE_INFINITY;
  const emptyB = b === "" || b === Number.NEGATIVE_INFINITY;
  if (emptyA && emptyB) return 0;
  if (emptyA) return 1;
  if (emptyB) return -1;
  if (typeof a === "number" && typeof b === "number") return a - b;
  return String(a).localeCompare(String(b), undefined, { numeric: true, sensitivity: "base" });
}

function sortedVariants(variants, sortKey, sortDir) {
  const mult = sortDir === "desc" ? -1 : 1;
  return [...variants].sort((left, right) => {
    let cmp = compareSortValues(sortValue(left, sortKey), sortValue(right, sortKey));
    if (cmp === 0 && sortKey !== "pos") {
      cmp = compareSortValues(sortValue(left, "pos"), sortValue(right, "pos"));
    }
    if (cmp === 0) {
      cmp = String(left.variant_id || "").localeCompare(String(right.variant_id || ""));
    }
    return mult * cmp;
  });
}

function renderSourceBadges(variant) {
  const badges = [];
  if (variant.source_exome) {
    badges.push('<span class="src-badge src-exome" title="gnomAD exomes">E</span>');
  }
  if (variant.source_genome) {
    badges.push('<span class="src-badge src-genome" title="gnomAD genomes">G</span>');
  }
  return badges.length ? badges.join("") : "—";
}

function renderVepCell(variant) {
  const label = variant.vep_annotation;
  if (!label) return "—";
  const category = variant.vep_category || "other";
  return `<span class="vep-ann vep-${esc(category)}"><span class="vep-dot" aria-hidden="true"></span>${esc(label)}</span>`;
}

function renderFlags(flags) {
  if (!flags?.length) return "—";
  return flags.map((flag) => `<span class="variant-flag">${esc(flag)}</span>`).join(" ");
}

function renderVariantIdCell(variant) {
  if (variant.variant_page) {
    return `<a class="variant-id-link" href="${esc(variant.variant_page)}" target="_blank" rel="noopener noreferrer">${esc(variant.variant_id)}</a>`;
  }
  return `<span class="variant-id-text">${esc(variant.variant_id)}</span>`;
}

function renderReportCell(variant, sampleId) {
  return `
    <div class="variant-report-actions">
      <button type="button" class="report-pdf-btn" data-report-kind="reference" data-variant-id="${esc(variant.variant_id)}" title="Baseline variant report (English PDF)">Baseline PDF</button>
      <button type="button" class="report-pdf-btn" data-report-kind="sample" data-variant-id="${esc(variant.variant_id)}" data-sample-id="${esc(sampleId)}" title="Sample vs baseline report (English PDF)">Sample PDF</button>
    </div>
  `;
}

function renderVariantRow(variant, sampleId) {
  const germline = variant.germline_classification
    ? `<span class="germline-link">${esc(variant.germline_classification)}</span>`
    : "—";

  return `
    <tr>
      <td class="col-variant-id">${renderVariantIdCell(variant)}</td>
      <td class="col-source">${renderSourceBadges(variant)}</td>
      <td class="col-hgvs">${variant.hgvs_consequence ? esc(variant.hgvs_consequence) : "—"}</td>
      <td class="col-vep">${renderVepCell(variant)}</td>
      <td class="col-lof">${variant.lof_curation ? esc(variant.lof_curation) : "—"}</td>
      <td class="col-germline">${germline}</td>
      <td class="col-flags">${renderFlags(variant.flags)}</td>
      <td class="col-num">${formatInt(variant.allele_count)}</td>
      <td class="col-num">${formatInt(variant.allele_number)}</td>
      <td class="col-num">${formatAf(variant.allele_frequency)}</td>
      <td class="col-num">${formatInt(variant.homozygote_count)}</td>
      <td class="col-report">${renderReportCell(variant, sampleId)}</td>
    </tr>
  `;
}

function renderTableHead(sortKey, sortDir) {
  return TABLE_COLUMNS.map((col) => {
    if (!col.sortable) {
      return `<th scope="col" class="${col.className}">${esc(col.label)}</th>`;
    }
    const active = col.key === sortKey;
    const ariaSort = active ? (sortDir === "asc" ? "ascending" : "descending") : "none";
    const indicator = active ? (sortDir === "asc" ? " ▲" : " ▼") : "";
    return (
      `<th scope="col" class="${col.className} col-sortable" aria-sort="${ariaSort}">` +
      `<button type="button" class="table-sort-btn" data-sort-key="${esc(col.key)}">${esc(col.label)}${indicator}</button>` +
      `</th>`
    );
  }).join("");
}

function renderVariantTable(variants, sampleId, sortKey, sortDir) {
  const rows = sortedVariants(variants, sortKey, sortDir);
  const body = rows.map((variant) => renderVariantRow(variant, sampleId)).join("");

  return `
    <div class="sample-variant-table-wrap">
      <table class="sample-variant-table">
        <thead>
          <tr>${renderTableHead(sortKey, sortDir)}</tr>
        </thead>
        <tbody>${body}</tbody>
      </table>
    </div>
  `;
}

function refreshVariantTableBody() {
  const body = document.getElementById("sample-bin-detail-body");
  if (!body || !tableState.variants.length) return;
  body.innerHTML = renderVariantTable(
    tableState.variants,
    tableState.sampleId,
    tableState.sortKey,
    tableState.sortDir
  );
}

function renderDetailPayload(payload) {
  const header = document.getElementById("sample-bin-detail-header");
  const body = document.getElementById("sample-bin-detail-body");
  const panel = document.getElementById("sample-bin-detail");
  if (!header || !body || !panel) return;

  const range = `${Number(payload.bin_start).toLocaleString()} – ${Number(payload.bin_end).toLocaleString()}`;
  const countNote = payload.truncated
    ? `Showing ${payload.returned_count} of ${payload.total_count} variants`
    : `${payload.total_count} variant${payload.total_count === 1 ? "" : "s"}`;

  header.innerHTML = `
    <strong>${esc(payload.sample_id)}</strong>
    <span class="sample-bin-detail-range">${esc(payload.chrom)}:${range}</span>
    <span class="sample-bin-detail-count">${esc(countNote)}</span>
  `;

  if (!payload.variants?.length) {
    tableState = { variants: [], sampleId: "", sortKey: "pos", sortDir: "asc" };
    body.innerHTML = `<p class="sample-bin-detail-empty">No variants in this bin.</p>`;
  } else {
    tableState = {
      variants: payload.variants,
      sampleId: payload.sample_id,
      sortKey: "pos",
      sortDir: "asc",
    };
    body.innerHTML = renderVariantTable(
      tableState.variants,
      tableState.sampleId,
      tableState.sortKey,
      tableState.sortDir
    );
  }

  panel.hidden = false;
}

function clearSelection() {
  if (selectedHit) {
    selectedHit.classList.remove("sample-bin-selected");
    selectedHit = null;
  }
}

export function clearSampleBinDetail() {
  clearSelection();
  tableState = { variants: [], sampleId: "", sortKey: "pos", sortDir: "asc" };
  const panel = document.getElementById("sample-bin-detail");
  const header = document.getElementById("sample-bin-detail-header");
  const body = document.getElementById("sample-bin-detail-body");
  if (panel) panel.hidden = true;
  if (header) header.innerHTML = "";
  if (body) body.innerHTML = "";
}

function markSelected(hit) {
  clearSelection();
  selectedHit = hit;
  hit.classList.add("sample-bin-selected");
}

async function onSampleBinClick(hit) {
  if (!context) return;
  const binStart = Number(hit.dataset.binStart);
  const binEnd = Number(hit.dataset.binEnd);
  if (!Number.isFinite(binStart) || !Number.isFinite(binEnd)) return;

  markSelected(hit);

  const header = document.getElementById("sample-bin-detail-header");
  const body = document.getElementById("sample-bin-detail-body");
  const panel = document.getElementById("sample-bin-detail");
  if (!header || !body || !panel) return;

  panel.hidden = false;
  header.innerHTML = `<strong>${esc(context.getSampleId())}</strong> <span class="sample-bin-detail-range">Loading bin ${binStart.toLocaleString()} – ${binEnd.toLocaleString()}…</span>`;
  body.innerHTML = `<p class="sample-bin-detail-empty">Loading variants…</p>`;

  try {
    const payload = await fetchSampleBinVariants(
      context.getChrom(),
      context.getSampleId(),
      binStart,
      binEnd
    );
    renderDetailPayload(payload);
    markSelected(hit);
  } catch (err) {
    header.innerHTML = `<strong>${esc(context.getSampleId())}</strong> <span class="sample-bin-detail-range">Bin load failed</span>`;
    body.innerHTML = `<p class="sample-bin-detail-empty">${esc(err.message)}</p>`;
  }
}

export function setupSampleBinDetail(ctx, root = document.getElementById("tracks-section")) {
  context = ctx;
  if (!root) return;

  if (root._sampleBinClick) {
    root.removeEventListener("click", root._sampleBinClick);
  }
  if (root._sampleReportClick) {
    root.removeEventListener("click", root._sampleReportClick);
  }
  if (root._tableSortClick) {
    root.removeEventListener("click", root._tableSortClick);
  }

  root._sampleBinClick = (ev) => {
    if (ev.target.closest?.(".report-pdf-btn, .table-sort-btn")) return;
    const hit = ev.target.closest?.(".sample-bin-hit");
    if (!hit || !root.contains(hit)) return;
    ev.preventDefault();
    onSampleBinClick(hit);
  };

  root._sampleReportClick = (ev) => {
    const btn = ev.target.closest?.(".report-pdf-btn");
    if (!btn || !root.contains(btn)) return;
    ev.preventDefault();
    ev.stopPropagation();
    const kind = btn.dataset.reportKind;
    const variantId = btn.dataset.variantId;
    const sampleId = btn.dataset.sampleId || context?.getSampleId?.();
    if (!variantId || !kind) return;
    if (kind === "sample" && !sampleId) return;
    downloadReportPdf(btn, variantId, kind, sampleId);
  };

  root._tableSortClick = (ev) => {
    const btn = ev.target.closest?.(".table-sort-btn");
    if (!btn || !root.contains(btn)) return;
    ev.preventDefault();
    ev.stopPropagation();
    const key = btn.dataset.sortKey;
    if (!key || !tableState.variants.length) return;
    if (tableState.sortKey === key) {
      tableState.sortDir = tableState.sortDir === "asc" ? "desc" : "asc";
    } else {
      tableState.sortKey = key;
      tableState.sortDir = key === "variant_id" || key === "germline_classification" ? "asc" : "desc";
    }
    refreshVariantTableBody();
  };

  root.addEventListener("click", root._sampleBinClick);
  root.addEventListener("click", root._sampleReportClick);
  root.addEventListener("click", root._tableSortClick);
}
