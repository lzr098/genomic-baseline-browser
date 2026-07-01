import { fetchConfig, fetchSampleBinVariants, fetchSampleViewportVariants, fetchVariantGnomad } from "./api.js?v=15";

let selectedHit = null;
let context = null;
let config = null;
let tableState = {
  variants: [],
  sampleId: "",
  sortKey: "pos",
  sortDir: "asc",
  mode: "viewport", // "viewport" or "bin"
  rangeLabel: "",
};

const TABLE_COLUMNS = [
  { key: "locus", label: "Locus", className: "col-locus", sortable: true },
  { key: "gene_name", label: "Gene", className: "col-gene", sortable: true },
  { key: "location", label: "Location", className: "col-location", sortable: true },
  { key: "consequence", label: "Consequence", className: "col-consequence", sortable: true },
  { key: "variant", label: "Ref → Alt", className: "col-variant", sortable: false },
  { key: "vep_annotation", label: "VEP", className: "col-vep", sortable: true },
  { key: "germline_classification", label: "ClinVar", className: "col-clinvar", sortable: true },
  { key: "allele_frequency", label: "gnomAD AF", className: "col-af", sortable: true },
  { key: "flags", label: "Flags", className: "col-flags", sortable: false },
  { key: "report", label: "Actions", className: "col-report", sortable: false },
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
  if (n < 0.001) return n.toExponential(2);
  return n.toFixed(4);
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

function geneName(variant) {
  const ctx = variant.gene_context;
  if (ctx && ctx.gene_name) return ctx.gene_name;
  return "—";
}

function geneLocation(variant) {
  const ctx = variant.gene_context;
  if (!ctx) return "intergenic";
  const loc = ctx.location || "";
  if (loc === "intron" && ctx.distance_to_exon != null) {
    return `intron (${ctx.distance_to_exon.toLocaleString()} bp to exon)`;
  }
  return loc;
}

function consequence(variant) {
  if (variant.vep_annotation) return variant.vep_annotation;
  if (variant.location_consequence) return variant.location_consequence;
  return "—";
}

function vepCategoryClass(variant) {
  const cat = variant.vep_category || "other";
  return `vep-${esc(cat)}`;
}

function sortValue(variant, key) {
  switch (key) {
    case "locus":
      return Number(variant.pos) || 0;
    case "gene_name":
      return geneName(variant);
    case "location":
      return geneLocation(variant);
    case "consequence":
      return consequence(variant);
    case "variant":
      return `${variant.ref || ""}→${variant.alt || ""}`;
    case "vep_annotation":
      return variant.vep_annotation || "";
    case "germline_classification":
      return variant.germline_classification || "";
    case "allele_frequency":
      return variant.allele_frequency == null ? Number.NEGATIVE_INFINITY : Number(variant.allele_frequency);
    case "flags":
      return (variant.flags || []).join(",");
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

function renderVariantCell(variant) {
  const ref = variant.ref || "—";
  const alt = variant.alt || "—";
  return `${esc(ref)} <span style="color:#999">→</span> ${esc(alt)}`;
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

function renderAfCell(variant, index) {
  const value = formatAf(variant.allele_frequency);
  const page = variant.variant_page;
  if (!page) return value;
  if (variant._querying) {
    return `<span class="gnomad-query-link" aria-busy="true">…</span>`;
  }
  if (value === "—") {
    return `<button type="button" class="gnomad-query-link query-af-btn" data-variant-idx="${index}" title="Query gnomAD for this variant">Query ↗</button>`;
  }
  return `<a class="gnomad-af-link" href="${esc(page)}" target="_blank" rel="noopener noreferrer" title="Open in gnomAD">${value} ↗</a>`;
}

function renderReportCell(variant, sampleId) {
  const gnomadConfigured = config?.gnomad_exomes_vcf_dir_configured ?? false;
  const disabled = gnomadConfigured ? "" : " disabled";
  const title = gnomadConfigured
    ? "Reference / baseline variant report (gnomAD-based)"
    : "Reports require a local gnomAD exome VCF. Set GNOMAD_EXOMES_VCF_DIR and restart.";
  return `
    <div class="variant-report-actions">
      <button type="button" class="report-pdf-btn"${disabled} data-report-kind="reference" data-variant-id="${esc(variant.variant_id)}" title="${title}">Baseline</button>
      <button type="button" class="report-pdf-btn"${disabled} data-report-kind="sample" data-variant-id="${esc(variant.variant_id)}" data-sample-id="${esc(sampleId)}" title="Sample-specific variant report">Sample</button>
    </div>
  `;
}

function renderVariantRow(variant, sampleId, index) {
  const clin = variant.germline_classification
    ? `<span class="germline-link">${esc(variant.germline_classification)}</span>`
    : "—";

  return `
    <tr>
      <td class="col-locus"><a class="variant-id-link" href="${esc(variant.variant_page || "#")}" target="_blank" rel="noopener noreferrer">${esc(variant.locus || `${variant.chrom}:${variant.pos}`)}</a></td>
      <td class="col-gene">${esc(geneName(variant))}</td>
      <td class="col-location">${esc(geneLocation(variant))}</td>
      <td class="col-consequence">${esc(consequence(variant))}</td>
      <td class="col-variant">${renderVariantCell(variant)}</td>
      <td class="col-vep">${renderVepCell(variant)}</td>
      <td class="col-clinvar">${clin}</td>
      <td class="col-af">${renderAfCell(variant, index)}</td>
      <td class="col-flags">${renderFlags(variant.flags)}</td>
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
  const body = rows.map((variant, index) => renderVariantRow(variant, sampleId, index)).join("");

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

function renderDetailHeader() {
  const header = document.getElementById("sample-bin-detail-header");
  if (!header) return;

  const { mode, rangeLabel, variants, sampleId, truncated, totalCount } = tableState;
  const countNote = totalCount === 0
    ? "0 variants"
    : truncated
      ? `Showing ${variants.length} of ${totalCount} variants`
      : `${totalCount} variant${totalCount === 1 ? "" : "s"}`;
  const modeLabel = mode === "bin" ? "bin" : "current view";

  header.innerHTML = `
    <strong>${esc(sampleId)}</strong>
    <span class="sample-bin-detail-range">${esc(rangeLabel)}</span>
    <span class="sample-bin-detail-count">${esc(countNote)} in ${esc(modeLabel)}</span>
    ${mode === "bin" ? `<button type="button" id="show-viewport-variants" class="report-pdf-btn">Show all in view</button>` : ""}
  `;

  const showAllBtn = document.getElementById("show-viewport-variants");
  if (showAllBtn) {
    showAllBtn.addEventListener("click", () => loadViewportVariants());
  }
}

function renderDetailBody() {
  const body = document.getElementById("sample-bin-detail-body");
  if (!body) return;
  if (!tableState.variants.length) {
    body.innerHTML = `<p class="sample-bin-detail-empty">No variants in this region.</p>`;
  } else {
    body.innerHTML = renderVariantTable(
      tableState.variants,
      tableState.sampleId,
      tableState.sortKey,
      tableState.sortDir
    );
  }
}

function renderDetailPanel() {
  const panel = document.getElementById("sample-bin-detail");
  if (!panel) return;
  panel.hidden = false;
  renderDetailHeader();
  renderDetailBody();
}

function setTableState(updates) {
  tableState = { ...tableState, ...updates };
  renderDetailPanel();
}

function clearSelection() {
  if (selectedHit) {
    selectedHit.classList.remove("sample-bin-selected");
    selectedHit = null;
  }
}

function markSelected(hit) {
  clearSelection();
  selectedHit = hit;
  hit.classList.add("sample-bin-selected");
}

export function clearSampleBinDetail() {
  clearSelection();
  tableState = { variants: [], sampleId: "", sortKey: "pos", sortDir: "asc", mode: "viewport", rangeLabel: "" };
  const panel = document.getElementById("sample-bin-detail");
  const header = document.getElementById("sample-bin-detail-header");
  const body = document.getElementById("sample-bin-detail-body");
  if (panel) panel.hidden = true;
  if (header) header.innerHTML = "";
  if (body) body.innerHTML = "";
}

async function loadViewportVariants() {
  if (!context) return;
  const sampleId = context.getSampleId();
  const chrom = context.getChrom();
  const { start, end } = context.getViewport() || { start: 1, end: 1 };
  if (!sampleId || !chrom || end < start) return;

  setTableState({
    mode: "viewport",
    rangeLabel: `${chrom}:${start.toLocaleString()} – ${end.toLocaleString()}`,
    sampleId,
    variants: [],
    totalCount: 0,
    truncated: false,
  });

  try {
    const payload = await fetchSampleViewportVariants(chrom, sampleId, start, end, 500);
    setTableState({
      variants: payload.variants || [],
      totalCount: payload.total_count || 0,
      truncated: payload.truncated || false,
      sampleId: payload.sample_id || sampleId,
    });
  } catch (err) {
    setTableState({ variants: [], totalCount: 0, truncated: false });
    const body = document.getElementById("sample-bin-detail-body");
    if (body) body.innerHTML = `<p class="sample-bin-detail-empty">${esc(err.message)}</p>`;
  }
}

async function onSampleBinClick(hit) {
  if (!context) return;
  const binStart = Number(hit.dataset.binStart);
  const binEnd = Number(hit.dataset.binEnd);
  if (!Number.isFinite(binStart) || !Number.isFinite(binEnd)) return;

  markSelected(hit);

  const sampleId = context.getSampleId();
  const chrom = context.getChrom();
  setTableState({
    mode: "bin",
    rangeLabel: `${chrom}:${binStart.toLocaleString()} – ${binEnd.toLocaleString()}`,
    sampleId,
    variants: [],
    totalCount: 0,
    truncated: false,
  });

  try {
    const payload = await fetchSampleBinVariants(chrom, sampleId, binStart, binEnd, 100);
    setTableState({
      variants: payload.variants || [],
      totalCount: payload.total_count || 0,
      truncated: payload.truncated || false,
      sampleId: payload.sample_id || sampleId,
    });
    markSelected(hit);
  } catch (err) {
    setTableState({ variants: [], totalCount: 0, truncated: false });
    const body = document.getElementById("sample-bin-detail-body");
    if (body) body.innerHTML = `<p class="sample-bin-detail-empty">${esc(err.message)}</p>`;
  }
}

async function onQueryAfClick(btn) {
  const idx = Number(btn.dataset.variantIdx);
  if (!Number.isFinite(idx) || idx < 0 || idx >= tableState.variants.length) return;
  const variant = tableState.variants[idx];
  if (!variant || variant.allele_frequency != null || variant._querying) return;

  variant._querying = true;
  refreshVariantTableBody();

  try {
    const data = await fetchVariantGnomad(variant.variant_id);
    variant.allele_frequency = data.af ?? null;
    variant.allele_count = data.exome?.ac ?? data.genome?.ac ?? null;
    variant.allele_number = data.exome?.an ?? data.genome?.an ?? null;
    variant.homozygote_count = null;
    variant._querying = false;
  } catch (err) {
    variant._querying = false;
    variant._queryError = err.message || "Query failed";
  }
  refreshVariantTableBody();
}

export function setupSampleBinDetail(ctx, root = document.getElementById("tracks-section")) {
  context = ctx;
  if (!root) return;

  // Load config once; it is used to decide whether report buttons are enabled.
  fetchConfig().then((cfg) => {
    config = cfg;
    refreshVariantTableBody();
  }).catch(() => {
    config = { gnomad_exomes_vcf_dir_configured: false };
    refreshVariantTableBody();
  });

  if (root._sampleBinClick) root.removeEventListener("click", root._sampleBinClick);
  if (root._sampleReportClick) root.removeEventListener("click", root._sampleReportClick);
  if (root._tableSortClick) root.removeEventListener("click", root._tableSortClick);
  if (root._queryAfClick) root.removeEventListener("click", root._queryAfClick);

  root._sampleBinClick = (ev) => {
    if (ev.target.closest?.(".report-pdf-btn, .table-sort-btn, #show-viewport-variants, .query-af-btn")) return;
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
    if (btn.disabled) {
      window.alert("Reports require a local gnomAD exome VCF. Set GNOMAD_EXOMES_VCF_DIR and restart.");
      return;
    }
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
      tableState.sortDir = key === "locus" || key === "gene_name" || key === "location" ? "asc" : "desc";
    }
    refreshVariantTableBody();
    renderDetailHeader();
  };

  root._queryAfClick = (ev) => {
    const btn = ev.target.closest?.(".query-af-btn");
    if (!btn || !root.contains(btn)) return;
    ev.preventDefault();
    ev.stopPropagation();
    onQueryAfClick(btn);
  };

  root.addEventListener("click", root._sampleBinClick);
  root.addEventListener("click", root._sampleReportClick);
  root.addEventListener("click", root._tableSortClick);
  root.addEventListener("click", root._queryAfClick);

  // Load the viewport-wide variant list by default
  loadViewportVariants();
}

export function refreshViewportVariants() {
  loadViewportVariants();
}
