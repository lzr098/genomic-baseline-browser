import { LABEL_WIDTH, TRACK_PAD_X, xp } from "./coords.js?v=8";
import {
  IDEOGRAM_HEIGHT,
  ideogramLabel,
  renderIdeogramSvg,
  updateIdeogramLabelText,
  updateIdeogramMarker,
} from "./ideogram.js?v=8";
import { setupBinTooltips, tipAttr } from "./bin-tooltips.js?v=8";

export { setupBinTooltips };

const UCSC_RULER = "#666666";
const EXON_FILL = "#697e15";
const EXON_STROKE = "#4f5f10";
const BROWSER_FONT = "Helvetica, Arial, sans-serif";
const TRACK_SUBGRID_STROKE = "#C8DDF5";
const CLINVAR_PLP = "#C62828";
const CLINVAR_VUS = "#F9A825";
const CLINVAR_BENIGN = "#7CB342";

const SAMPLE_KNOWN = "#B8B8B8";
const SAMPLE_NOVEL = "#D84315";

const HEIGHT = {
  ideogram: IDEOGRAM_HEIGHT,
  ruler: 42,
  histogram: 66,
  clinvar: 132,
  gene: 84,
  sample: 132,
};

function esc(text) {
  return String(text ?? "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

function renderLabel(name, hint, { swatchClass = "" } = {}) {
  const swatchExtra = name ? swatchClass : " track-swatch-ruler";
  const swatch = `<span class="track-swatch${swatchExtra}"></span>`;
  const body = name
    ? `<div class="track-label">
        <span class="track-label-name">${esc(name)}</span>
        <div class="track-label-tooltip" role="tooltip">${hint}</div>
      </div>`
    : "";
  return `<div class="track-label-row">${swatch}${body}</div>`;
}

function svgOpen(trackWidth, height) {
  return `<svg viewBox="0 0 ${trackWidth} ${height}" preserveAspectRatio="none" font-family="${BROWSER_FONT}" xmlns="http://www.w3.org/2000/svg">`;
}

function subgridStepBp(ruler) {
  const main = Number(ruler?.step_bp) || 10_000;
  return Math.max(1, Math.floor(main / 10));
}

/** UCSC-style faint vertical guides at one-tenth of the main ruler spacing. */
function renderSubgridLines(ruler, start, end, trackWidth, height) {
  const step = subgridStepBp(ruler);
  const parts = [];
  let pos = Math.ceil(start / step) * step;
  for (; pos <= end; pos += step) {
    const x = xp(pos, start, end, trackWidth);
    parts.push(
      `<line class="track-subgrid-line" x1="${x.toFixed(2)}" y1="0" x2="${x.toFixed(2)}" y2="${height}" ` +
        `stroke="${TRACK_SUBGRID_STROKE}" stroke-width="0.35"/>`
    );
  }
  return parts.join("");
}

export function renderRulerSvg(ruler, start, end, trackWidth, height = HEIGHT.ruler) {
  const lineY = 14;
  const pad = TRACK_PAD_X;
  const parts = [];
  parts.push(svgOpen(trackWidth, height));
  parts.push(
    `<line x1="${pad}" y1="${lineY}" x2="${trackWidth - pad}" y2="${lineY}" stroke="${UCSC_RULER}" stroke-width="0.8"/>`
  );
  for (const tick of ruler?.ticks || []) {
    const x = xp(tick.pos, start, end, trackWidth);
    parts.push(
      `<line x1="${x.toFixed(2)}" y1="${lineY - 5}" x2="${x.toFixed(2)}" y2="${lineY + 5}" stroke="${UCSC_RULER}" stroke-width="0.8"/>`
    );
    if (tick.label) {
      parts.push(
        `<text x="${x.toFixed(2)}" y="${lineY + 22}" font-size="11" fill="${UCSC_RULER}" text-anchor="middle">${esc(tick.label)}</text>`
      );
    }
  }
  parts.push("</svg>");
  return parts.join("");
}

export function renderHistogramSvg(track, start, end, trackWidth, height = HEIGHT.histogram, ruler = null) {
  const bins = track.bins || [];
  const baseline = height - 8;
  const top = 10;
  const hAvail = baseline - top;
  const PEAK_FRAC = 0.9;
  const vtotal = Math.max(1, bins.reduce((sum, b) => sum + (b.value || 0), 0));
  const vmax = Math.max(1, ...bins.map((b) => b.value || 0));
  // Share of total counts, then scale so the tallest bar uses PEAK_FRAC of track height.
  const peakScale = (PEAK_FRAC * vtotal) / vmax;
  const parts = [];
  parts.push(svgOpen(trackWidth, height));
  if (ruler) parts.push(renderSubgridLines(ruler, start, end, trackWidth, height));
  for (const bin of bins) {
    const x0 = xp(bin.start, start, end, trackWidth);
    const x1 = xp(bin.end, start, end, trackWidth);
    const w = Math.max(1, x1 - x0);
    const h = (bin.value / vtotal) * hAvail * peakScale;
    parts.push(
      `<rect ${tipAttr(gnomadBinTipLines(bin, track))} x="${x0.toFixed(2)}" y="${(baseline - h).toFixed(2)}" width="${w.toFixed(2)}" height="${h.toFixed(2)}" fill="${esc(track.fill)}" fill-opacity="0.9"></rect>`
    );
  }
  parts.push("</svg>");
  return parts.join("");
}

function formatBinRange(bin) {
  return `${bin.start.toLocaleString()} – ${bin.end.toLocaleString()}`;
}

function gnomadBinTipLines(bin, track) {
  const lines = [
    formatBinRange(bin),
    `Variants: ${Number(bin.value || 0).toLocaleString()}`,
  ];
  if (bin.snv_count != null) lines.push(`SNVs: ${Number(bin.snv_count).toLocaleString()}`);
  if (bin.indel_count != null) lines.push(`Indels: ${Number(bin.indel_count).toLocaleString()}`);
  if (bin.rare_count != null) lines.push(`Rare: ${Number(bin.rare_count).toLocaleString()}`);
  if (bin.common_count != null) lines.push(`Common: ${Number(bin.common_count).toLocaleString()}`);
  if (track.resolution_bp != null) {
    lines.push(`Bin: ${Number(track.resolution_bp).toLocaleString()} bp`);
  }
  return lines;
}

function clinvarBinTipLines(bin) {
  return [
    formatBinRange(bin),
    `P/LP: ${Number(bin.plp || 0).toLocaleString()}`,
    `VUS / conflicting: ${Number(bin.vus || 0).toLocaleString()}`,
    `Benign: ${Number(bin.benign || 0).toLocaleString()}`,
    `Total: ${Number(bin.total || 0).toLocaleString()}`,
  ];
}

function sampleBinTipLines(bin) {
  return [
    formatBinRange(bin),
    `Known: ${Number(bin.known || 0).toLocaleString()}`,
    `Novel: ${Number(bin.novel || 0).toLocaleString()}`,
    `Total: ${Number(bin.total || 0).toLocaleString()}`,
  ];
}

export function renderClinvarStackedSvg(
  track,
  start,
  end,
  trackWidth,
  height = HEIGHT.histogram,
  ruler = null
) {
  const bins = track.bins || [];
  const baseline = height - 8;
  const top = 10;
  const hAvail = baseline - top;
  const PEAK_FRAC = 0.9;
  const vtotal = Math.max(1, bins.reduce((sum, b) => sum + (b.total || 0), 0));
  const vmax = bins.length ? Math.max(1, ...bins.map((b) => b.total || 0)) : 1;
  const peakScale = (PEAK_FRAC * vtotal) / vmax;
  const parts = [];
  parts.push(svgOpen(trackWidth, height));
  if (ruler) parts.push(renderSubgridLines(ruler, start, end, trackWidth, height));

  for (const bin of bins) {
    const x0 = xp(bin.start, start, end, trackWidth);
    const x1 = xp(bin.end, start, end, trackWidth);
    const w = Math.max(1, x1 - x0);
    const hTotal = (bin.total / vtotal) * hAvail * peakScale;
    if (hTotal <= 0) continue;

    const tip = tipAttr(clinvarBinTipLines(bin));
    let yTop = baseline;
    const segments = [
      { count: bin.plp || 0, fill: CLINVAR_PLP },
      { count: bin.vus || 0, fill: CLINVAR_VUS },
      { count: bin.benign || 0, fill: CLINVAR_BENIGN },
    ];
    for (const seg of segments) {
      if (seg.count <= 0) continue;
      const hSeg = (seg.count / bin.total) * hTotal;
      yTop -= hSeg;
      parts.push(
        `<rect x="${x0.toFixed(2)}" y="${yTop.toFixed(2)}" width="${w.toFixed(2)}" height="${hSeg.toFixed(2)}" ` +
          `fill="${seg.fill}" fill-opacity="0.92" pointer-events="none"></rect>`
      );
    }
    parts.push(
      `<rect ${tip} x="${x0.toFixed(2)}" y="${(baseline - hTotal).toFixed(2)}" width="${w.toFixed(2)}" height="${hTotal.toFixed(2)}" fill="transparent"></rect>`
    );
  }

  parts.push("</svg>");
  return parts.join("");
}

export function renderGeneSvg(track, start, end, trackWidth, height = HEIGHT.gene, ruler = null) {
  const features = track.features || [];
  const cy = height / 2 + 2;
  const exonH = 30;
  const parts = [];
  parts.push(svgOpen(trackWidth, height));
  if (ruler) parts.push(renderSubgridLines(ruler, start, end, trackWidth, height));
  for (const gene of features) {
    const xLeft = xp(Math.max(gene.start, start), start, end, trackWidth);
    const xRight = xp(Math.min(gene.end, end), start, end, trackWidth);
    if (xRight - xLeft < 1) continue;

    parts.push(
      `<line x1="${xLeft.toFixed(2)}" y1="${cy}" x2="${xRight.toFixed(2)}" y2="${cy}" stroke="#000" stroke-width="1.1"/>`
    );

    for (const ex of gene.exons || []) {
      if (ex.end < start || ex.start > end) continue;
      const x0 = xp(Math.max(ex.start, start), start, end, trackWidth);
      const x1 = xp(Math.min(ex.end, end), start, end, trackWidth);
      const w = Math.max(2, x1 - x0);
      parts.push(
        `<rect x="${x0.toFixed(2)}" y="${(cy - exonH / 2).toFixed(2)}" width="${w.toFixed(2)}" height="${exonH}" fill="${EXON_FILL}" stroke="${EXON_STROKE}" stroke-width="0.6"><title>${esc(gene.gene_name)} exon ${ex.start.toLocaleString()}-${ex.end.toLocaleString()}</title></rect>`
      );
    }

    const direction = gene.strand === "-" ? -1 : 1;
    let x = xLeft + 10;
    while (x < xRight - 10) {
      const pts =
        direction > 0
          ? `${x},${cy - 4} ${x + 7},${cy} ${x},${cy + 4}`
          : `${x},${cy - 4} ${x - 7},${cy} ${x},${cy + 4}`;
      parts.push(`<polygon points="${pts}" fill="#000"/>`);
      x += 18;
    }

    const labelX = Math.min(xRight + 8, trackWidth - 48);
    if (xRight - xLeft > 22) {
      parts.push(
        `<text x="${labelX.toFixed(2)}" y="${cy + 5}" font-size="11" font-weight="normal" fill="#000">${esc(gene.gene_name)}</text>`
      );
    }
  }
  parts.push("</svg>");
  return parts.join("");
}

export function renderSampleStackedSvg(
  track,
  start,
  end,
  trackWidth,
  height = HEIGHT.sample,
  ruler = null
) {
  const bins = track.bins || [];
  const baseline = height - 8;
  const top = 10;
  const hAvail = baseline - top;
  const PEAK_FRAC = 0.9;
  const vtotal = Math.max(1, bins.reduce((sum, b) => sum + (b.total || 0), 0));
  const vmax = bins.length ? Math.max(1, ...bins.map((b) => b.total || 0)) : 1;
  const peakScale = (PEAK_FRAC * vtotal) / vmax;
  const parts = [];
  parts.push(svgOpen(trackWidth, height));
  if (ruler) parts.push(renderSubgridLines(ruler, start, end, trackWidth, height));

  for (const bin of bins) {
    const x0 = xp(bin.start, start, end, trackWidth);
    const x1 = xp(bin.end, start, end, trackWidth);
    const w = Math.max(1, x1 - x0);
    const hTotal = (bin.total / vtotal) * hAvail * peakScale;
    if (hTotal <= 0) continue;

    const tip = tipAttr(sampleBinTipLines(bin), "sample-bin-hit");
    let yTop = baseline;
    const segments = [
      { count: bin.known || 0, fill: SAMPLE_KNOWN },
      { count: bin.novel || 0, fill: SAMPLE_NOVEL },
    ];
    for (const seg of segments) {
      if (seg.count <= 0) continue;
      const hSeg = (seg.count / bin.total) * hTotal;
      yTop -= hSeg;
      parts.push(
        `<rect x="${x0.toFixed(2)}" y="${yTop.toFixed(2)}" width="${w.toFixed(2)}" height="${hSeg.toFixed(2)}" ` +
          `fill="${seg.fill}" fill-opacity="0.92" pointer-events="none"></rect>`
      );
    }
    parts.push(
      `<rect ${tip} data-bin-start="${bin.start}" data-bin-end="${bin.end}" ` +
        `x="${x0.toFixed(2)}" y="${(baseline - hTotal).toFixed(2)}" width="${w.toFixed(2)}" height="${hTotal.toFixed(2)}" fill="transparent"></rect>`
    );
  }

  parts.push("</svg>");
  return parts.join("");
}

export function renderTrackCanvas(track, start, end, trackWidth, ruler = null) {
  if (track.type === "gene") {
    return { html: renderGeneSvg(track, start, end, trackWidth, HEIGHT.gene, ruler), height: HEIGHT.gene };
  }
  if (track.type === "clinvar_stacked") {
    return {
      html: renderClinvarStackedSvg(track, start, end, trackWidth, HEIGHT.clinvar, ruler),
      height: HEIGHT.clinvar,
    };
  }
  if (track.type === "sample_stacked" || track.type === "sample_variant") {
    return {
      html: renderSampleStackedSvg(track, start, end, trackWidth, HEIGHT.sample, ruler),
      height: HEIGHT.sample,
    };
  }
  return {
    html: renderHistogramSvg(track, start, end, trackWidth, HEIGHT.histogram, ruler),
    height: HEIGHT.histogram,
  };
}

export function renderIdeogramSection(viewport) {
  const { chrom, start, end, chrom_length, ideogram_cytobands = [] } = viewport;
  const bands = Array.isArray(ideogram_cytobands) && ideogram_cytobands.length ? ideogram_cytobands : [];
  const label = ideogramLabel(chrom, bands, start, end);
  const viewCenter = Math.floor((start + end) / 2);
  const ideogram = renderIdeogramSvg(bands, chrom_length, viewCenter);

  return `
    <div class="ideogram-panel" id="ideogram-panel">
      <div class="ideogram-fixed-wrap">
        <div class="ideogram-inner-row">
          <div class="ideogram-chrom-label">
            <span id="ideogram-label-text" class="ideogram-label-text" data-default-label="${esc(label)}">${esc(label)}</span>
            <span class="ideogram-label-hint" title="Hover band for name; click to jump; drag to select range">ⓘ</span>
          </div>
          <div class="ideogram-canvas" id="ideogram-track">
            ${ideogram.svg}
            <div class="ideogram-band-labels" aria-hidden="true">${ideogram.labels}</div>
          </div>
        </div>
      </div>
    </div>
  `;
}

export function renderTracksSection(viewport, layout) {
  const { start, end, tracks = [], ruler } = viewport;
  const { trackWidth, scrollContentWidth, scrollOffset = 0 } = layout;

  const labelRows = [];
  const canvasRows = [];

  labelRows.push(renderLabel(null, ""));
  tracks.forEach((track) => {
    let countLabel;
    if (track.type === "gene") {
      countLabel = `${track.feature_count ?? 0} genes`;
    } else if (track.type === "sample_stacked" || track.type === "sample_variant") {
      const novelN = track.novel_count ?? 0;
      countLabel = `${track.variant_count ?? 0} variants · ${novelN} novel`;
    } else if (track.type === "histogram" || track.type === "clinvar_stacked") {
      const bins = track.bins || [];
      const resBp =
        track.resolution_bp != null ? Number(track.resolution_bp).toLocaleString() : "?";
      countLabel = `${bins.length} bins · ${resBp} bp`;
    } else {
      countLabel = track.label || track.type || "track";
    }
    const hint = `${esc(track.sublabel)}<br/>${esc(countLabel)}`;
    const isSample = track.type === "sample_stacked" || track.type === "sample_variant";
    labelRows.push(
      renderLabel(track.label, hint, {
        swatchClass: isSample ? " track-swatch-sample" : "",
      })
    );
  });

  canvasRows.push({
    html: renderRulerSvg(ruler, start, end, trackWidth, HEIGHT.ruler),
    height: HEIGHT.ruler,
  });
  tracks.forEach((track) => {
    canvasRows.push(renderTrackCanvas(track, start, end, trackWidth, ruler));
  });

  const scrollLabelHtml = labelRows
    .map((row, i) => {
      const h = canvasRows[i].height;
      return `<div class="label-row" data-track-index="${i}" style="height:${h}px">${row}</div>`;
    })
    .join("");

  const canvasHtml = canvasRows
    .map(
      (row, i) => `
      <div class="canvas-row" data-track-index="${i}" style="min-height:${row.height}px;width:${scrollContentWidth}px">
        <div class="canvas-frame" style="width:${trackWidth}px;margin-left:${scrollOffset}px;min-height:${row.height}px">
          ${row.html}
        </div>
      </div>`
    )
    .join("");

  const panEnabled = scrollContentWidth > trackWidth;

  const sampleDetailHtml = `
        <div id="sample-bin-detail" class="sample-bin-detail" hidden>
          <div class="sample-bin-detail-header" id="sample-bin-detail-header"></div>
          <div class="sample-bin-detail-body" id="sample-bin-detail-body"></div>
        </div>`;

  return `
    <div class="tracks-layout" id="tracks-section">
      <div class="labels-column" style="width:${LABEL_WIDTH}px">${scrollLabelHtml}</div>
      <div class="tracks-main">
        <div class="canvas-scroll-viewport" id="pan-viewport">
          <div class="canvas-scroll-inner" id="pan-scroll-inner" style="width:${scrollContentWidth}px">
            ${canvasHtml}
          </div>
        </div>
        ${sampleDetailHtml}
        <div class="pan-scrollbar-row ${panEnabled ? "" : "is-disabled"}">
          <input type="range" id="pan-slider" class="pan-slider"
            min="0" max="1000" value="0" step="1"
            aria-label="Pan along chromosome" ${panEnabled ? "" : "disabled"} />
        </div>
      </div>
    </div>
  `;
}

export function renderBrowserPanel(viewport, layout) {
  const { chrom, start, end, assembly, window_bp, sample_id: sampleId } = viewport;
  const sampleSuffix = sampleId
    ? ` · <span class="panel-sample-id">sample ${esc(sampleId)}</span>`
    : "";

  return `
    <div class="browser-panel">
      <p class="panel-title">Human ${esc(assembly)} — ${esc(chrom)} exome variant landscape${sampleSuffix}</p>
      <p class="panel-subtitle" id="panel-subtitle">${start.toLocaleString()} – ${end.toLocaleString()} (${window_bp.toLocaleString()} bp)</p>
      ${renderIdeogramSection(viewport)}
      ${renderTracksSection(viewport, layout)}
    </div>
  `;
}

export function updateIdeogramFromViewport(viewport) {
  const { start, end, chrom, chrom_length, ideogram_cytobands = [] } = viewport;
  const viewCenter = Math.floor((start + end) / 2);
  updateIdeogramMarker(viewCenter, ideogram_cytobands, chrom_length);
  const label = ideogramLabel(chrom, ideogram_cytobands, start, end);
  updateIdeogramLabelText(label);
}

export function updatePanelSubtitle(viewport) {
  const el = document.getElementById("panel-subtitle");
  if (!el) return;
  const { start, end, window_bp } = viewport;
  el.textContent = `${start.toLocaleString()} – ${end.toLocaleString()} (${window_bp.toLocaleString()} bp)`;
}
