export const IDEOGRAM_VIEW_WIDTH = 1000;
export const IDEOGRAM_HEIGHT = 27;

const VIEWPORT_RED = "#CC0000";
const ACEN_PINCH_RED = "#8B0000";
const ACEN_PINCH_RED_Q = "#A52A2A";
const BAND_STROKE = "#000000";
const BAND_STROKE_W = 0.45;
const MIN_LABEL_VIEWBOX_W = 22;
const MIN_LABEL_CHAR_W = 6.2;
const BROWSER_FONT = "Helvetica, Arial, sans-serif";

function esc(text) {
  return String(text ?? "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

/** UCSC cytoBand 0-based coordinate → ideogram x (full chromosome fills 0..width). */
export function ideogramX(pos, chromLength, trackWidth = IDEOGRAM_VIEW_WIDTH) {
  const len = Math.max(1, Number(chromLength));
  return (Number(pos) / len) * trackWidth;
}

/** Ideogram x → 1-based genomic position (matches track viewport coords). */
export function ideogramPosAtX(svgX, chromLength, trackWidth = IDEOGRAM_VIEW_WIDTH) {
  const len = Math.max(1, Number(chromLength));
  const w = Math.max(1, trackWidth);
  const pos = Math.round((svgX / w) * len);
  return Math.max(1, Math.min(len, pos));
}

export function cytobandAt(bands, pos) {
  for (const band of bands) {
    if (pos >= band.start && pos < band.end) return band.name;
  }
  return bands.length ? bands[bands.length - 1].name : "";
}

export function ideogramLabel(chrom, bands, viewStart, viewEnd) {
  const center = Math.floor((viewStart + viewEnd) / 2);
  const band = cytobandAt(bands, center);
  return band ? `${chrom} (${band})` : chrom;
}

export function readIdeogramMeta(svg) {
  return {
    chromLength: Number(svg.dataset.chromLength),
  };
}

function bandLabelFits(bw, name) {
  const need = Math.max(MIN_LABEL_VIEWBOX_W, String(name).length * MIN_LABEL_CHAR_W);
  return bw >= need;
}

function bandTitle(band) {
  return `${band.name} (${band.start.toLocaleString()}-${band.end.toLocaleString()})`;
}

function renderIdeogramBandLabelsHtml(bandGeom, trackWidth) {
  return bandGeom
    .filter(({ w: bw, band }) => bandLabelFits(bw, band.name))
    .map(({ x0, w: bw, band }) => {
      const leftPct = ((x0 + bw / 2) / trackWidth) * 100;
      return `<span class="ideogram-band-label" style="left:${leftPct.toFixed(4)}%" title="${esc(bandTitle(band))}">${esc(band.name)}</span>`;
    })
    .join("");
}

/** Draw band edges: rects get top/bottom + verticals; acen gets bowtie outline only. */
function appendIdeogramOutlines(parts, bands, acenSet, acenIndices, junctionBp, len, w, height, barY) {
  const yTop = barY;
  const yBot = barY + height;
  const midY = barY + height / 2;
  const s = `class="ideogram-outline" stroke="${BAND_STROKE}" stroke-width="${BAND_STROKE_W}" fill="none" shape-rendering="crispEdges" pointer-events="none"`;

  for (let i = 0; i < bands.length; i++) {
    if (acenSet.has(i)) continue;

    const band = bands[i];
    const x0 = ideogramX(band.start, len, w);
    const x1 = ideogramX(band.end, len, w);
    const beforeAcen = acenIndices.length >= 2 && i === acenIndices[0] - 1;
    const afterAcen = acenIndices.length >= 2 && i === acenIndices[acenIndices.length - 1] + 1;

    if (beforeAcen) {
      const left = bands[acenIndices[0]];
      const right = bands[acenIndices[acenIndices.length - 1]];
      const x0p = ideogramX(left.start, len, w);
      const x1q = ideogramX(right.end, len, w);
      const cx = ideogramX(junctionBp, len, w);

      parts.push(`<line x1="${x0.toFixed(2)}" y1="${yTop}" x2="${x0p.toFixed(2)}" y2="${yTop}" ${s}/>`);
      parts.push(`<line x1="${x0.toFixed(2)}" y1="${yBot}" x2="${x0p.toFixed(2)}" y2="${yBot}" ${s}/>`);
      if (x0 > 0.01) {
        parts.push(`<line x1="${x0.toFixed(2)}" y1="${yTop}" x2="${x0.toFixed(2)}" y2="${yBot}" ${s}/>`);
      }

      parts.push(
        `<polyline points="${x0p.toFixed(2)},${yTop} ${cx.toFixed(2)},${midY.toFixed(2)} ${x0p.toFixed(2)},${yBot}" ${s}/>`
      );
      parts.push(
        `<polyline points="${cx.toFixed(2)},${midY.toFixed(2)} ${x1q.toFixed(2)},${yTop} ${x1q.toFixed(2)},${yBot}" ${s}/>`
      );
      continue;
    }

    if (afterAcen) {
      parts.push(`<line x1="${x0.toFixed(2)}" y1="${yTop}" x2="${x1.toFixed(2)}" y2="${yTop}" ${s}/>`);
      parts.push(`<line x1="${x0.toFixed(2)}" y1="${yBot}" x2="${x1.toFixed(2)}" y2="${yBot}" ${s}/>`);
      parts.push(`<line x1="${x0.toFixed(2)}" y1="${yTop}" x2="${x0.toFixed(2)}" y2="${yBot}" ${s}/>`);
      continue;
    }

    parts.push(`<line x1="${x0.toFixed(2)}" y1="${yTop}" x2="${x1.toFixed(2)}" y2="${yTop}" ${s}/>`);
    parts.push(`<line x1="${x0.toFixed(2)}" y1="${yBot}" x2="${x1.toFixed(2)}" y2="${yBot}" ${s}/>`);
    if (x0 > 0.01) {
      parts.push(`<line x1="${x0.toFixed(2)}" y1="${yTop}" x2="${x0.toFixed(2)}" y2="${yBot}" ${s}/>`);
    }
  }

  parts.push(`<line x1="0" y1="${yTop}" x2="0" y2="${yBot}" ${s}/>`);
  parts.push(`<line x1="${w}" y1="${yTop}" x2="${w}" y2="${yBot}" ${s}/>`);
}

/** Full-chromosome ideogram; marker updated separately. */
export function renderIdeogramSvg(bands, chromLength, viewCenter) {
  const w = IDEOGRAM_VIEW_WIDTH;
  const height = IDEOGRAM_HEIGHT;
  const barY = 0;
  const len = Math.max(1, Number(chromLength));

  const parts = [];
  parts.push(
    `<svg class="ideogram-svg" viewBox="0 0 ${w} ${height}" preserveAspectRatio="none" ` +
      `data-chrom-length="${len}" font-family="${BROWSER_FONT}" xmlns="http://www.w3.org/2000/svg">`
  );

  const acenIndices = [];
  const bandGeom = [];

  bands.forEach((band, i) => {
    if ((band.stain || "") === "acen") acenIndices.push(i);
  });

  const acenSet = new Set(acenIndices);
  const junctionBp =
    acenIndices.length >= 2 ? Number(bands[acenIndices[0]].end) : null;

  bands.forEach((band, i) => {
    if (acenSet.has(i)) return;

    const x0 = ideogramX(band.start, len, w);
    const x1 = ideogramX(band.end, len, w);
    const bw = Math.max(0.2, x1 - x0);
    const fill = band.fill || "#FFFFFF";
    parts.push(
      `<rect class="cytoband-band ideogram-band" data-band-start="${band.start}" data-band-end="${band.end}" data-band-name="${esc(band.name)}" ` +
        `x="${x0.toFixed(2)}" y="${barY}" width="${bw.toFixed(2)}" height="${height}" fill="${esc(fill)}" ` +
        `stroke="none" shape-rendering="crispEdges">` +
        `<title>${esc(bandTitle(band))}</title></rect>`
    );
    bandGeom.push({ x0, w: bw, band });
  });

  if (acenIndices.length >= 2) {
    const left = bands[acenIndices[0]];
    const right = bands[acenIndices[acenIndices.length - 1]];
    const x0p = ideogramX(left.start, len, w);
    const x1q = ideogramX(right.end, len, w);
    const cx = ideogramX(junctionBp, len, w);
    const midY = barY + height / 2;

    // p11.1 / q11.1: bowtie trapezoids only — no surrounding rectangle.
    parts.push(
      `<polygon class="cytoband-band ideogram-band ideogram-acen-p" data-band-name="${esc(left.name)}" ` +
        `points="${x0p.toFixed(2)},${barY} ${cx.toFixed(2)},${midY.toFixed(2)} ${x0p.toFixed(2)},${barY + height}" ` +
        `fill="${ACEN_PINCH_RED}" stroke="none" shape-rendering="geometricPrecision">` +
        `<title>${esc(bandTitle(left))}</title></polygon>`
    );
    parts.push(
      `<polygon class="cytoband-band ideogram-band ideogram-acen-q" data-band-name="${esc(right.name)}" ` +
        `points="${cx.toFixed(2)},${midY.toFixed(2)} ${x1q.toFixed(2)},${barY} ${x1q.toFixed(2)},${barY + height}" ` +
        `fill="${ACEN_PINCH_RED_Q}" stroke="none" shape-rendering="geometricPrecision">` +
        `<title>${esc(bandTitle(right))}</title></polygon>`
    );

    bandGeom.push({ x0: x0p, w: cx - x0p, band: left });
    bandGeom.push({ x0: cx, w: x1q - cx, band: right });
  }

  appendIdeogramOutlines(parts, bands, acenSet, acenIndices, junctionBp, len, w, height, barY);

  parts.push(
    `<rect id="ideogram-selection-fill" x="0" y="${barY}" width="0" height="${height}" ` +
      `fill="#0C5EC4" fill-opacity="0.28" visibility="hidden" pointer-events="none"/>`
  );
  parts.push(
    `<line id="ideogram-selection-a" x1="0" y1="${barY}" x2="0" y2="${barY + height}" ` +
      `stroke="#0C5EC4" stroke-width="1.2" visibility="hidden" pointer-events="none"/>`
  );
  parts.push(
    `<line id="ideogram-selection-b" x1="0" y1="${barY}" x2="0" y2="${barY + height}" ` +
      `stroke="#0C5EC4" stroke-width="1.2" visibility="hidden" pointer-events="none"/>`
  );

  const markerX = ideogramX(viewCenter, len, w);
  parts.push(
    `<line id="ideogram-marker" x1="${markerX.toFixed(2)}" y1="${barY}" x2="${markerX.toFixed(2)}" y2="${barY + height}" ` +
      `stroke="${VIEWPORT_RED}" stroke-width="1.5" pointer-events="none"/>`
  );

  parts.push("</svg>");
  return {
    svg: parts.join(""),
    labels: renderIdeogramBandLabelsHtml(bandGeom, w),
  };
}

export function updateIdeogramMarker(viewCenter, bands, chromLength) {
  const svg = document.querySelector(".ideogram-svg");
  if (!svg) return;
  const len = Math.max(1, Number(chromLength ?? svg.dataset.chromLength));
  const x = ideogramX(viewCenter, len, IDEOGRAM_VIEW_WIDTH);
  const line = document.getElementById("ideogram-marker");
  if (line) {
    line.setAttribute("x1", x.toFixed(2));
    line.setAttribute("x2", x.toFixed(2));
  }
}

export function updateIdeogramLabelText(text) {
  const el = document.getElementById("ideogram-label-text");
  if (!el) return;
  el.textContent = text;
  el.dataset.defaultLabel = text;
}

export function ideogramClientXToSvgX(svg, clientX) {
  const pt = svg.createSVGPoint();
  pt.x = clientX;
  pt.y = svg.getBoundingClientRect().top + 1;
  const ctm = svg.getScreenCTM();
  if (!ctm) return 0;
  const local = pt.matrixTransform(ctm.inverse());
  return Math.max(0, Math.min(IDEOGRAM_VIEW_WIDTH, local.x));
}

function clampPos(pos, chromLength) {
  return Math.max(1, Math.min(chromLength, pos));
}

function setSelectionLines(x0, x1) {
  const lineA = document.getElementById("ideogram-selection-a");
  const lineB = document.getElementById("ideogram-selection-b");
  const fill = document.getElementById("ideogram-selection-fill");
  if (!lineA || !lineB) return;
  const left = Math.max(0, Math.min(x0, x1));
  const right = Math.max(0, Math.max(x0, x1));
  const span = right - left;

  const placeLine = (line, x) => {
    line.setAttribute("x1", x.toFixed(2));
    line.setAttribute("x2", x.toFixed(2));
    line.setAttribute("visibility", "visible");
  };

  if (span < 1) {
    placeLine(lineA, left);
    lineB.setAttribute("visibility", "hidden");
    if (fill) fill.setAttribute("visibility", "hidden");
    return;
  }

  if (fill) {
    fill.setAttribute("x", left.toFixed(2));
    fill.setAttribute("width", span.toFixed(2));
    fill.setAttribute("visibility", "visible");
  }
  placeLine(lineA, left);
  placeLine(lineB, right);
}

function hideSelectionLines() {
  const lineA = document.getElementById("ideogram-selection-a");
  const lineB = document.getElementById("ideogram-selection-b");
  const fill = document.getElementById("ideogram-selection-fill");
  if (lineA) lineA.setAttribute("visibility", "hidden");
  if (lineB) lineB.setAttribute("visibility", "hidden");
  if (fill) fill.setAttribute("visibility", "hidden");
}

/**
 * Click to center viewport; left-drag shows vertical lines and selects a range.
 */
export function setupIdeogramDrag(trackEl, onNavigate) {
  if (!trackEl) return;
  if (trackEl._dragAbort) trackEl._dragAbort.abort();

  const ac = new AbortController();
  trackEl._dragAbort = ac;
  const { signal } = ac;

  let dragging = false;
  let anchorSvgX = 0;

  const svg = () => trackEl.querySelector(".ideogram-svg");

  const onMouseDown = (ev) => {
    if (ev.button !== 0) return;
    const s = svg();
    if (!s) return;
    dragging = true;
    trackEl.classList.add("ideogram-dragging");
    anchorSvgX = ideogramClientXToSvgX(s, ev.clientX);
    setSelectionLines(anchorSvgX, anchorSvgX);
    ev.preventDefault();
  };

  const onMouseMove = (ev) => {
    if (!dragging) return;
    const s = svg();
    if (!s) return;
    const x = ideogramClientXToSvgX(s, ev.clientX);
    setSelectionLines(anchorSvgX, x);
  };

  const onMouseUp = async (ev) => {
    if (!dragging) return;
    dragging = false;
    trackEl.classList.remove("ideogram-dragging");
    const s = svg();
    if (!s) return;
    const { chromLength } = readIdeogramMeta(s);
    const endSvgX = ideogramClientXToSvgX(s, ev.clientX);
    hideSelectionLines();

    const x0 = Math.min(anchorSvgX, endSvgX);
    const x1 = Math.max(anchorSvgX, endSvgX);
    const minDragPx = 6;

    if (x1 - x0 < minDragPx) {
      const pos = clampPos(ideogramPosAtX(anchorSvgX, chromLength), chromLength);
      const windowBp = onNavigate.getWindowBp?.() ?? chromLength;
      const half = Math.floor(windowBp / 2);
      let start = Math.max(1, pos - half);
      let end = Math.min(chromLength, start + windowBp - 1);
      if (end - start + 1 < windowBp) start = Math.max(1, end - windowBp + 1);
      await onNavigate({ start, end, isClick: true });
      return;
    }

    let start = clampPos(ideogramPosAtX(x0, chromLength), chromLength);
    let end = clampPos(ideogramPosAtX(x1, chromLength), chromLength);
    if (start > end) [start, end] = [end, start];
    if (start === end) end = Math.min(chromLength, start + 1);
    await onNavigate({ start, end, isClick: false });
  };

  trackEl.addEventListener("mousedown", onMouseDown, { signal });
  window.addEventListener("mousemove", onMouseMove, { signal });
  window.addEventListener("mouseup", onMouseUp, { signal });

  trackEl.addEventListener(
    "mouseover",
    (ev) => {
      const band = ev.target.closest?.(".cytoband-band");
      if (!band?.dataset.bandName) return;
      updateIdeogramLabelText(band.dataset.bandName);
    },
    { signal }
  );
  trackEl.addEventListener(
    "mouseleave",
    () => {
      const labelEl = document.getElementById("ideogram-label-text");
      if (labelEl?.dataset.defaultLabel) updateIdeogramLabelText(labelEl.dataset.defaultLabel);
    },
    { signal }
  );
}
