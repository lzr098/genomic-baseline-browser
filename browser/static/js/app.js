import { fetchChromMeta, fetchManifest, fetchViewport } from "./api.js?v=8";
import {
  LABEL_WIDTH,
  centerViewport,
  clampViewport,
  formatBp,
  formatBpExact,
  panContentWidth,
  panScrollFromStart,
  startFromPanScroll,
} from "./coords.js?v=8";
import { clearSampleBinDetail, setupSampleBinDetail } from "./sample-detail.js?v=12";
import { renderBrowserPanel, renderTracksSection, setupBinTooltips, updateIdeogramFromViewport, updatePanelSubtitle } from "./tracks.js?v=10";
import { setupIdeogramDrag } from "./ideogram.js?v=8";

const state = {
  manifest: null,
  chromMeta: null,
  chrom: "chr21",
  start: 1,
  end: 46_709_983,
  chromLength: 46_709_983,
  zoomSteps: [],
  zoomIndex: 0,
  sampleId: null,
};

let suppressPanSync = false;
let ideogramChrom = null;
let panScrollTimer = null;
let resizeObserver = null;

const els = {
  root: document.getElementById("browser-root"),
  chromButtons: document.getElementById("chrom-buttons"),
  locationInput: document.getElementById("location-input"),
  locationGo: document.getElementById("location-go"),
  zoomIn: document.getElementById("zoom-in"),
  zoomOut: document.getElementById("zoom-out"),
  windowBp: document.getElementById("window-bp"),
  status: document.getElementById("status"),
};

function setStatus(text, isError = false) {
  els.status.textContent = text;
  els.status.style.color = isError ? "#b71c1c" : "";
}

function updateLocationInput() {
  els.locationInput.value = `${state.chrom}:${state.start.toLocaleString("en-US")}-${state.end.toLocaleString("en-US")}`;
}

function updateWindowBpLabel() {
  const windowBp = state.end - state.start + 1;
  const label =
    windowBp >= state.chromLength
      ? `${formatBpExact(windowBp)} (full chr)`
      : formatBpExact(windowBp);
  els.windowBp.textContent = label;
}

function getLayoutMetrics() {
  const rootWidth = els.root.clientWidth || window.innerWidth;
  const trackWidth = Math.max(320, rootWidth - LABEL_WIDTH);
  const windowBp = state.end - state.start + 1;
  const scrollContentWidth = panContentWidth(trackWidth, windowBp, state.chromLength);
  const scrollOffset = panScrollFromStart(
    state.start,
    trackWidth,
    scrollContentWidth,
    windowBp,
    state.chromLength
  );
  return { trackWidth, scrollContentWidth, windowBp, scrollOffset };
}

function findZoomIndex(windowBp) {
  const steps = state.zoomSteps || [];
  if (!steps.length) return 0;
  let best = 0;
  let bestDiff = Infinity;
  steps.forEach((step, i) => {
    const diff = Math.abs(step - windowBp);
    if (diff < bestDiff) {
      bestDiff = diff;
      best = i;
    }
  });
  return best;
}

function applyZoomIndex(index) {
  const steps = state.zoomSteps || [];
  if (!steps.length) return;
  state.zoomIndex = Math.max(0, Math.min(index, steps.length - 1));
  const windowBp = steps[state.zoomIndex];
  const center = Math.floor((state.start + state.end) / 2);
  const vp = centerViewport(center, windowBp, state.chromLength);
  state.start = vp.start;
  state.end = vp.end;
}

function normalizeChrom(value) {
  let c = value.trim();
  if (!c.toLowerCase().startsWith("chr")) c = `chr${c}`;
  const suffix = c.slice(3);
  if (suffix.toUpperCase() === "X") return "chrX";
  if (suffix.toUpperCase() === "Y") return "chrY";
  return `chr${suffix}`;
}

function parseLocation(text) {
  let raw = text.trim().replace(/,/g, "");
  if (!raw) throw new Error("empty location");

  let chrom = state.chrom;
  if (raw.includes(":")) {
    const idx = raw.indexOf(":");
    chrom = normalizeChrom(raw.slice(0, idx));
    raw = raw.slice(idx + 1);
  }

  if (raw.includes("-")) {
    const [a, b] = raw.split("-", 2);
    const start = parseInt(a, 10);
    const end = parseInt(b, 10);
    if (!Number.isFinite(start) || !Number.isFinite(end)) throw new Error("invalid range");
    return { chrom, start, end };
  }

  const pos = parseInt(raw, 10);
  if (!Number.isFinite(pos)) throw new Error("invalid position");
  const windowBp = state.end - state.start + 1;
  return { chrom, ...centerViewport(pos, windowBp, state.chromLength) };
}

async function loadChrom(chrom) {
  state.chrom = chrom;
  state.chromMeta = await fetchChromMeta(chrom);
  state.chromLength = state.chromMeta.length;
  state.zoomSteps = state.chromMeta.zoom_steps_bp || [];
  state.zoomIndex = 0;
  state.start = 1;
  state.end = state.chromLength;
  renderChromButtons();
}

function syncPanControls(trackWidth, scrollContentWidth, windowBp) {
  const viewport = document.getElementById("pan-viewport");
  const slider = document.getElementById("pan-slider");
  if (!viewport || !slider) return;

  const panEnabled = scrollContentWidth > trackWidth;
  slider.disabled = !panEnabled;

  if (!panEnabled) {
    slider.value = "0";
    viewport.scrollLeft = 0;
    return;
  }

  suppressPanSync = true;
  const scrollLeft = panScrollFromStart(
    state.start,
    trackWidth,
    scrollContentWidth,
    windowBp,
    state.chromLength
  );
  viewport.scrollLeft = scrollLeft;
  const maxScroll = scrollContentWidth - trackWidth;
  slider.value = String(maxScroll > 0 ? Math.round((scrollLeft / maxScroll) * 1000) : 0);
  requestAnimationFrame(() => {
    suppressPanSync = false;
  });
}

function setupPanControls(layout) {
  const viewport = document.getElementById("pan-viewport");
  const slider = document.getElementById("pan-slider");
  if (!viewport || !slider) return;

  const { trackWidth, scrollContentWidth, windowBp } = layout;
  syncPanControls(trackWidth, scrollContentWidth, windowBp);

  const onPanChange = async () => {
    if (suppressPanSync) return;
    if (scrollContentWidth <= trackWidth) return;

    let scrollLeft = viewport.scrollLeft;
    if (document.activeElement === slider) {
      const maxScroll = scrollContentWidth - trackWidth;
      scrollLeft = (parseInt(slider.value, 10) / 1000) * maxScroll;
      suppressPanSync = true;
      viewport.scrollLeft = scrollLeft;
      requestAnimationFrame(() => {
        suppressPanSync = false;
      });
    }

    const newStart = startFromPanScroll(
      scrollLeft,
      trackWidth,
      scrollContentWidth,
      windowBp,
      state.chromLength
    );
    if (newStart === state.start) return;

    state.start = newStart;
    state.end = Math.min(state.chromLength, newStart + windowBp - 1);
    clearTimeout(panScrollTimer);
    panScrollTimer = setTimeout(() => {
      refreshViewport({ tracksOnly: true });
    }, 100);
  };

  viewport.onscroll = onPanChange;
  slider.oninput = onPanChange;
}

function setupResizeObserver() {
  if (resizeObserver) resizeObserver.disconnect();
  resizeObserver = new ResizeObserver(() => {
    clearTimeout(panScrollTimer);
    panScrollTimer = setTimeout(() => refreshViewport({ tracksOnly: true }), 150);
  });
  resizeObserver.observe(els.root);
}

function setupIdeogramInteraction() {
  const track = document.getElementById("ideogram-track");
  if (!track) return;

  const onNavigate = async ({ start, end, isClick = false }) => {
    state.start = start;
    state.end = end;
    state.zoomIndex = findZoomIndex(end - start + 1);
    updateLocationInput();
    await refreshViewport({ tracksOnly: !isClick });
  };
  onNavigate.getWindowBp = () => state.end - state.start + 1;

  setupIdeogramDrag(track, onNavigate);
}

function setupTrackHover() {
  const section = document.getElementById("tracks-section");
  if (!section) return;

  section.querySelectorAll(".label-row").forEach((labelRow) => {
    const idx = labelRow.dataset.trackIndex;
    const canvasRow = section.querySelector(`.canvas-row[data-track-index="${idx}"]`);
    if (!canvasRow) return;

    const onEnter = () => {
      canvasRow.classList.add("track-highlight");
      labelRow.classList.add("track-label-hover");
    };
    const onLeave = (ev) => {
      const related = ev.relatedTarget;
      if (related && (labelRow.contains(related) || canvasRow.contains(related))) return;
      canvasRow.classList.remove("track-highlight");
      labelRow.classList.remove("track-label-hover");
    };

    labelRow.addEventListener("mouseenter", onEnter);
    labelRow.addEventListener("mouseleave", onLeave);
    canvasRow.addEventListener("mouseenter", onEnter);
    canvasRow.addEventListener("mouseleave", onLeave);
  });
}

function renderTracksArea(viewport, layout) {
  const section = document.getElementById("tracks-section");
  if (section) {
    section.outerHTML = renderTracksSection(viewport, layout);
  }
}

function sampleDetailContext() {
  return {
    getChrom: () => state.chrom,
    getSampleId: () => state.sampleId,
  };
}

async function refreshViewport({ tracksOnly = false } = {}) {
  setStatus("Loading…");
  try {
    const viewport = await fetchViewport(state.chrom, state.start, state.end, state.sampleId);
    const layout = getLayoutMetrics();
    const chromChanged = ideogramChrom !== viewport.chrom;
    const canReuseIdeogram = tracksOnly && !chromChanged && document.getElementById("ideogram-panel");

    clearSampleBinDetail();

    if (canReuseIdeogram) {
      updateIdeogramFromViewport(viewport);
      updatePanelSubtitle(viewport);
      renderTracksArea(viewport, layout);
      setupPanControls(layout);
      setupTrackHover();
      setupBinTooltips();
      setupSampleBinDetail(sampleDetailContext());
    } else {
      els.root.innerHTML = renderBrowserPanel(viewport, layout);
      ideogramChrom = viewport.chrom;
      setupIdeogramInteraction();
      setupPanControls(layout);
      setupTrackHover();
      setupBinTooltips();
      setupSampleBinDetail(sampleDetailContext());
    }

    if (!resizeObserver) setupResizeObserver();
    updateLocationInput();
    updateWindowBpLabel();
    state.zoomIndex = findZoomIndex(viewport.window_bp);
    const sampleTrack = viewport.tracks?.find(
      (t) => t.type === "sample_stacked" || t.type === "sample_variant"
    );
    const sampleNote = sampleTrack
      ? ` · ${sampleTrack.label} ${sampleTrack.variant_count} vars`
      : "";
    const res = viewport.resolutions || {};
    const gnomadKb = res.gnomad_bp != null ? res.gnomad_bp / 1000 : "?";
    setStatus(
      `${viewport.chrom} · ${formatBp(viewport.window_bp)} · gnomAD ${gnomadKb}kb bins${sampleNote}`
    );
  } catch (err) {
    els.root.innerHTML = `<div class="error">${err.message}</div>`;
    setStatus(err.message, true);
  }
}

function chromSortKey(chrom) {
  const suffix = chrom.replace(/^chr/i, "");
  if (suffix.toUpperCase() === "X") return 23;
  if (suffix.toUpperCase() === "Y") return 24;
  const n = parseInt(suffix, 10);
  return Number.isFinite(n) ? n : 99;
}

function renderChromButtons() {
  const chroms = (state.manifest?.chromosomes || [])
    .map((c) => c.chrom)
    .sort((a, b) => chromSortKey(a) - chromSortKey(b));
  els.chromButtons.innerHTML = "";
  for (const chrom of chroms) {
    const btn = document.createElement("button");
    btn.type = "button";
    btn.textContent = chrom.replace("chr", "");
    btn.className = chrom === state.chrom ? "active" : "";
    btn.addEventListener("click", async () => {
      if (chrom === state.chrom) return;
      await loadChrom(chrom);
      await refreshViewport();
    });
    els.chromButtons.appendChild(btn);
  }
}

function wireControls() {
  els.locationGo.addEventListener("click", async () => {
    try {
      const parsed = parseLocation(els.locationInput.value);
      const chromChanged = parsed.chrom !== state.chrom;
      if (chromChanged) {
        await loadChrom(parsed.chrom);
      }
      const vp = clampViewport(parsed.start, parsed.end, state.chromLength);
      state.start = vp.start;
      state.end = vp.end;
      state.zoomIndex = findZoomIndex(state.end - state.start + 1);
      await refreshViewport({ tracksOnly: !chromChanged });
    } catch (err) {
      setStatus(`Location error: ${err.message}`, true);
    }
  });

  els.locationInput.addEventListener("keydown", (ev) => {
    if (ev.key === "Enter") els.locationGo.click();
  });

  els.zoomIn.addEventListener("click", async () => {
    if (state.zoomIndex >= state.zoomSteps.length - 1) return;
    applyZoomIndex(state.zoomIndex + 1);
    await refreshViewport({ tracksOnly: true });
  });

  els.zoomOut.addEventListener("click", async () => {
    if (state.zoomIndex <= 0) return;
    applyZoomIndex(state.zoomIndex - 1);
    await refreshViewport({ tracksOnly: true });
  });
}

async function init() {
  wireControls();
  try {
    state.manifest = await fetchManifest();
    state.sampleId = state.manifest.default_sample || null;
    const defaultChrom = state.manifest.default_chrom || "chr21";
    await loadChrom(defaultChrom);
    await refreshViewport();
  } catch (err) {
    els.root.innerHTML = `<div class="error">${err.message}</div>`;
    setStatus(err.message, true);
  }
}

init();
