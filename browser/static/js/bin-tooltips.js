/** Multi-line HTML tooltips for histogram / stacked bin rects. */

let tipEl = null;

function ensureTipEl() {
  if (!tipEl) {
    tipEl = document.createElement("div");
    tipEl.id = "bin-tooltip";
    tipEl.className = "bin-tooltip";
    tipEl.hidden = true;
    document.body.appendChild(tipEl);
  }
  return tipEl;
}

export function tipAttr(lines, extraClass = "") {
  const classes = ["bin-hit", extraClass].filter(Boolean).join(" ");
  return `class="${classes}" data-tip="${encodeURIComponent(JSON.stringify(lines))}"`;
}

function parseTip(el) {
  try {
    return JSON.parse(decodeURIComponent(el.dataset.tip || "[]"));
  } catch {
    return [];
  }
}

function positionTip(ev) {
  const tip = ensureTipEl();
  const pad = 12;
  let x = ev.clientX + pad;
  let y = ev.clientY + pad;
  const rect = tip.getBoundingClientRect();
  if (x + rect.width > window.innerWidth - 8) x = ev.clientX - rect.width - pad;
  if (y + rect.height > window.innerHeight - 8) y = ev.clientY - rect.height - pad;
  tip.style.left = `${Math.max(8, x)}px`;
  tip.style.top = `${Math.max(8, y)}px`;
}

function showTip(el, ev) {
  const lines = parseTip(el);
  if (!lines.length) return;
  const tip = ensureTipEl();
  tip.innerHTML = lines.map((line) => `<div class="bin-tooltip-line">${escapeHtml(line)}</div>`).join("");
  tip.hidden = false;
  positionTip(ev);
}

function hideTip() {
  if (tipEl) tipEl.hidden = true;
}

function escapeHtml(text) {
  return String(text ?? "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

export function setupBinTooltips(root = document.getElementById("tracks-section")) {
  if (!root) return;

  if (root._binTipOver) {
    root.removeEventListener("mouseover", root._binTipOver);
    root.removeEventListener("mousemove", root._binTipMove);
    root.removeEventListener("mouseout", root._binTipOut);
  }

  root._binTipOver = (ev) => {
    const hit = ev.target.closest?.(".bin-hit");
    if (!hit) return;
    showTip(hit, ev);
  };
  root._binTipMove = (ev) => {
    const hit = ev.target.closest?.(".bin-hit");
    if (!hit || tipEl?.hidden) return;
    positionTip(ev);
  };
  root._binTipOut = (ev) => {
    const from = ev.target.closest?.(".bin-hit");
    const to = ev.relatedTarget?.closest?.(".bin-hit");
    if (from && from === to) return;
    hideTip();
  };

  root.addEventListener("mouseover", root._binTipOver);
  root.addEventListener("mousemove", root._binTipMove);
  root.addEventListener("mouseout", root._binTipOut);
}
