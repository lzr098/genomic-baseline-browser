export const TRACK_WIDTH = 1048;
export const LABEL_WIDTH = 196;
/** Horizontal inset so ruler labels and track edges are not clipped. */
export const TRACK_PAD_X = 28;

export function xp(pos, start, end, trackW = TRACK_WIDTH, padX = TRACK_PAD_X) {
  const innerW = Math.max(1, trackW - 2 * padX);
  return padX + ((pos - start) / Math.max(1, end - start + 1)) * innerW;
}

export function posAtX(x, start, end, trackW = TRACK_WIDTH, padX = TRACK_PAD_X) {
  const innerW = Math.max(1, trackW - 2 * padX);
  const local = Math.max(0, Math.min(innerW, x - padX));
  return start + (local / innerW) * (end - start + 1);
}

export function formatBp(n) {
  if (n >= 1_000_000) {
    const mb = n / 1_000_000;
    return Number.isInteger(mb) ? `${mb} Mb` : `${mb.toFixed(2)} Mb`;
  }
  if (n >= 1_000) {
    const kb = n / 1_000;
    return Number.isInteger(kb) ? `${kb} kb` : `${kb.toFixed(1)} kb`;
  }
  return `${n} bp`;
}

export function formatBpExact(n) {
  return `${n.toLocaleString("en-US")} bp`;
}

export function clampViewport(start, end, chromLength) {
  let s = Math.max(1, Math.floor(start));
  let e = Math.min(chromLength, Math.floor(end));
  if (s >= e) {
    e = Math.min(chromLength, s + 1);
  }
  return { start: s, end: e };
}

export function centerViewport(center, windowBp, chromLength) {
  const half = Math.floor(windowBp / 2);
  let start = Math.max(1, Math.floor(center - half));
  let end = Math.min(chromLength, start + windowBp - 1);
  if (end - start + 1 < windowBp) {
    start = Math.max(1, end - windowBp + 1);
  }
  return clampViewport(start, end, chromLength);
}

/** Pixel width of the pannable canvas (full chromosome at current zoom). */
export function panContentWidth(trackWidth, windowBp, chromLength) {
  if (windowBp >= chromLength) return trackWidth;
  return Math.max(trackWidth, Math.ceil(trackWidth * (chromLength / windowBp)));
}

/** Map horizontal scroll offset to 1-based viewport start. */
export function startFromPanScroll(scrollLeft, trackWidth, contentWidth, windowBp, chromLength) {
  const maxStart = Math.max(1, chromLength - windowBp + 1);
  if (contentWidth <= trackWidth || maxStart <= 1) return 1;
  const maxScroll = contentWidth - trackWidth;
  const fraction = Math.max(0, Math.min(1, scrollLeft / maxScroll));
  return Math.max(1, Math.min(maxStart, Math.round(1 + fraction * (maxStart - 1))));
}

/** Map viewport start to horizontal scroll offset. */
export function panScrollFromStart(start, trackWidth, contentWidth, windowBp, chromLength) {
  const maxStart = Math.max(1, chromLength - windowBp + 1);
  if (contentWidth <= trackWidth || maxStart <= 1) return 0;
  const maxScroll = contentWidth - trackWidth;
  const fraction = (start - 1) / (maxStart - 1);
  return fraction * maxScroll;
}
