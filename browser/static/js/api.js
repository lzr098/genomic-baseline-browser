export async function fetchManifest() {
  const res = await fetch("/api/manifest");
  if (!res.ok) throw new Error(`manifest ${res.status}`);
  return res.json();
}

export async function fetchChromMeta(chrom) {
  const res = await fetch(`/api/chrom/${encodeURIComponent(chrom)}/meta`);
  if (!res.ok) throw new Error(`meta ${res.status}`);
  return res.json();
}

export async function fetchViewport(chrom, start, end, sample = null) {
  const params = new URLSearchParams({ start: String(start), end: String(end) });
  if (sample) params.set("sample", sample);
  const res = await fetch(`/api/chrom/${encodeURIComponent(chrom)}/viewport?${params}`);
  if (!res.ok) {
    const err = await res.json().catch(() => ({}));
    throw new Error(err.detail || `viewport ${res.status}`);
  }
  return res.json();
}

export async function fetchSampleBinVariants(chrom, sampleId, binStart, binEnd, limit = 100) {
  const params = new URLSearchParams({
    bin_start: String(binStart),
    bin_end: String(binEnd),
    limit: String(limit),
  });
  const res = await fetch(
    `/api/chrom/${encodeURIComponent(chrom)}/sample/${encodeURIComponent(sampleId)}/bin-variants?${params}`
  );
  if (!res.ok) {
    const err = await res.json().catch(() => ({}));
    throw new Error(err.detail || `bin-variants ${res.status}`);
  }
  return res.json();
}

export async function fetchSampleViewportVariants(chrom, sampleId, start, end, limit = 500) {
  const params = new URLSearchParams({
    start: String(start),
    end: String(end),
    limit: String(limit),
  });
  const res = await fetch(
    `/api/chrom/${encodeURIComponent(chrom)}/sample/${encodeURIComponent(sampleId)}/variants?${params}`
  );
  if (!res.ok) {
    const err = await res.json().catch(() => ({}));
    throw new Error(err.detail || `variants ${res.status}`);
  }
  return res.json();
}
