// AnalogFinder static client. Loads per-dimension uint8 score grids and recomputes the weighted overall score
// client-side (same formula as src/similarity.py: weights renormalised over dimensions with data).
"use strict";

const VIRIDIS = [[68,1,84],[72,40,120],[62,73,137],[49,104,142],[38,130,142],[31,158,137],[53,183,121],[110,206,88],[181,222,43],[253,231,37]];
const NODATA = 255;
const MARS_VS_EARTH = [
  ["Gravity", "0.38 g vs 1 g: soil mechanics, dust lofting, slope stability and rover traction all scale differently"],
  ["Pressure", "~6 mbar CO2 vs ~1013 mbar air: no liquid water stability, different aeolian thresholds and heat transfer"],
  ["Radiation", "no global magnetic field and thin atmosphere: surface UV and cosmic-ray dose far above any Earth site"],
  ["Temperature", "Mars surface nights reach ~175–185 K at these sites (model); Earth land rarely gets below ~220 K outside the poles"],
];

const S = { meta: null, target: null, layer: "overall", weights: {}, dims: {}, earth: {}, overlay: null, list: "discovery", markers: null, sel: null };

const $ = (id) => document.getElementById(id);
const fmt = (v, d = 2) => (v === null || v === undefined || Number.isNaN(v) ? "—" : Number(v).toFixed(d));

async function loadBin(path) {
  const r = await fetch(path);
  if (!r.ok) throw new Error(`${path}: ${r.status}`);
  return new Uint8Array(await r.arrayBuffer());
}

function viridis(t) {
  const x = Math.max(0, Math.min(1, t)) * (VIRIDIS.length - 1);
  const i = Math.min(Math.floor(x), VIRIDIS.length - 2), f = x - i;
  const a = VIRIDIS[i], b = VIRIDIS[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}

async function loadTarget(t) {
  S.target = t;
  const dims = S.meta.dims;
  const arrs = await Promise.all(dims.map((d) => loadBin(`data/${t}_${d}.bin`)));
  S.dims = Object.fromEntries(dims.map((d, i) => [d, arrs[i]]));
  const fp = S.meta.targets[t];
  $("target-desc").textContent = `${fp.lat.toFixed(1)}°, ${fp.lon_e360.toFixed(1)}°E (Mars, east-positive 0–360). Fingerprint = median & IQR over a ${fp.window_deg}°×${fp.window_deg}° window.`;
  renderLists();
  render();
  if (S.sel) showCell(S.sel.lat, S.sel.lon, S.sel.title);
}

function overallAt(i) {
  let num = 0, den = 0;
  for (const d of S.meta.dims) {
    const v = S.dims[d][i];
    if (v === NODATA) continue;
    const w = S.weights[d];
    num += w * (v / 254); den += w;
  }
  return den > 0 ? num / den : NaN;
}

function valueAt(i) {
  if (S.layer === "overall") return overallAt(i);
  const v = S.dims[S.layer][i];
  return v === NODATA ? NaN : v / 254;
}

// Draw the equirectangular grid into a Web-Mercator-spaced canvas so the image overlay registers correctly.
function render() {
  const g = S.meta.grid;
  const latTop = g.lat_top, latBot = g.lat_top - g.rows / g.ppd;
  const merc = (lat) => Math.log(Math.tan(Math.PI / 4 + (lat * Math.PI) / 360));
  const yT = merc(latTop), yB = merc(latBot);
  const W = g.cols, H = Math.round(g.rows * 1.4);
  const cv = document.createElement("canvas");
  cv.width = W; cv.height = H;
  const ctx = cv.getContext("2d");
  const img = ctx.createImageData(W, H);
  // stretch colours to the visible range of the current layer for contrast (legend shows 0..1 of that stretch)
  const vals = new Float32Array(g.rows * g.cols);
  for (let i = 0; i < vals.length; i++) vals[i] = valueAt(i);
  for (let y = 0; y < H; y++) {
    const my = yT + ((yB - yT) * (y + 0.5)) / H;
    const lat = (2 * Math.atan(Math.exp(my)) - Math.PI / 2) * 180 / Math.PI;
    const r = Math.min(g.rows - 1, Math.max(0, Math.floor((latTop - lat) * g.ppd)));
    for (let x = 0; x < W; x++) {
      const v = vals[r * g.cols + x];
      const o = (y * W + x) * 4;
      if (Number.isNaN(v)) { img.data[o + 3] = 0; continue; }
      const c = viridis(v);
      img.data[o] = c[0]; img.data[o + 1] = c[1]; img.data[o + 2] = c[2]; img.data[o + 3] = 215;
    }
  }
  ctx.putImageData(img, 0, 0);
  const url = cv.toDataURL();
  const bounds = [[latBot, -180], [latTop, 180]];
  if (S.overlay) S.overlay.setUrl(url); else S.overlay = L.imageOverlay(url, bounds, { opacity: 0.85, interactive: false }).addTo(S.map);
  $("legend-label").textContent = S.layer === "overall" ? "overall similarity (weighted)" : `${S.meta.dim_labels[S.layer]} similarity`;
}

function cellIndex(lat, lon) {
  const g = S.meta.grid;
  const r = Math.floor((g.lat_top - lat) * g.ppd);
  let c = Math.floor((((lon + 180) % 360) + 360) % 360 * g.ppd);
  if (r < 0 || r >= g.rows) return -1;
  return r * g.cols + Math.min(c, g.cols - 1);
}

async function ensureEarth() {
  const keys = Object.keys(S.meta.earth_features);
  const missing = keys.filter((k) => !S.earth[k]);
  const arrs = await Promise.all(missing.map((k) => loadBin(`data/earth_${k}.bin`)));
  missing.forEach((k, i) => (S.earth[k] = arrs[i]));
}

function earthVal(k, i) {
  const f = S.meta.earth_features[k];
  const q = S.earth[k][i];
  return q === NODATA ? NaN : f.lo + (q / 254) * (f.hi - f.lo);
}

function radar(scores) {
  const dims = S.meta.dims, n = dims.length, R = 80;
  const pt = (k, r) => { const a = -Math.PI / 2 + (2 * Math.PI * k) / n; return [r * Math.cos(a), r * Math.sin(a)]; };
  let s = "";
  for (const lv of [0.25, 0.5, 0.75, 1]) s += `<polygon points="${dims.map((_, k) => pt(k, R * lv).join(",")).join(" ")}" fill="none" stroke="#2a3550"/>`;
  dims.forEach((d, k) => {
    const [x, y] = pt(k, R); const [tx, ty] = pt(k, R + 18);
    s += `<line x1="0" y1="0" x2="${x}" y2="${y}" stroke="#2a3550"/>`;
    s += `<text x="${tx}" y="${ty}" text-anchor="middle" dominant-baseline="middle">${S.meta.dim_labels[d].split(" ")[0]}</text>`;
  });
  const poly = dims.map((d, k) => pt(k, R * (Number.isNaN(scores[d]) ? 0 : scores[d])).join(",")).join(" ");
  s += `<polygon points="${poly}" fill="rgba(244,197,66,.35)" stroke="#f4c542" stroke-width="2"/>`;
  $("radar").innerHTML = s;
}

async function showCell(lat, lon, title) {
  S.sel = { lat, lon, title };
  const i = cellIndex(lat, lon);
  $("panel").classList.remove("hidden");
  $("p-title").textContent = title || "Selected cell";
  $("p-coords").textContent = `${lat.toFixed(2)}°, ${lon.toFixed(2)}° (Earth) — vs ${S.meta.targets[S.target].name}`;
  if (i < 0) { $("p-score").textContent = "outside grid"; return; }
  await ensureEarth();
  const scores = Object.fromEntries(S.meta.dims.map((d) => [d, S.dims[d][i] === NODATA ? NaN : S.dims[d][i] / 254]));
  const ov = overallAt(i);
  $("p-score").innerHTML = Number.isNaN(ov) ? "no data <small>(water, ice, or data gap)</small>" : `${fmt(ov)} <small>overall similarity (current weights)</small>`;
  radar(scores);
  const lab = S.meta.dim_labels;
  const valid = S.meta.dims.filter((d) => !Number.isNaN(scores[d]));
  const strong = valid.filter((d) => scores[d] >= 0.6).map((d) => `${lab[d]} (${fmt(scores[d])})`);
  const weak = valid.filter((d) => scores[d] < 0.35).map((d) => `${lab[d]} (${fmt(scores[d])})`);
  let html = "";
  if (!Number.isNaN(ov)) {
    html += `<p><b class="good">Strong match on:</b> ${strong.length ? strong.join(", ") : "none of the five dimensions"}.</p>`;
    html += `<p><b class="bad">Weak on:</b> ${weak.length ? weak.join(", ") : "nothing below 0.35"}.</p>`;
    if (scores.aridity < 0.6) html += `<p class="small muted">Vegetated ground lowers the aridity score, but life is also an asset: biosignature and life-detection instrument tests need biology to find.</p>`;
  }
  html += `<p><b>Earth can't reproduce:</b></p><ul class="small">${MARS_VS_EARTH.map(([a, b]) => `<li><b>${a}:</b> ${b}</li>`).join("")}</ul>`;
  $("p-text").innerHTML = html;
  const fp = S.meta.targets[S.target].features;
  const rows = Object.entries(S.meta.earth_features).map(([k, f]) => {
    const m = fp[k];
    const mars = m && m.median !== null ? (m.iqr ? `${fmt(m.median)} [${fmt(m.q25)}, ${fmt(m.q75)}]` : `${fmt(m.median)}${m.source ? "*" : ""}`) : "—";
    return `<tr><td>${f.label}${f.unit ? ` <small class="muted">${f.unit}</small>` : ""}</td><td>${mars}</td><td>${fmt(earthVal(k, i))}</td></tr>`;
  });
  $("p-table").querySelector("tbody").innerHTML = rows.join("") +
    `<tr><td colspan="3" class="small muted">* Mars thermal extremes from a 1-D thermal model driven by the site's TES thermal inertia and albedo. Thermal inertia is compared as within-planet percentiles (Mars: TES TI; Earth: apparent TI from MODIS).</td></tr>`;
}

function renderLists() {
  const t = S.target;
  const cands = (S.meta.candidates[t] || []).filter((c) => c.list === S.list);
  $("cands").innerHTML = cands.map((c, k) =>
    `<li data-k="${k}"><span class="sc">${fmt(c.mean_score)}</span>${c.country || "?"} <small class="muted">${fmt(c.lat, 1)}, ${fmt(c.lon, 1)} · ${Math.round(c.area_km2).toLocaleString()} km²${c.near_known ? " · near " + c.near_known : ""}</small></li>`).join("");
  $("cands").querySelectorAll("li").forEach((li) => li.addEventListener("click", () => {
    const c = cands[+li.dataset.k];
    S.map.setView([c.lat, c.lon], 6);
    showCell(c.lat, c.lon, `Candidate #${c.rank} — ${c.country || ""}`);
  }));
  const val = S.meta.validation[t] || {};
  $("known").innerHTML = Object.values(val).map((v, k) =>
    `<li data-k="${k}"><span class="sc">${v.percentile === null ? "NaN" : fmt(v.percentile, 0) + "th pct"}</span>${v.name}${v.note ? ` <small class="muted">(${v.note})</small>` : ""}</li>`).join("");
  const vals = Object.values(val);
  $("known").querySelectorAll("li").forEach((li) => li.addEventListener("click", () => {
    const v = vals[+li.dataset.k];
    S.map.setView([v.lat, v.lon], 6);
    showCell(v.lat, v.lon, `Known analog — ${v.name}`);
  }));
  if (S.markers) S.markers.remove();
  S.markers = L.layerGroup().addTo(S.map);
  cands.slice(0, 25).forEach((c) => L.marker([c.lat, c.lon], { icon: L.divIcon({ className: "cand-marker", iconSize: [12, 12] }), title: `#${c.rank} ${c.country}` })
    .on("click", () => showCell(c.lat, c.lon, `Candidate #${c.rank} — ${c.country || ""}`)).addTo(S.markers));
  vals.forEach((v) => L.marker([v.lat, v.lon], { icon: L.divIcon({ className: "known-marker", iconSize: [10, 10] }), title: v.name })
    .on("click", () => showCell(v.lat, v.lon, `Known analog — ${v.name}`)).addTo(S.markers));
}

function buildControls() {
  const sel = $("target");
  sel.innerHTML = Object.entries(S.meta.targets).map(([k, v]) => `<option value="${k}">${v.name}</option>`).join("");
  sel.addEventListener("change", () => loadTarget(sel.value));
  const layers = ["overall", ...S.meta.dims];
  $("layers").innerHTML = layers.map((l) => `<button data-l="${l}" class="${l === S.layer ? "on" : ""}">${l === "overall" ? "Overall" : S.meta.dim_labels[l].split(" ")[0]}</button>`).join("");
  $("layers").querySelectorAll("button").forEach((b) => b.addEventListener("click", () => {
    S.layer = b.dataset.l;
    $("layers").querySelectorAll("button").forEach((x) => x.classList.toggle("on", x === b));
    render();
  }));
  const drawWeights = () => {
    $("weights").innerHTML = S.meta.dims.map((d) => `<label class="w"><span>${S.meta.dim_labels[d]}</span><input type="range" min="0" max="1" step="0.05" value="${S.weights[d]}" data-d="${d}"><output>${fmt(S.weights[d])}</output></label>`).join("");
    $("weights").querySelectorAll("input").forEach((inp) => inp.addEventListener("input", () => {
      S.weights[inp.dataset.d] = +inp.value;
      inp.nextElementSibling.textContent = fmt(+inp.value);
      if (S.layer === "overall") render();
      if (S.sel) showCell(S.sel.lat, S.sel.lon, S.sel.title);
    }));
  };
  drawWeights();
  $("reset-weights").addEventListener("click", () => { S.weights = { ...S.meta.weights }; drawWeights(); render(); });
  document.querySelectorAll(".tabs button").forEach((b) => b.addEventListener("click", () => {
    S.list = b.dataset.list;
    document.querySelectorAll(".tabs button").forEach((x) => x.classList.toggle("on", x === b));
    renderLists();
  }));
  $("close-panel").addEventListener("click", () => { $("panel").classList.add("hidden"); S.sel = null; });
}

async function main() {
  S.meta = await (await fetch("data/meta.json")).json();
  S.weights = { ...S.meta.weights };
  S.map = L.map("map", { worldCopyJump: true, minZoom: 2, maxZoom: 9 }).setView([20, 0], 2);
  L.tileLayer("https://{s}.basemaps.cartocdn.com/dark_nolabels/{z}/{x}/{y}{r}.png", {
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OSM</a> &copy; <a href="https://carto.com/">CARTO</a>', subdomains: "abcd" }).addTo(S.map);
  S.map.createPane("labels"); S.map.getPane("labels").style.zIndex = 650; S.map.getPane("labels").style.pointerEvents = "none";
  L.tileLayer("https://{s}.basemaps.cartocdn.com/dark_only_labels/{z}/{x}/{y}{r}.png", { pane: "labels", subdomains: "abcd" }).addTo(S.map);
  S.map.on("click", (e) => showCell(e.latlng.lat, ((e.latlng.lng + 540) % 360) - 180));
  buildControls();
  await loadTarget(Object.keys(S.meta.targets)[0]);
}

main().catch((e) => { document.body.insertAdjacentHTML("beforeend", `<pre style="position:fixed;bottom:0;left:0;background:#400;color:#fff;padding:8px;z-index:9999">${e}</pre>`); });
