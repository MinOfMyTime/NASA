// AnalogFinder static client. Loads per-dimension uint8 score grids and recomputes the weighted overall score
// client-side (same rule as src/similarity.py: weights renormalised, but a cell missing any weighted dimension is n/a).
"use strict";

const VIRIDIS = [[68,1,84],[72,40,120],[62,73,137],[49,104,142],[38,130,142],[31,158,137],[53,183,121],[110,206,88],[181,222,43],[253,231,37]];
const NODATA = 255;
const MARS_VS_EARTH = [
  ["Gravity", "0.38 g vs 1 g: soil mechanics, dust lofting, slope stability and rover traction all scale differently"],
  ["Pressure", "~6 mbar CO2 vs ~1013 mbar air: no liquid water stability, different aeolian thresholds and heat transfer"],
  ["Radiation", "no global magnetic field and thin atmosphere: surface UV and cosmic-ray dose far above any Earth site"],
];
const REV_DIMS = { topo: ["log_relief", "slope_deg", "log_rough", "tpi_m"], thermo: ["pct_ti", "albedo"], mineral: ["comp_mafic", "comp_silica", "comp_fines"] };

const S = { meta: null, target: null, layer: "overall", weights: {}, dims: {}, gap: null, earth: {}, mars: {}, overlay: null,
  list: "discovery", markers: null, sel: null, mode: "default", marsMap: null, marsOverlay: null, marsMarkers: null, tourTimer: null };

const $ = (id) => document.getElementById(id);
const fmt = (v, d = 2) => (v === null || v === undefined || Number.isNaN(v) ? "n/a" : Number(v).toFixed(d));

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
  const arrs = await Promise.all([...dims.map((d) => loadBin(`data/${t}_${d}.bin`)), loadBin(`data/${t}_gap.bin`)]);
  S.dims = Object.fromEntries(dims.map((d, i) => [d, arrs[i]]));
  S.gap = arrs[dims.length];
  const fp = S.meta.targets[t];
  $("target-desc").textContent = `${fp.lat.toFixed(1)}°, ${fp.lon_e360.toFixed(1)}°E (Mars, east-positive 0–360). Fingerprint = median & IQR over a ${fp.window_deg}°×${fp.window_deg}° window.`;
  renderLists();
  render();
  if (S.sel) showCell(S.sel.lat, S.sel.lon, S.sel.title);
}

function overallAt(i) {
  let num = 0, den = 0;
  for (const d of S.meta.dims) {
    const w = S.weights[d];
    const v = S.dims[d][i];
    if (v === NODATA) { if (w > 0) return NaN; continue; }
    num += w * (v / 254); den += w;
  }
  return den > 0 ? num / den : NaN;
}

function valueAt(i) {
  if (S.layer === "overall") return overallAt(i);
  if (S.layer === "gap") return S.gap[i] === NODATA ? NaN : S.gap[i] / 254;
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
  $("legend-label").textContent = S.layer === "overall" ? "overall similarity (weighted; n/a where any dimension is missing)"
    : S.layer === "gap" ? "Mars gap: absolute thermal-extremes match (not weighted)" : `${S.meta.dim_labels[S.layer]} similarity`;
}

function cellIndex(lat, lon) {
  const g = S.meta.grid;
  const r = Math.floor((g.lat_top - lat) * g.ppd);
  const c = Math.floor((((lon + 180) % 360) + 360) % 360 * g.ppd);
  if (r < 0 || r >= g.rows) return -1;
  return r * g.cols + Math.min(c, g.cols - 1);
}

async function ensureEarth() {
  const keys = Object.keys(S.meta.earth_features);
  const missing = keys.filter((k) => !S.earth[k]);
  const arrs = await Promise.all(missing.map((k) => loadBin(`data/earth_${k}.bin`)));
  missing.forEach((k, i) => (S.earth[k] = arrs[i]));
}

async function ensureMars() {
  const keys = Object.keys(S.meta.mars.features);
  const missing = keys.filter((k) => !S.mars[k]);
  const arrs = await Promise.all(missing.map((k) => loadBin(`data/mars_${k}.bin`)));
  missing.forEach((k, i) => (S.mars[k] = arrs[i]));
}

function earthVal(k, i) {
  const f = S.meta.earth_features[k];
  const q = S.earth[k][i];
  return q === NODATA ? NaN : f.lo + (q / 254) * (f.hi - f.lo);
}
function marsVal(k, i) {
  const f = S.meta.mars.features[k];
  const q = S.mars[k][i];
  return q === NODATA ? NaN : f.lo + (q / 254) * (f.hi - f.lo);
}

function radar(scores) {
  const dims = S.meta.dims, n = dims.length, R = 80;
  const pt = (k, r) => { const a = -Math.PI / 2 + (2 * Math.PI * k) / n; return [r * Math.cos(a), r * Math.sin(a)]; };
  let s = "";
  for (const lv of [0.25, 0.5, 0.75, 1]) s += `<polygon points="${dims.map((_, k) => pt(k, R * lv).join(",")).join(" ")}" fill="none" stroke="#2a3550"/>`;
  dims.forEach((d, k) => {
    const [x, y] = pt(k, R); const [tx, ty] = pt(k, R + 18);
    const lab = d === "extremes" ? "Extremes (rel.)" : S.meta.dim_labels[d].split(" ")[0];
    s += `<line x1="0" y1="0" x2="${x}" y2="${y}" stroke="#2a3550"/>`;
    s += `<text x="${tx}" y="${ty}" text-anchor="middle" dominant-baseline="middle">${lab}</text>`;
  });
  const poly = dims.map((d, k) => pt(k, R * (Number.isNaN(scores[d]) ? 0 : scores[d])).join(",")).join(" ");
  s += `<polygon points="${poly}" fill="rgba(244,197,66,.35)" stroke="#f4c542" stroke-width="2"/>`;
  $("radar").innerHTML = s;
}

function gapBars(i) {
  const fp = S.meta.targets[S.target].features;
  const dtE = earthVal("diurnal_range_k", i), tminE = earthVal("lst_min_k", i);
  const dtM = fp.diurnal_range_k.median, tminM = fp.lst_min_k.median;
  const gapScore = S.gap[i] === NODATA ? NaN : S.gap[i] / 254;
  const ratio = Number.isNaN(dtE) ? NaN : Math.min(1, dtE / dtM);
  const cold = Number.isNaN(tminE) ? NaN : Math.max(0, Math.min(1, (320 - tminE) / (320 - tminM)));
  const bar = (label, frac, text) => `<div class="row"><span>${label}</span><div class="track"><div class="fill" style="width:${Number.isNaN(frac) ? 0 : frac * 100}%"></div></div><span>${text}</span></div>`;
  $("p-gap").innerHTML = `<div class="gap"><b>Mars gap</b> <span class="muted">(absolute, not in the weighted score; score ${fmt(gapScore)})</span>
    ${bar("Diurnal range", ratio, `${fmt(dtE, 0)} / ${fmt(dtM, 0)} K`)}
    ${bar("Coldest night", cold, `${fmt(tminE, 0)} vs ${fmt(tminM, 0)} K`)}
    <span class="muted">Earth reaches ${Number.isNaN(ratio) ? "n/a" : Math.round(ratio * 100) + "%"} of the Mars diurnal range here; nights are ${Number.isNaN(tminE) ? "n/a" : Math.round(tminE - tminM) + " K warmer"} than the Mars target (1-D model).</span></div>`;
}

async function showCell(lat, lon, title, extra) {
  S.sel = { lat, lon, title, extra };
  const i = cellIndex(lat, lon);
  $("panel").classList.remove("hidden");
  $("mars-wrap").classList.add("hidden");
  $("p-title").textContent = title || "Selected cell";
  $("p-coords").textContent = `${lat.toFixed(2)}°, ${lon.toFixed(2)}° (Earth) — vs ${S.meta.targets[S.target].name}`;
  $("p-mode").textContent = S.mode === "default" ? "" : `Mission mode: ${S.meta.mission_modes[S.mode].label}` + (extra && extra.mode ? ` · listed under this mode` : "");
  if (i < 0) { $("p-score").textContent = "outside grid"; return; }
  await ensureEarth();
  const scores = Object.fromEntries(S.meta.dims.map((d) => [d, S.dims[d][i] === NODATA ? NaN : S.dims[d][i] / 254]));
  const ov = overallAt(i);
  const missing = S.meta.dims.filter((d) => Number.isNaN(scores[d]) && S.weights[d] > 0);
  $("p-score").innerHTML = Number.isNaN(ov)
    ? `n/a <small>${missing.length ? `(${missing.map((d) => S.meta.dim_labels[d]).join(", ")} not scoreable here: ${missing.includes("mineral") ? "no usable TIR mineral signal" : "water, ice, polar ATI gap or data gap"})` : "(water, ice, or data gap)"}</small>`
    : `${fmt(ov)} <small>overall similarity (current weights)</small>`;
  radar(scores);
  gapBars(i);
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
  html += `<p><b>Earth can't reproduce:</b></p><ul class="small">${MARS_VS_EARTH.map(([a, b]) => `<li><b>${a}:</b> ${b}</li>`).join("")}<li><b>Temperature:</b> see the Mars gap bars above.</li></ul>`;
  $("p-text").innerHTML = html;
  const fp = S.meta.targets[S.target].features;
  const rows = Object.entries(S.meta.earth_features).map(([k, f]) => {
    const m = fp[k];
    const mars = m && m.median !== null ? (m.iqr ? `${fmt(m.median)} [${fmt(m.q25)}, ${fmt(m.q75)}]` : `${fmt(m.median)}${m.source ? "*" : ""}`) : "—";
    return `<tr><td>${f.label}${f.unit ? ` <small class="muted">${f.unit}</small>` : ""}</td><td>${mars}</td><td>${fmt(earthVal(k, i))}</td></tr>`;
  });
  $("p-table").querySelector("tbody").innerHTML = rows.join("") +
    `<tr><td colspan="3" class="small muted">* Mars thermal extremes from a 1-D thermal model driven by the site's TES thermal inertia and albedo. Thermal inertia and composition are compared as within-planet percentiles. "Thermal extremes (rel.)" in the radar is Earth-relative: how Mars-like this cell is compared with the rest of Earth's land.</td></tr>`;
}

// ---------------------------------------------------------------- reverse search: where on Mars?
function marsCellIndex(lat, lonE360) {
  const g = S.meta.mars.grid;
  const r = Math.floor((g.lat_top - lat) * g.ppd);
  const c = Math.floor(((lonE360 % 360) + 360) % 360 * g.ppd);
  return Math.min(Math.max(r, 0), g.rows - 1) * g.cols + Math.min(c, g.cols - 1);
}

function hav(lat1, lon1, lat2, lon2) {
  const p = Math.PI / 180, a = Math.sin((lat2 - lat1) * p / 2) ** 2 + Math.cos(lat1 * p) * Math.cos(lat2 * p) * Math.sin((lon2 - lon1) * p / 2) ** 2;
  return 2 * Math.asin(Math.sqrt(a)) / p;
}

async function reverseSearch() {
  if (!S.sel) return;
  await ensureEarth(); await ensureMars();
  const ei = cellIndex(S.sel.lat, S.sel.lon);
  const g = S.meta.mars.grid, scale = S.meta.mars.reverse_scale;
  const target = {};
  for (const k of Object.keys(S.meta.mars.features)) target[k] = earthVal(k, ei);
  const n = g.rows * g.cols;
  const ov = new Float32Array(n);
  const wd = { topo: S.weights.topo, thermo: S.weights.thermo, mineral: S.weights.mineral };
  for (let i = 0; i < n; i++) {
    let num = 0, den = 0, bad = false;
    for (const [d, feats] of Object.entries(REV_DIMS)) {
      let d2 = 0, m = 0;
      for (const k of feats) {
        const mv = marsVal(k, i), ev = target[k];
        if (Number.isNaN(mv) || Number.isNaN(ev)) { bad = true; break; }
        const z = (mv - ev) / scale[k]; d2 += z * z; m++;
      }
      if (bad) break;
      const s = Math.exp(-(d2 / m) / 2);
      num += wd[d] * s; den += wd[d];
    }
    ov[i] = bad || den === 0 ? NaN : num / den;
  }
  // draw (equirectangular; the Mars map uses EPSG:4326 so no warp needed)
  const cv = document.createElement("canvas"); cv.width = g.cols; cv.height = g.rows;
  const ctx = cv.getContext("2d"); const img = ctx.createImageData(g.cols, g.rows);
  for (let i = 0; i < n; i++) { const o = i * 4; if (Number.isNaN(ov[i])) { img.data[o + 3] = 0; continue; } const c = viridis(ov[i]); img.data[o] = c[0]; img.data[o + 1] = c[1]; img.data[o + 2] = c[2]; img.data[o + 3] = 170; }
  ctx.putImageData(img, 0, 0);
  $("mars-wrap").classList.remove("hidden");
  if (!S.marsMap) {
    S.marsMap = L.map("mars-map", { crs: L.CRS.EPSG4326, minZoom: -1.5, maxZoom: 4, zoomSnap: 0.25, attributionControl: false });
    L.imageOverlay("data/mars_hillshade.png", [[-90, 0], [90, 360]]).addTo(S.marsMap);
    S.marsMap.fitBounds([[-90, 0], [90, 360]]);
    S.marsMarkers = L.layerGroup().addTo(S.marsMap);
  }
  if (S.marsOverlay) S.marsOverlay.setUrl(cv.toDataURL()); else S.marsOverlay = L.imageOverlay(cv.toDataURL(), [[-90, 0], [90, 360]]).addTo(S.marsMap);
  setTimeout(() => { S.marsMap.invalidateSize(); S.marsMap.fitBounds([[-90, 0], [90, 360]]); }, 50);
  // top 3 distinct locations (>= 15 deg apart), labelled by nearest named feature
  const order = Array.from(ov.keys()).filter((i) => !Number.isNaN(ov[i])).sort((a, b) => ov[b] - ov[a]);
  const picks = [];
  for (const i of order) {
    const lat = g.lat_top - (Math.floor(i / g.cols) + 0.5) / g.ppd, lon = (i % g.cols + 0.5) / g.ppd;
    if (picks.every((p) => hav(lat, lon, p.lat, p.lon) >= 15)) picks.push({ lat, lon, s: ov[i] });
    if (picks.length === 3) break;
  }
  S.marsMarkers.clearLayers();
  $("mars-top").innerHTML = picks.map((p, k) => {
    let best = null, bd = 1e9;
    for (const f of S.meta.mars.named) { const d = hav(p.lat, p.lon, f.lat, f.lon_e360); if (d < bd) { bd = d; best = f; } }
    L.marker([p.lat, p.lon], { icon: L.divIcon({ className: "mars-marker", iconSize: [12, 12] }) }).addTo(S.marsMarkers);
    return `<li><span class="sc">${fmt(p.s)}</span>${best.name} <small class="muted">(${p.lat.toFixed(1)}°, ${p.lon.toFixed(1)}°E; ${bd.toFixed(0)}° from the feature centre)</small></li>`;
  }).join("");
}

// ---------------------------------------------------------------- lists, modes, controls
function passesMode(c) {
  const f = S.meta.mission_modes[S.mode].filters || {};
  const i = cellIndex(c.lat, c.lon);
  if (f.slope_deg_max !== undefined && S.earth.slope_deg && !(earthVal("slope_deg", i) <= f.slope_deg_max)) return false;
  if (f.ndvi_min !== undefined && S.earth.ndvi && !(earthVal("ndvi", i) >= f.ndvi_min && earthVal("ndvi", i) <= f.ndvi_max)) return false;
  return true;
}

function badge(c) {
  const p = c.robust_pool_pct, t = c.robust_top25_pct;
  if (p === null || p === undefined) return `<span class="badge">n/a</span>`;
  const cls = p >= 80 ? "hi" : p >= 50 ? "mid" : "lo";
  return `<span class="badge ${cls}" title="stays a candidate in ${p}% of 500 weight samples; in the top 25 in ${t}%">${Math.round(p)}% (${Math.round(t)})</span>`;
}

function renderLists() {
  const t = S.target;
  const all = (S.meta.candidates[t] || []).filter((c) => c.list === S.list);
  const cands = all.filter(passesMode);
  $("cands").innerHTML = cands.map((c, k) =>
    `<li data-k="${k}">${badge(c)}<span class="sc">${fmt(c.mean_score)}</span>${c.country || "?"} <small class="muted">${fmt(c.lat, 1)}, ${fmt(c.lon, 1)} · ${Math.round(c.area_km2).toLocaleString()} km²${c.near_known ? " · near " + c.near_known : ""}</small></li>`).join("")
    + (cands.length < all.length ? `<li class="muted small">${all.length - cands.length} candidates hidden by the ${S.meta.mission_modes[S.mode].label} filter</li>` : "");
  $("cands").querySelectorAll("li[data-k]").forEach((li) => li.addEventListener("click", () => {
    const c = cands[+li.dataset.k];
    S.map.setView([c.lat, c.lon], 6);
    showCell(c.lat, c.lon, `Candidate #${c.rank} — ${c.country || ""}`, { mode: S.mode });
  }));
  const fillSites = (id, val, prefix) => {
    const vals = Object.values(val || {});
    $(id).innerHTML = vals.map((v, k) =>
      `<li data-k="${k}"><span class="sc">${v.percentile === null ? "n/a" : fmt(v.percentile, 0) + "th pct"}</span>${v.name}${v.note ? ` <small class="muted">(${v.note})</small>` : ""}</li>`).join("");
    $(id).querySelectorAll("li").forEach((li) => li.addEventListener("click", () => {
      const v = vals[+li.dataset.k];
      S.map.setView([v.lat, v.lon], 6);
      showCell(v.lat, v.lon, `${prefix} — ${v.name}`);
    }));
    return vals;
  };
  const known = fillSites("known", S.meta.validation[t], "Known analog");
  const blind = fillSites("blind", S.meta.blind_validation[t], "Blind-set analog");
  if (S.markers) S.markers.remove();
  S.markers = L.layerGroup().addTo(S.map);
  cands.slice(0, 25).forEach((c) => L.marker([c.lat, c.lon], { icon: L.divIcon({ className: "cand-marker", iconSize: [12, 12] }), title: `#${c.rank} ${c.country}` })
    .on("click", () => showCell(c.lat, c.lon, `Candidate #${c.rank} — ${c.country || ""}`, { mode: S.mode })).addTo(S.markers));
  [...known, ...blind].forEach((v) => L.marker([v.lat, v.lon], { icon: L.divIcon({ className: "known-marker", iconSize: [10, 10] }), title: v.name })
    .on("click", () => showCell(v.lat, v.lon, `Known analog — ${v.name}`)).addTo(S.markers));
}

function drawWeights() {
  $("weights").innerHTML = S.meta.dims.map((d) => `<label class="w"><span>${d === "extremes" ? "Extremes (relative)" : S.meta.dim_labels[d]}</span><input type="range" min="0" max="1" step="0.05" value="${S.weights[d]}" data-d="${d}"><output>${fmt(S.weights[d])}</output></label>`).join("");
  $("weights").querySelectorAll("input").forEach((inp) => inp.addEventListener("input", () => {
    S.weights[inp.dataset.d] = +inp.value;
    inp.nextElementSibling.textContent = fmt(+inp.value);
    if (S.layer === "overall") render();
    if (S.sel) showCell(S.sel.lat, S.sel.lon, S.sel.title, S.sel.extra);
  }));
}

function setMode(m) {
  S.mode = m;
  const mm = S.meta.mission_modes[m];
  S.weights = { ...(mm.weights || S.meta.weights) };
  $("mode-note").textContent = mm.note || "";
  drawWeights();
  ensureEarth().then(() => { renderLists(); render(); if (S.sel) showCell(S.sel.lat, S.sel.lon, S.sel.title, S.sel.extra); });
}

function buildControls() {
  const sel = $("target");
  sel.innerHTML = Object.entries(S.meta.targets).map(([k, v]) => `<option value="${k}">${v.name}</option>`).join("");
  sel.addEventListener("change", () => loadTarget(sel.value));
  const mode = $("mode");
  mode.innerHTML = Object.entries(S.meta.mission_modes).map(([k, v]) => `<option value="${k}">${v.label}</option>`).join("");
  mode.addEventListener("change", () => setMode(mode.value));
  const layers = ["overall", ...S.meta.dims, "gap"];
  $("layers").innerHTML = layers.map((l) => `<button data-l="${l}" class="${l === S.layer ? "on" : ""}">${l === "overall" ? "Overall" : l === "gap" ? "Mars gap" : l === "extremes" ? "Extremes (rel.)" : S.meta.dim_labels[l].split(" ")[0]}</button>`).join("");
  $("layers").querySelectorAll("button").forEach((b) => b.addEventListener("click", () => {
    S.layer = b.dataset.l;
    $("layers").querySelectorAll("button").forEach((x) => x.classList.toggle("on", x === b));
    render();
  }));
  drawWeights();
  $("reset-weights").addEventListener("click", () => { S.weights = { ...S.meta.weights }; drawWeights(); render(); });
  document.querySelectorAll(".tabs button").forEach((b) => b.addEventListener("click", () => {
    S.list = b.dataset.list;
    document.querySelectorAll(".tabs button").forEach((x) => x.classList.toggle("on", x === b));
    renderLists();
  }));
  $("close-panel").addEventListener("click", () => { $("panel").classList.add("hidden"); S.sel = null; });
  $("reverse-btn").addEventListener("click", reverseSearch);
  $("blind-btn").addEventListener("click", () => { renderBlind(); $("blind-panel").classList.remove("hidden"); });
  $("close-blind").addEventListener("click", () => $("blind-panel").classList.add("hidden"));
  $("tour-btn").addEventListener("click", tour);
}

// ---------------------------------------------------------------- blind-predictions panel
function renderBlind() {
  const bp = S.meta.blind_predictions || {};
  const tn = S.meta.targets;
  let html = `<p>Hematite, carbonates, clays and sulfates were <b>never inputs</b> to the global score (ASTER's five TIR bands cannot resolve them). After ranking, we looked at the top candidates with NASA's EMIT imaging spectrometer (60 m) and compared each candidate box with ${"20,000"} random pixels from the same EMIT scenes. Whatever came out is shown, including misses.</p>`;
  for (const [t, rows] of Object.entries(bp)) {
    if (!rows.length) continue;
    html += `<h3>${tn[t].name} — expected: ${rows[0].expected.join(", ")}</h3><table><tr><th>Candidate</th>${rows[0].expected.map((m) => `<th>${m}<br><small>box vs scene baseline</small></th>`).join("")}</tr>`;
    for (const r of rows) {
      html += `<tr><td>#${r.rank} ${r.country} (${r.lat.toFixed(1)}, ${r.lon.toFixed(1)})</td>` + r.expected.map((m) => {
        if (r.emit_status !== "ok") return `<td class="muted">n/a (${r.emit_status})</td>`;
        const a = r.fraction[m] * 100, b = r.baseline[m] * 100;
        const cls = a >= 5 && a > 1.5 * b ? "hit" : "miss";
        return `<td><span class="${cls}">${a.toFixed(1)}%</span> vs ${b.toFixed(1)}%</td>`;
      }).join("") + "</tr>";
    }
    html += "</table>";
  }
  html += `<p class="muted small">Green = mineral present in ≥5% of candidate pixels and enriched ≥1.5× over the scene baseline. Source: <code>src/case_study.py</code>, EMIT L2B MIN v001.</p>`;
  $("blind-content").innerHTML = html;
}

// ---------------------------------------------------------------- story tour
async function tour() {
  if (S.tourTimer) { clearTimeout(S.tourTimer); S.tourTimer = null; $("caption").classList.add("hidden"); $("tour-btn").textContent = "▶ Tour (60 s)"; return; }
  $("tour-btn").textContent = "■ Stop tour";
  const pick = (t, n) => (S.meta.candidates[t] || []).filter((c) => c.list === "discovery")[n];
  const steps = [
    { t: "meridiani", c: pick("meridiani", 0), z: 6, cap: "Meridiani's #1 look-alike: a flat, dust-free, mafic plain. EMIT later found hematite here, Meridiani's signature mineral, without it ever being an input." },
    { t: "jezero", c: pick("jezero", 0), z: 6, cap: "Jezero's best under-studied twin: matching relief, thermal inertia and TIR composition, thousands of km from any classic analog site." },
    { t: "gale", c: pick("gale", 0), z: 6, cap: "Gale's: layered mountain-desert terrain. Round 1 put Gale's twins in Greenland; a polar artifact in apparent thermal inertia was found and fixed." },
    { t: S.target, c: S.meta.reverse_demo, z: 7, cap: "Any place works both ways: click the Jordan Valley and ask where on Mars it fits best.", reverse: true },
  ].filter((s) => s.c);
  let k = 0;
  const next = async () => {
    if (k >= steps.length) { $("caption").classList.add("hidden"); S.tourTimer = null; $("tour-btn").textContent = "▶ Tour (60 s)"; renderBlind(); $("blind-panel").classList.remove("hidden"); return; }
    const s = steps[k++];
    if (s.t !== S.target) { $("target").value = s.t; await loadTarget(s.t); }
    S.map.flyTo([s.c.lat, s.c.lon], s.z, { duration: 3 });
    $("caption").textContent = s.cap; $("caption").classList.remove("hidden");
    setTimeout(() => showCell(s.c.lat, s.c.lon, s.c.name || `Candidate #${s.c.rank} — ${s.c.country}`).then(() => { if (s.reverse) reverseSearch(); }), 3200);
    S.tourTimer = setTimeout(next, 13000);
  };
  next();
}

async function main() {
  S.meta = await (await fetch("data/meta.json")).json();
  S.weights = { ...S.meta.weights };
  S.map = L.map("map", { worldCopyJump: true, minZoom: 2, maxZoom: 9 }).setView([20, 0], 2);
  const esri = "https://server.arcgisonline.com/ArcGIS/rest/services/Canvas";
  L.tileLayer(`${esri}/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}`, {
    attribution: "Basemap &copy; Esri, HERE, Garmin, &copy; OpenStreetMap contributors", maxZoom: 16 }).addTo(S.map);
  S.map.createPane("labels"); S.map.getPane("labels").style.zIndex = 650; S.map.getPane("labels").style.pointerEvents = "none";
  L.tileLayer(`${esri}/World_Dark_Gray_Reference/MapServer/tile/{z}/{y}/{x}`, { pane: "labels", maxZoom: 16 }).addTo(S.map);
  S.map.on("click", (e) => showCell(e.latlng.lat, ((e.latlng.lng + 540) % 360) - 180));
  buildControls();
  await loadTarget(Object.keys(S.meta.targets)[0]);
  ensureEarth();
}

main().catch((e) => { document.body.insertAdjacentHTML("beforeend", `<pre style="position:fixed;bottom:0;left:0;background:#400;color:#fff;padding:8px;z-index:9999">${e}</pre>`); });
