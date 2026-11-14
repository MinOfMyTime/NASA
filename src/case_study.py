"""Zoom-in case studies for candidate analogs: ASTER GED 100 m + EMIT L2B minerals (60 m) + Copernicus GLO-30 DEM.

  python src/case_study.py <target> [n=3]      # top-n under-studied candidates of results/candidates_<target>.csv

Writes results/analog_cards/<target>_<rank>.md and figures. All data are fetched on demand (EMIT/ASTER via the
Earthdata token, GLO-30 from the public AWS bucket); nothing is filled in where coverage is missing.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import download_earth as de  # noqa: E402
from features_earth import mineral_indices  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CARDS = ROOT / "results" / "analog_cards"
EMIT_MATRIX_URL = "https://raw.githubusercontent.com/emit-sds/emit-sds-l2b/develop/data/mineral_grouping_matrix_20230503.csv"
EMIT_10 = ["Calcite", "Chlorite", "Dolomite", "Goethite", "Gypsum", "Hematite", "Illite+Muscovite", "Kaolinite",
           "Montmorillonite", "Vermiculite"]
HALF = 0.15  # half box size, deg (~33 km box)


# ---------------------------------------------------------------- ASTER GED 100 m
def ag100_box(lat, lon, half=HALF):
    import h5py

    de.ag100(lat, lon, half=half)
    out = {}
    stack = []
    for f in sorted((de.RAW / "AG100").glob("*.h5")):
        with h5py.File(f) as h:
            glat = h["Geolocation/Latitude"][()]
            glon = h["Geolocation/Longitude"][()]
            top, left = float(np.round(glat.max())), float(np.round(glon.min()))
            if not (top - 1 <= lat + half and top >= lat - half and left <= lon + half and left + 1 >= lon - half):
                continue
            em = h["Emissivity/Mean"][()].astype(np.float32)
            ndvi = h["NDVI/Mean"][()].astype(np.float32)
            lw = h["Land Water Map/LWmap"][()]
        ny, nx = ndvi.shape
        la = top - (np.arange(ny) + 0.5) / ny
        lo = left + (np.arange(nx) + 0.5) / nx
        rs = (la <= lat + half) & (la >= lat - half)
        cs = (lo >= lon - half) & (lo <= lon + half)
        em = em[:, rs][:, :, cs] * 0.001
        em[(em < 0.5) | (em > 1.0)] = np.nan
        nd = ndvi[rs][:, cs] * 0.01
        nd[ndvi[rs][:, cs] == -9999] = np.nan
        land = lw[rs][:, cs] == 0
        stack.append((em, nd, land))
    if not stack:
        return None
    em = np.concatenate([s[0].reshape(5, -1) for s in stack], axis=1)
    nd = np.concatenate([s[1].ravel() for s in stack])
    land = np.concatenate([s[2].ravel() for s in stack])
    qi, mi, con = mineral_indices(*em)
    m = land & np.isfinite(qi)
    out = {"n_px_100m": int(m.sum()), "ndvi_median": float(np.nanmedian(nd[m])),
           "quartz_index_median": float(np.nanmedian(qi[m])), "mafic_index_median": float(np.nanmedian(mi[m])),
           "tir_contrast_median": float(np.nanmedian(con[m])),
           "bare_fraction": float(np.mean(nd[m] < 0.15))}
    # keep the largest single tile window for the figure
    out["_fig"] = stack[0]
    return out


# ---------------------------------------------------------------- EMIT L2B minerals
def emit_box(lat, lon, half=HALF):
    import netCDF4

    try:
        granules = de.emit(lat, lon, half=half)
    except Exception as e:  # noqa: BLE001
        return {"status": f"EMIT search/download failed: {e}"}
    if not granules:
        return {"status": "no EMIT L2B MIN granules intersect this box"}
    mpath = de.RAW / "emit_mineral_grouping_matrix_20230503.csv"
    if not mpath.exists():
        import urllib.request

        mpath.write_bytes(urllib.request.urlopen(EMIT_MATRIX_URL, timeout=60).read())
    mat = pd.read_csv(mpath).set_index("Index")
    counts = {k: 0 for k in EMIT_10}
    depth = {k: 0.0 for k in EMIT_10}
    n_valid = 0
    used = []
    for g in granules[:4]:
        link = [u for u in g.data_links() if "_MIN_" in u][0]
        f = de.RAW / "EMITL2BMIN" / Path(link).name
        if not f.exists():
            continue
        ds = netCDF4.Dataset(f)
        gt = ds.geotransform  # GLT grid geotransform (lon0, dx, 0, lat0, 0, dy)
        gx = ds["location/glt_x"][:].filled(0).astype(int)
        gy = ds["location/glt_y"][:].filled(0).astype(int)
        lon_g = gt[0] + (np.arange(gx.shape[1]) + 0.5) * gt[1]
        lat_g = gt[3] + (np.arange(gx.shape[0]) + 0.5) * gt[5]
        rs = (lat_g <= lat + half) & (lat_g >= lat - half)
        cs = (lon_g >= lon - half) & (lon_g <= lon + half)
        sx, sy = gx[rs][:, cs], gy[rs][:, cs]
        ok = (sx > 0) & (sy > 0)  # GLT is 1-based; 0 = no swath pixel
        if not ok.any():
            ds.close()
            continue
        used.append(f.name)
        for grp in ("group_1", "group_2"):
            mid = ds[f"{grp}_mineral_id"][:].filled(0)[sy[ok] - 1, sx[ok] - 1].astype(int)
            bd = ds[f"{grp}_band_depth"][:].filled(0)[sy[ok] - 1, sx[ok] - 1].astype(float)
            for k in EMIT_10:
                w = mat[k].reindex(mid).fillna(0).to_numpy()
                hit = (w > 0) & (bd > 0)
                counts[k] += int(hit.sum())
                depth[k] += float((bd * w)[hit].sum())
        n_valid += int(ok.sum())
        ds.close()
    if n_valid == 0:
        return {"status": "EMIT granules found but none cover the box after GLT orthorectification"}
    return {"status": "ok", "granules": used, "n_px_60m": n_valid,
            "fraction_px": {k: counts[k] / n_valid for k in EMIT_10},
            "mean_band_depth_x_weight": {k: depth[k] / n_valid for k in EMIT_10}}


# ---------------------------------------------------------------- Copernicus GLO-30
def glo30_box(lat, lon, half=HALF):
    import os

    os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    import rasterio
    from rasterio.windows import from_bounds

    la, lo = int(np.floor(lat)), int(np.floor(lon))
    ns, ew = ("N" if la >= 0 else "S"), ("E" if lo >= 0 else "W")
    name = f"Copernicus_DSM_COG_10_{ns}{abs(la):02d}_00_{ew}{abs(lo):03d}_00_DEM"
    url = f"/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif"
    try:
        with rasterio.open(url) as d:
            w = from_bounds(lon - half, lat - half, lon + half, lat + half, d.transform)
            z = d.read(1, window=w, boundless=True, fill_value=np.nan).astype(np.float32)
            res_deg = d.res
    except Exception as e:  # noqa: BLE001
        return {"status": f"GLO-30 read failed ({name}): {e}"}
    dy = res_deg[1] * 111.32e3
    dx = res_deg[0] * 111.32e3 * np.cos(np.radians(lat))
    gy, gx = np.gradient(z, dy, dx)
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    return {"status": "ok", "tile": name, "relief_m": float(np.nanmax(z) - np.nanmin(z)),
            "elev_median_m": float(np.nanmedian(z)), "slope30m_median_deg": float(np.nanmedian(slope)),
            "slope30m_p90_deg": float(np.nanpercentile(slope, 90)),
            "frac_slope_gt_15deg": float(np.nanmean(slope > 15)), "_z": z, "_slope": slope}


# ---------------------------------------------------------------- card
MISSION = {
    "topo": "rover mobility and landing-hazard tests on terrain with Mars-like relief and slope statistics",
    "thermo": "thermal-inertia-based landing-site certification, regolith/bedrock discrimination with TIR cameras",
    "mineral": "TIR/VNIR spectrometer field validation (TES/THEMIS/Mini-TES-style instruments) on a matching composition mix",
    "aridity": "contamination-controlled sampling and life-detection limits in near-sterile ground",
}


def card(target, fp, row, val_dims, ag, emit, dem, figpath):
    dims = ["topo", "thermo", "mineral", "aridity", "extremes"]
    labels = {"topo": "Topography", "thermo": "Thermophysics", "mineral": "Mineralogy (TIR)", "aridity": "Aridity",
              "extremes": "Thermal extremes"}
    s = {d: row[f"s_{d}"] for d in dims}
    strong = [labels[d] for d in dims if s[d] >= 0.6]
    weak = [labels[d] for d in dims if s[d] < 0.35]
    L = [f"# Analog Card — {row['country']} ({row['lat']:.2f}, {row['lon']:.2f}) vs {fp['name']}", "",
         f"Under-studied candidate #{int(row['rank'])} for {fp['name']} · region area {row['area_km2']:.0f} km² · "
         f"mean overall score **{row['mean_score']:.2f}** (top 0.5% of Earth's land) · nearest listed analog "
         f"{row['nearest_known_km']:.0f} km away.", "",
         f"![case study]({figpath.name})", "", "## What matches, and how well (0–1)", "",
         "| Dimension | Score | Mars target (median [IQR]) | Earth region |", "|---|---|---|---|"]
    for d in dims:
        L.append(f"| {labels[d]} | {s[d]:.2f} | {val_dims[d][0]} | {val_dims[d][1]} |")
    L += ["", f"**Strong:** {', '.join(strong) or '—'}. **Weak:** {', '.join(weak) or '—'}.", "",
          "## Zoom-in evidence (≈33 km box)", ""]
    if ag:
        L.append(f"- **ASTER GED 100 m** ({ag['n_px_100m']:,} px): bare fraction (NDVI<0.15) {ag['bare_fraction']:.0%}, "
                 f"median NDVI {ag['ndvi_median']:.2f}, quartz index {ag['quartz_index_median']:.3f}, "
                 f"mafic index {ag['mafic_index_median']:.3f}, TIR contrast {ag['tir_contrast_median']:.3f}.")
    else:
        L.append("- **ASTER GED 100 m**: no tiles cover this box.")
    if emit.get("status") == "ok":
        top = [kv for kv in sorted(emit["fraction_px"].items(), key=lambda kv: -kv[1])[:4] if kv[1] >= 0.005]
        found = ", ".join(f"{k} in {v:.0%} of pixels" for k, v in top) or "none of the ten core minerals in >=0.5% of pixels"
        L.append(f"- **EMIT L2B minerals 60 m** ({emit['n_px_60m']:,} px from {len(emit['granules'])} scene(s)): "
                 + found
                 + ". (Pixel fraction where the mineral appears in group-1 or group-2 identification with positive band depth.)")
    else:
        L.append(f"- **EMIT**: {emit.get('status')}.")
    if dem.get("status") == "ok":
        L.append(f"- **Copernicus GLO-30 DEM**: relief {dem['relief_m']:.0f} m, median 30 m slope {dem['slope30m_median_deg']:.1f}°, "
                 f"90th pct {dem['slope30m_p90_deg']:.1f}°, {dem['frac_slope_gt_15deg']:.0%} of the box steeper than 15° (landing-hazard proxy).")
    else:
        L.append(f"- **DEM**: {dem.get('status')}.")
    L += ["", "## Where the analogy breaks", "",
          "- **Gravity** 1 g vs 0.38 g; **pressure** ~1013 mbar vs ~6 mbar (no stable liquid water on Mars, different aeolian thresholds).",
          "- **Radiation**: Earth's magnetic field and atmosphere shield the surface; Mars surface UV/GCR dose is orders of magnitude higher.",
          f"- **Temperature** (thermal-extremes score {s['extremes']:.2f}): Mars target {val_dims['extremes'][0]} vs Earth region "
          f"{val_dims['extremes'][1]} (coldest monthly-mean night LST).",
          ]
    if s["aridity"] < 0.6:
        L.append("- **Biology**: vegetation/soil life present — a drawback for sterility, an asset for biosignature-detection tests.")
    if s["thermo"] < 0.35:
        L.append("- **Thermophysics** mismatch: Earth apparent thermal inertia/albedo differ from the target (moisture, vegetation, or dust cover).")
    L += ["", "## Mission tests it could support", ""]
    for d in ["topo", "thermo", "mineral", "aridity"]:
        if s[d] >= 0.5:
            L.append(f"- {MISSION[d]}.")
    if dem.get("status") == "ok" and dem["frac_slope_gt_15deg"] < 0.05:
        L.append("- landing-zone and EDL-sensor tests (gentle slopes over most of the box).")
    if emit.get("status") == "ok" and any(emit["fraction_px"][k] > 0.05 for k in ["Hematite", "Goethite", "Gypsum", "Kaolinite", "Montmorillonite"]):
        L.append("- ISRU and aqueous-mineral detection (hydrated/oxidised minerals mapped by EMIT).")
    L += ["", "_Generated by `src/case_study.py` from the candidates CSV, ASTER GED AG100, EMIT L2B MIN and Copernicus GLO-30._"]
    return "\n".join(L)


def figure(path, title, ag, dem):
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    if ag and "_fig" in ag:
        em, nd, _ = ag["_fig"]
        qi, mi, _ = mineral_indices(*em)
        axs[0].imshow(qi, cmap="cividis", vmin=np.nanpercentile(qi, 2), vmax=np.nanpercentile(qi, 98))
        axs[0].set_title("ASTER 100 m quartz index")
        axs[1].imshow(mi, cmap="viridis", vmin=np.nanpercentile(mi, 2), vmax=np.nanpercentile(mi, 98))
        axs[1].set_title("ASTER 100 m mafic index")
    if dem.get("status") == "ok":
        axs[2].imshow(dem["_slope"], cmap="magma", vmin=0, vmax=30)
        axs[2].set_title("GLO-30 slope (0–30°)")
    for a in axs:
        a.set_xticks([])
        a.set_yticks([])
    fig.suptitle(title)
    fig.savefig(path, dpi=80, bbox_inches="tight")
    plt.close(fig)


def main(target, n=3):
    de.login()
    CARDS.mkdir(parents=True, exist_ok=True)
    fps = json.load(open(ROOT / "data" / "processed" / "mars_fingerprints.json"))
    fp = fps[target]
    df = pd.read_csv(ROOT / "results" / f"candidates_{target}.csv")
    disc = df[df.list == "discovery"].head(n)
    E = {k: np.load(ROOT / "data" / "processed" / "earth" / f"{k}.npy", mmap_mode="r")
         for k in ["log_relief", "slope_deg", "pct_ti", "albedo", "comp_mafic", "comp_silica", "comp_fines", "ndvi",
                   "diurnal_range_k", "lst_min_k"]}
    from geo import Grid

    g = Grid(20, -180.0)
    summary = []
    for _, row in disc.iterrows():
        r, c = g.rowcol(row.lat, row.lon)
        sl = (slice(r - 2, r + 3), slice(c - 2, c + 3))
        ev = {k: float(np.nanmedian(np.asarray(v[sl]))) for k, v in E.items()}
        f = fp["features"]
        mv = lambda k, d=2: f"{f[k]['median']:.{d}f} [{f[k]['q25']:.{d}f}, {f[k]['q75']:.{d}f}]"  # noqa: E731
        val_dims = {
            "topo": (f"relief {10 ** f['log_relief']['median'] - 1:.0f} m, slope {f['slope_deg']['median']:.1f}°",
                     f"relief {10 ** ev['log_relief'] - 1:.0f} m, slope {ev['slope_deg']:.1f}°"),
            "thermo": (f"TI pct {mv('pct_ti')}, albedo {mv('albedo')}", f"ATI pct {ev['pct_ti']:.2f}, albedo {ev['albedo']:.2f}"),
            "mineral": (f"mafic/silica/fines {f['comp_mafic']['median']:.2f}/{f['comp_silica']['median']:.2f}/{f['comp_fines']['median']:.2f}",
                        f"{ev['comp_mafic']:.2f}/{ev['comp_silica']:.2f}/{ev['comp_fines']:.2f}"),
            "aridity": ("NDVI 0", f"NDVI {ev['ndvi']:.2f}"),
            "extremes": (f"ΔT {f['diurnal_range_k']['median']:.0f} K, min {f['lst_min_k']['median']:.0f} K (model)",
                         f"ΔT {ev['diurnal_range_k']:.0f} K, min {ev['lst_min_k']:.0f} K"),
        }
        print(f"== {target} #{int(row['rank'])} {row.country} ({row.lat:.2f},{row.lon:.2f})")
        ag = ag100_box(row.lat, row.lon)
        emit = emit_box(row.lat, row.lon)
        dem = glo30_box(row.lat, row.lon)
        figp = CARDS / f"{target}_{int(row['rank'])}.png"
        figure(figp, f"{target} candidate #{int(row['rank'])} — {row.country}", ag, dem)
        md = card(target, fp, row, val_dims, ag, emit, dem, figp)
        (CARDS / f"{target}_{int(row['rank'])}.md").write_text(md)
        summary.append({"target": target, "rank": int(row["rank"]), "emit": emit.get("status"),
                        "ag100": bool(ag), "dem": dem.get("status")})
        print("   ", summary[-1])
    return summary


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 3)
