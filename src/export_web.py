"""Export compact web assets to web/data/.

  <target>_<dim>.bin   uint8 0.25 deg grid (rows = web_lat_range north->south, 1440 cols from -180E);
                        value = round(score*254), 255 = no data
  earth_<feat>.bin     uint8 0.25 deg Earth feature values for the side panel (linear scale/offset in meta.json)
  meta.json            grid geometry, targets + Mars fingerprints, default weights, known-analog validation,
                        candidates, feature encodings, provenance
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo import block_reduce_nan  # noqa: E402
from similarity import DIMS, SCALE_KEY  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
WEB = ROOT / "web" / "data"

PANEL_FEATURES = {  # name: (label, unit, lo, hi)
    "log_relief": ("Local relief (log10 1+m, 18.5 km window)", "log10 m", 0.0, 4.0),
    "slope_deg": ("Mean slope (3.7 km baseline)", "deg", 0.0, 30.0),
    "pct_ti": ("Thermal inertia percentile (Earth: apparent TI)", "0-1", 0.0, 1.0),
    "albedo": ("Broadband albedo", "", 0.0, 0.8),
    "comp_mafic": ("Composition: mafic share", "", 0.0, 1.0),
    "comp_silica": ("Composition: silica-rich share", "", 0.0, 1.0),
    "comp_fines": ("Composition: dust/fines share", "", 0.0, 1.0),
    "ndvi": ("NDVI (ASTER GED mean)", "", -0.2, 1.0),
    "diurnal_range_k": ("Diurnal LST range (day-night)", "K", 0.0, 120.0),
    "lst_min_k": ("Annual min LST (monthly night)", "K", 180.0, 320.0),
    "lst_max_k": ("Annual max LST (monthly day)", "K", 230.0, 345.0),
    "log_rough": ("Roughness (log10 1+m, detrended, 18.5 km window)", "log10 m", 0.0, 4.0),
    "tpi_m": ("Topographic position index (33 km)", "m", -1500.0, 1500.0),
    "extremes_rel": ("Thermal extremes, relative to Earth land (percentile)", "0-1", 0.0, 1.0),
}


MARS_PANEL = ["log_relief", "slope_deg", "log_rough", "tpi_m", "pct_ti", "albedo", "comp_mafic", "comp_silica", "comp_fines", "ndvi"]
MARS_RANGES = {"log_rough": (0.0, 4.0), "tpi_m": (-1500.0, 1500.0), "elevation_m": (-9000.0, 22000.0)}

MISSION_MODES = {  # preset weights + candidate filters (applied client-side to the 0.25 deg Earth features)
    "default": {"label": "Default (balanced)", "weights": None, "filters": {}},
    "landing": {"label": "Landing / EDL test", "weights": {"topo": 0.4, "thermo": 0.35, "mineral": 0.1, "aridity": 0.1, "extremes": 0.05},
                "filters": {"slope_deg_max": 2.0}, "note": "mean 3.7 km slope < 2 deg; GLO-30 slope check (<5% of box > 15 deg) where an Analog Card exists"},
    "mobility": {"label": "Rover mobility", "weights": {"topo": 0.4, "thermo": 0.4, "mineral": 0.05, "aridity": 0.1, "extremes": 0.05}, "filters": {}},
    "isru": {"label": "ISRU / aqueous minerals", "weights": {"topo": 0.1, "thermo": 0.15, "mineral": 0.6, "aridity": 0.1, "extremes": 0.05},
             "filters": {}, "note": "EMIT hydrated/oxidised mineral fractions are shown on the Analog Cards where EMIT coverage exists"},
    "biosig": {"label": "Biosignature instrument tests", "weights": {"topo": 0.2, "thermo": 0.25, "mineral": 0.35, "aridity": 0.15, "extremes": 0.05},
               "filters": {"ndvi_min": 0.05, "ndvi_max": 0.3}, "note": "aridity preference inverted: sparse life present (NDVI 0.05-0.3)"},
    "spectro": {"label": "Spectrometer calibration", "weights": {"topo": 0.05, "thermo": 0.2, "mineral": 0.65, "aridity": 0.1, "extremes": 0.0}, "filters": {}},
}


def hillshade(z, ppd, radius_km, az=315.0, alt=45.0):
    """Simple hillshade of a lat/lon DEM (pixel size from the planet radius)."""
    lats = 90.0 - (np.arange(z.shape[0]) + 0.5) / ppd
    dy = 2 * np.pi * radius_km * 1000 / 360 / ppd
    dx = dy * np.clip(np.cos(np.radians(lats)), 0.05, None)[:, None]
    zz = np.nan_to_num(z, nan=float(np.nanmedian(z)))
    gy, gx = np.gradient(zz, dy, axis=0), np.gradient(zz, axis=1) / dx
    slope = np.arctan(2.0 * np.hypot(gx, gy))  # x2 vertical exaggeration
    aspect = np.arctan2(-gx, gy)
    a, b = np.radians(az), np.radians(alt)
    hs = np.sin(b) * np.cos(slope) + np.cos(b) * np.sin(slope) * np.cos(a - aspect)
    return np.clip(hs, 0, 1)


def write_png(path, img):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.imsave(path, img, cmap="gray", vmin=0, vmax=1)


def blind_predictions(fps):
    p = ROOT / "results" / "analog_cards" / "summary.json"
    rows = json.load(open(p)) if p.exists() else []
    return {t: [r for r in rows if r["target"] == t] for t in fps}


def export_mars(cfg, fps):
    from mars_features_named import MARS_FEATURES

    z = np.load(PROC / "mars_features_4ppd.npz")
    out = {}
    for k in MARS_PANEL:
        lo, hi = MARS_RANGES.get(k) or (PANEL_FEATURES[k][2], PANEL_FEATURES[k][3])
        enc(z[k].astype(np.float32), lo, hi).tofile(WEB / f"mars_{k}.bin")
        out[k] = {"lo": lo, "hi": hi}
    write_png(WEB / "mars_hillshade.png", hillshade(z["elevation_m"], 4, 3396.0))
    return {"grid": {"ppd": 4, "lat_top": 90.0, "lon_left": 0.0, "rows": 720, "cols": 1440, "nodata": 255},
            "features": out, "named": [{"name": n, "lat": la, "lon_e360": lo} for n, la, lo in MARS_FEATURES],
            "reverse_scale": {k: cfg["similarity"]["min_scale"][SCALE_KEY[k]] for k in MARS_PANEL if k in SCALE_KEY}}


def to_web(a, cfg):
    ppd = cfg["grids"]["earth_ppd"]
    f = ppd // cfg["grids"]["web_ppd"]
    lat_lo, lat_hi = cfg["grids"]["web_lat_range"]
    r0 = int((90 - lat_hi) * ppd)
    r1 = int((90 - lat_lo) * ppd)
    return block_reduce_nan(np.asarray(a[r0:r1], dtype=np.float32), f)


def enc(a, lo, hi):
    q = np.round((np.clip(a, lo, hi) - lo) / (hi - lo) * 254)
    q[~np.isfinite(a)] = 255
    return q.astype(np.uint8)


def clean(o):
    """JSON-safe: NaN/inf -> null (allow_nan=False would otherwise raise)."""
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    return o


def main():
    cfg = yaml.safe_load(open(ROOT / "config.yaml"))
    fps = json.load(open(PROC / "mars_fingerprints.json"))
    val = json.load(open(PROC / "validation.json"))
    WEB.mkdir(parents=True, exist_ok=True)
    shape = None
    for t in fps:
        for d in DIMS:
            w = to_web(np.load(PROC / "similarity" / f"{t}_{d}.npy", mmap_mode="r"), cfg)
            shape = w.shape
            enc(w, 0.0, 1.0).tofile(WEB / f"{t}_{d}.bin")
        gap = to_web(np.load(PROC / "similarity" / f"{t}_gap.npy", mmap_mode="r"), cfg)
        enc(gap, 0.0, 1.0).tofile(WEB / f"{t}_gap.bin")
    feats = {}
    valid05 = np.load(PROC / "earth" / "valid_mask.npy").astype(bool)
    for k, (label, unit, lo, hi) in PANEL_FEATURES.items():
        a05 = np.array(np.load(PROC / "earth" / f"{k}.npy", mmap_mode="r"), dtype=np.float32)
        a05[~valid05] = np.nan  # average only valid 0.05 deg cells, like the score grids
        a = to_web(a05, cfg)
        enc(a, lo, hi).tofile(WEB / f"earth_{k}.bin")
        feats[k] = {"label": label, "unit": unit, "lo": lo, "hi": hi}
    cands = {}
    for t in fps:
        df = pd.read_csv(ROOT / "results" / f"candidates_{t}.csv")
        cands[t] = json.loads(df.to_json(orient="records"))
    meta = {
        "grid": {"ppd": cfg["grids"]["web_ppd"], "lat_top": cfg["grids"]["web_lat_range"][1],
                 "lon_left": -180.0, "rows": shape[0], "cols": shape[1], "nodata": 255},
        "dims": DIMS,
        "dim_labels": {"topo": "Topography", "thermo": "Thermophysics", "mineral": "Mineralogy (TIR)",
                       "aridity": "Aridity (bare ground)", "extremes": "Thermal extremes"},
        "weights": cfg["weights"],
        "targets": fps,
        "earth_features": feats,
        "validation": val["development"],
        "blind_validation": val["blind"],
        "null_test": val["null"],
        "blind_refs": {k: v["ref"] for k, v in cfg["blind_analogs"].items()},
        "candidates": cands,
        "year": cfg["year"],
        "mars": export_mars(cfg, fps),
        "mission_modes": MISSION_MODES,
        "blind_predictions": blind_predictions(fps),
        "lat_hist": {tag: json.load(open(ROOT / "results" / f"lat_hist_{tag}.json")) for tag in ("before", "after")
                     if (ROOT / "results" / f"lat_hist_{tag}.json").exists()},
        "reverse_demo": {"name": "Jordan Valley / Dead Sea", "lat": 31.5, "lon": 35.5},
        "mars_constants": {"gravity_g": 0.38, "pressure_mbar": 6, "earth_pressure_mbar": 1013},
    }
    json.dump(clean(meta), open(WEB / "meta.json", "w"), allow_nan=False,
              separators=(",", ":"))
    tot = sum(p.stat().st_size for p in WEB.iterdir())
    print(f"wrote {len(list(WEB.iterdir()))} files, {tot / 1e6:.1f} MB, grid {shape}")


if __name__ == "__main__":
    main()
