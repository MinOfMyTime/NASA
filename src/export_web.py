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
from similarity import DIMS  # noqa: E402

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
}


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
    feats = {}
    for k, (label, unit, lo, hi) in PANEL_FEATURES.items():
        a = to_web(np.load(PROC / "earth" / f"{k}.npy", mmap_mode="r"), cfg)
        valid = to_web(np.load(PROC / "earth" / "valid_mask.npy", mmap_mode="r").astype(np.float32), cfg)
        a[~(valid >= 0.5)] = np.nan
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
        "validation": val,
        "candidates": cands,
        "year": cfg["year"],
    }
    json.dump(clean(meta), open(WEB / "meta.json", "w"), allow_nan=False,
              separators=(",", ":"))
    tot = sum(p.stat().st_size for p in WEB.iterdir())
    print(f"wrote {len(list(WEB.iterdir()))} files, {tot / 1e6:.1f} MB, grid {shape}")


if __name__ == "__main__":
    main()
