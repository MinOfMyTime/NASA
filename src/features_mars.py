"""Mars features on the MOLA MEGDR grid (16 ppd, east-positive 0..360) + target fingerprints.

Outputs
  data/processed/mars_quantiles.npz   area-weighted percentile tables for Mars-relative features
  data/processed/mars_fingerprints.json   per-target distribution (median, IQR, n) of every comparison feature
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import io_mars as im  # noqa: E402
from geo import Grid, resample_nearest  # noqa: E402
from stats import area_weights, quantile_table, robust_summary, to_percentile  # noqa: E402
from topo import topo_features  # noqa: E402
from mars_thermal_model import diurnal_cycle, extremes  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"


def composition_vector(p_mafic, p_silica, p_fines):
    """Shared 3-class composition vector from planet-relative percentiles, normalised to sum 1."""
    s = p_mafic + p_silica + p_fines
    with np.errstate(invalid="ignore", divide="ignore"):
        return p_mafic / s, p_silica / s, p_fines / s


def build_mars(cfg) -> dict:
    tcfg = cfg["topography"]
    mola, gm = im.read_mola_megdr()
    feats = topo_features(mola, gm.lats, tcfg["pixel_km"], tcfg["window_px"], tcfg["tpi_window_px"])
    feats["elevation_m"] = mola

    alb, ga = im.read_tes_albedo()
    feats["albedo"] = resample_nearest(alb, ga, gm)
    ti, gt = im.read_tes_ti(drop_interpolated=True)
    feats["thermal_inertia"] = resample_nearest(ti, gt, gm)

    ori = cfg["mars"]["tes_mineral_orientation"]
    mins = {}
    for nm in im.TES_MINERALS:
        a, g = im.read_tes_mineral(nm, orientation=ori)
        mins[nm] = resample_nearest(a, g, gm)
    feats["tes_mafic"] = mins["plagioclase"] + mins["high_ca_pyroxene"] + mins["low_ca_pyroxene"] + mins["olivine"]
    feats["tes_silica"] = mins["quartz"] + mins["sheet_silicates_glass"]
    feats["tes_dust"] = mins["surface_dust"]
    feats["tes_hematite"] = mins["hematite"]

    # Planet-relative percentiles (area-weighted). TI is in physical units on Mars but only an *apparent* inertia
    # on Earth, so both are compared as within-planet percentiles.
    w = area_weights(gm.lats, gm.shape)
    q = {}
    for k in ("thermal_inertia", "tes_mafic", "tes_silica", "tes_dust"):
        q[k] = quantile_table(feats[k], w)
    feats["pct_ti"] = to_percentile(feats["thermal_inertia"], q["thermal_inertia"])
    pm = to_percentile(feats["tes_mafic"], q["tes_mafic"])
    ps = to_percentile(feats["tes_silica"], q["tes_silica"])
    pf = to_percentile(feats["tes_dust"], q["tes_dust"])
    feats["comp_mafic"], feats["comp_silica"], feats["comp_fines"] = composition_vector(pm, ps, pf)
    feats["ndvi"] = np.zeros_like(mola)  # no vegetation on Mars, by definition
    np.savez_compressed(PROC / "mars_quantiles.npz", **q)
    return feats, gm


COMPARE = {
    "topo": ["log_relief", "slope_deg", "log_rough", "tpi_m"],
    "thermo": ["pct_ti", "albedo"],
    "mineral": ["comp_mafic", "comp_silica", "comp_fines"],
    "aridity": ["ndvi"],
    "extremes": ["diurnal_range_k", "lst_min_k", "lst_max_k"],
}
DISPLAY_EXTRA = ["elevation_m", "thermal_inertia", "tes_mafic", "tes_silica", "tes_dust", "tes_hematite"]


def fingerprints(cfg, feats, gm) -> dict:
    out = {}
    size = cfg["mars"]["window_deg"]
    for key, t in cfg["mars"]["targets"].items():
        rows, cols = gm.window(t["lat"], t["lon_e360"], size)
        fp = {"name": t["name"], "lat": t["lat"], "lon_e360": t["lon_e360"], "window_deg": size, "features": {}}
        for k in [f for fs in COMPARE.values() for f in fs] + DISPLAY_EXTRA:
            if k in feats:
                fp["features"][k] = robust_summary(feats[k][np.ix_(rows, cols)])
        # Thermal extremes: no global orbital surface-temperature map was supplied, so the target window median TES
        # TI and albedo drive a 1-D thermal model (src/mars_thermal_model.py). Model-derived, single value per target.
        ti_med = fp["features"]["thermal_inertia"]["median"]
        alb_med = fp["features"]["albedo"]["median"]
        ext = extremes(ti_med, alb_med, t["lat"])
        for k in COMPARE["extremes"]:
            v = ext[k]
            fp["features"][k] = {"n": 1, "median": v, "q25": v, "q75": v, "iqr": 0.0,
                                 "source": "1-D thermal model (TES TI + albedo median, latitude)"}
        out[key] = fp
    return out


def main():
    cfg = yaml.safe_load(open(ROOT / "config.yaml"))
    PROC.mkdir(parents=True, exist_ok=True)
    feats, gm = build_mars(cfg)
    fps = fingerprints(cfg, feats, gm)
    json.dump(fps, open(PROC / "mars_fingerprints.json", "w"), indent=1)
    for k, fp in fps.items():
        print(f"== {fp['name']}")
        for f, s in fp["features"].items():
            if s["n"]:
                print(f"   {f:18s} n={s['n']:5d} median={s['median']:9.4g} IQR=[{s['q25']:.4g}, {s['q75']:.4g}]")
            else:
                print(f"   {f:18s} NO DATA")
    # small Mars grids for the web panel / quicklook (0.25 deg)
    g4 = Grid(4, lon0=0.0)
    np.savez_compressed(PROC / "mars_features_4ppd.npz",
                        **{k: resample_nearest(feats[k], gm, g4).astype(np.float32)
                           for k in ["elevation_m", "albedo", "thermal_inertia", "comp_mafic", "comp_silica", "comp_fines"]})


if __name__ == "__main__":
    main()
