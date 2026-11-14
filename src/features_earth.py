"""Earth features on the 0.05 deg CMG grid (7200 x 3600, -180..180E, north-up).

Inputs (data/raw/earth/): MYD11C3 (12 monthly HDF4), MCD43C3 (daily HDF4 samples), AG1KM (1 deg HDF5 tiles).
Outputs: data/processed/earth/<feature>.npy (float32) + quantile tables earth_quantiles.npz.

Scale factors / fills / QC (verified against the file attributes at read time, see read_sds):
  MYD11C3 LST_*_CMG uint16 * 0.02 K, fill 0. QC_* bits 0-1 mandatory QA (0 good, 1 other quality, 2/3 not
    produced) and bits 6-7 LST error (0: <=1 K, 1: <=2 K, 2: <=3 K, 3: >3 K). Keep QA<=1 and error<=2 K.
  MCD43C3 Albedo_WSA_shortwave int16 * 0.001, fill 32767. Albedo_Quality 0 (best) .. 3 kept (255 fill); Percent_Snow == 0 kept
    (snow-free albedo).
  AG1KM Emissivity/Mean bands 10-14 * 0.001, NDVI/Mean * 0.01, ASTER GDEM/ASTGDEM metres, Land Water Map/LWmap,
    fill -9999.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo import Grid  # noqa: E402
from stats import area_weights, quantile_table, to_percentile  # noqa: E402
from topo import topo_features  # noqa: E402
from features_mars import composition_vector  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "earth"
OUT = ROOT / "data" / "processed" / "earth"
G05 = Grid(20, lon0=-180.0)
G30 = Grid(30, lon0=-180.0)  # 1/30 deg = 3.706 km N-S, matching the MOLA MEGDR pixel (3.705 km)


# ----------------------------------------------------------------------------------------------------------- MODIS
def read_sds(path, name):
    """Read an HDF4 SDS applying its own scale_factor / add_offset / _FillValue attributes."""
    from pyhdf.SD import SD, SDC

    f = SD(str(path), SDC.READ)
    s = f.select(name)
    raw = s.get()
    att = s.attributes()
    f.end()
    return raw, att


def lst_qc_ok(qc: np.ndarray) -> np.ndarray:
    qc = qc.astype(np.uint8)
    mandatory = qc & 0b11
    lst_err = (qc >> 6) & 0b11
    return (mandatory <= 1) & (lst_err <= 1)


def scale_lst(raw: np.ndarray, att: dict) -> np.ndarray:
    fill = att.get("_FillValue", 0)
    sf = float(att.get("scale_factor", 0.02))
    off = float(att.get("add_offset", 0.0))
    out = raw.astype(np.float32) * sf + off  # MODIS LST convention: K = DN * scale_factor (+ 0 offset)
    out[raw == fill] = np.nan
    return out


def read_myd11c3(path):
    day_raw, da = read_sds(path, "LST_Day_CMG")
    night_raw, na = read_sds(path, "LST_Night_CMG")
    qd, _ = read_sds(path, "QC_Day")
    qn, _ = read_sds(path, "QC_Night")
    day = scale_lst(day_raw, da)
    night = scale_lst(night_raw, na)
    day[~lst_qc_ok(qd)] = np.nan
    night[~lst_qc_ok(qn)] = np.nan
    return day, night


def scale_albedo(raw, att):
    fill = att.get("_FillValue", 32767)
    sf = float(att.get("scale_factor", 0.001))
    off = float(att.get("add_offset", 0.0))
    out = (raw.astype(np.float32) - off) * sf if off else raw.astype(np.float32) * sf
    out[raw == fill] = np.nan
    out[(out < 0) | (out > 1)] = np.nan
    return out


def read_mcd43c3(path):
    raw, att = read_sds(path, "Albedo_WSA_shortwave")
    a = scale_albedo(raw, att)
    q, _ = read_sds(path, "Albedo_Quality")  # 0 best .. 3 magnitude inversions; 255 fill
    snow, _ = read_sds(path, "Percent_Snow")
    a[(q > 3) | (snow > 0)] = np.nan
    return a


def apparent_thermal_inertia(albedo, lst_day, lst_night, min_dt=1.0):
    """ATI = (1 - albedo) / (LST_day - LST_night). Cells with dT < min_dt K are undefined (NaN)."""
    dt = lst_day - lst_night
    with np.errstate(invalid="ignore", divide="ignore"):
        ati = (1.0 - albedo) / dt
    ati[~(dt >= min_dt)] = np.nan
    return ati.astype(np.float32)


def modis_features(cfg):
    files = sorted((RAW / "MYD11C3").glob("MYD11C3.A*.hdf"))
    assert len(files) == 12, f"need 12 MYD11C3 files, found {len(files)}"
    days, nights = [], []
    for f in files:
        d, n = read_myd11c3(f)
        days.append(d)
        nights.append(n)
        print(f"  {f.name}: day valid={np.isfinite(d).mean():.3f} [{np.nanmin(d):.1f},{np.nanmax(d):.1f}] K, "
              f"night valid={np.isfinite(n).mean():.3f} [{np.nanmin(n):.1f},{np.nanmax(n):.1f}] K")
    day = np.stack(days)
    night = np.stack(nights)
    alb_files = sorted((RAW / "MCD43C3").glob("MCD43C3.A*.hdf"))
    assert alb_files, "no MCD43C3 files"
    alb_stack = []
    for f in alb_files:
        a = read_mcd43c3(f)
        alb_stack.append(a)
        print(f"  {f.name}: albedo valid={np.isfinite(a).mean():.3f} median={np.nanmedian(a):.3f}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        albedo = np.nanmean(np.stack(alb_stack), axis=0).astype(np.float32)
        ati_m = np.stack([apparent_thermal_inertia(albedo, day[i], night[i]) for i in range(12)])
        feats = {
            "albedo": albedo,
            "ati": np.nanmedian(ati_m, axis=0).astype(np.float32),
            "diurnal_range_k": np.nanmean(day - night, axis=0).astype(np.float32),
            "lst_min_k": np.nanmin(night, axis=0).astype(np.float32),
            "lst_max_k": np.nanmax(day, axis=0).astype(np.float32),
            "lst_day_mean_k": np.nanmean(day, axis=0).astype(np.float32),
            "lst_night_mean_k": np.nanmean(night, axis=0).astype(np.float32),
        }
    return feats


# ----------------------------------------------------------------------------------------------------------- ASTER GED
def read_ag_tile(path):
    import h5py

    with h5py.File(path, "r") as h:
        glat = h["Geolocation/Latitude"][()]
        glon = h["Geolocation/Longitude"][()]
        # Geolocation arrays span the tile edges inclusively (e.g. -12.0..-11.0 over 100 px), so rebuild
        # pixel centres from the 1x1 deg tile extent.
        ny, nx = glat.shape
        top, left = float(np.round(glat.max())), float(np.round(glon.min()))
        lat = np.broadcast_to((top - (np.arange(ny) + 0.5) / ny)[:, None], (ny, nx))
        lon = np.broadcast_to((left + (np.arange(nx) + 0.5) / nx)[None, :], (ny, nx))
        em = h["Emissivity/Mean"][()].astype(np.float32)
        ndvi = h["NDVI/Mean"][()].astype(np.float32)
        dem = h["ASTER GDEM/ASTGDEM"][()].astype(np.float32)
        lw = h["Land Water Map/LWmap"][()]
    em[em == -9999] = np.nan
    em *= 0.001
    ndvi[ndvi == -9999] = np.nan
    ndvi *= 0.01
    dem[dem == -9999] = np.nan
    return lat, lon, em, ndvi, dem, lw


def lw_is_land(lw):
    """ASTER GED v3 land-water map: 0 = land, 1 = water.
    Verified empirically in scripts/inspect_earth.py against GDEM elevation (ocean pixels have elevation 0)."""
    return lw == LAND_CODE


LAND_CODE = 0


class Accum:
    def __init__(self, grid: Grid, names):
        self.g = grid
        n = grid.nrow * grid.ncol
        self.sum = {k: np.zeros(n, np.float64) for k in names}
        self.cnt = {k: np.zeros(n, np.float64) for k in names}

    def add(self, lat, lon, name, val):
        r, c = self.g.rowcol(lat, lon)
        idx = (r * self.g.ncol + c).ravel()
        v = val.ravel()
        ok = np.isfinite(v)
        n = self.g.nrow * self.g.ncol
        self.sum[name] += np.bincount(idx[ok], weights=v[ok], minlength=n)
        self.cnt[name] += np.bincount(idx[ok], minlength=n)

    def mean(self, name):
        with np.errstate(invalid="ignore", divide="ignore"):
            return (self.sum[name] / self.cnt[name]).reshape(self.g.shape).astype(np.float32)


def aster_features(cfg):
    files = sorted((RAW / "AG1KM").glob("*.h5"))
    print(f"  AG1KM tiles: {len(files)}")
    acc5 = Accum(G05, ["b10", "b11", "b12", "b13", "b14", "ndvi", "land", "dem5"])
    acc30 = Accum(G30, ["dem"])
    for i, f in enumerate(files):
        try:
            lat, lon, em, ndvi, dem, lw = read_ag_tile(f)
        except Exception as e:  # noqa: BLE001
            print(f"  WARN unreadable tile {f.name}: {e}")
            continue
        land = lw_is_land(lw)
        for b in range(5):
            e = np.where(land, em[b], np.nan)
            acc5.add(lat, lon, f"b1{b}", e)
        acc5.add(lat, lon, "ndvi", np.where(land, ndvi, np.nan))
        acc5.add(lat, lon, "land", land.astype(np.float32))
        acc5.add(lat, lon, "dem5", dem)
        acc30.add(lat, lon, "dem", dem)
        if i % 2000 == 0:
            print(f"   {i}/{len(files)}", flush=True)
    em = {f"emis_b1{b}": acc5.mean(f"b1{b}") for b in range(5)}
    land_frac = acc5.mean("land")
    return em, acc5.mean("ndvi"), land_frac, acc30.mean("dem"), acc5.mean("dem5")


def mineral_indices(b10, b11, b12, b13, b14):
    """ASTER TIR indices (Ninomiya, Fu & Cudahy 2005, Remote Sens. Environ. 99, 127-139):
    Quartz index QI = B11^2 / (B10 * B12); Mafic index MI = B12 * B14^3 / B13^4.
    Spectral contrast = max - min emissivity of B10-B14 (fine particulates / dust reduce TIR contrast)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        qi = b11 * b11 / (b10 * b12)
        mi = b12 * b14**3 / b13**4
    stack = np.stack([b10, b11, b12, b13, b14])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        contrast = np.nanmax(stack, 0) - np.nanmin(stack, 0)
    contrast[~np.all(np.isfinite(stack), 0)] = np.nan
    return qi.astype(np.float32), mi.astype(np.float32), contrast.astype(np.float32)


# ----------------------------------------------------------------------------------------------------------- main
def main():
    cfg = yaml.safe_load(open(ROOT / "config.yaml"))
    OUT.mkdir(parents=True, exist_ok=True)
    print("MODIS ...")
    feats = modis_features(cfg)
    print("ASTER GED ...")
    em, ndvi, land_frac, dem30, dem5 = aster_features(cfg)
    feats.update(em)
    feats["ndvi"] = ndvi
    feats["land_frac"] = land_frac
    feats["elevation_m"] = dem5

    print("Topography at matched scale (1/30 deg = 3.706 km) ...")
    t = cfg["topography"]
    z = np.where(np.isfinite(dem30), dem30, 0.0).astype(np.float32)  # no tile = ocean -> sea level
    topo30 = topo_features(z, G30.lats, t["pixel_km"], t["window_px"], t["tpi_window_px"])
    r, _ = G30.rowcol(G05.lats, np.zeros(G05.nrow))
    _, c = G30.rowcol(np.zeros(G05.ncol), G05.lons)
    for k, v in topo30.items():
        feats[k] = v[np.ix_(r, c)]

    # masks
    lat2d = np.broadcast_to(G05.lats[:, None], G05.shape)
    land = land_frac >= 0.5
    ice = ~(feats["lst_max_k"] > cfg["masks"]["ice_lst_max_k"])  # warmest month never above 0 C (or no LST)
    has = np.isfinite(feats["emis_b10"]) & np.isfinite(feats["lst_max_k"]) & np.isfinite(feats["ati"])
    valid = land & ~ice & has & (lat2d >= -60) & (lat2d <= 84)
    feats["land_mask"] = land
    feats["ice_mask"] = ice & land
    feats["valid_mask"] = valid

    w = area_weights(G05.lats, G05.shape)
    wv = np.where(valid, w, 0.0)
    q = {"ati": quantile_table(feats["ati"], wv)}
    feats["pct_ti"] = to_percentile(feats["ati"], q["ati"])

    qi, mi, contrast = mineral_indices(*(feats[f"emis_b1{b}"] for b in range(5)))
    bare = valid & (ndvi < cfg["masks"]["ndvi_max_for_mineral"])
    wb = np.where(bare, w, 0.0)
    for k, v in {"quartz_index": qi, "mafic_index": mi, "tir_contrast": contrast}.items():
        feats[k] = v
        q[k] = quantile_table(v, wb)
    pm = to_percentile(mi, q["mafic_index"])
    ps = to_percentile(qi, q["quartz_index"])
    pf = 1.0 - to_percentile(contrast, q["tir_contrast"])
    cm, cs, cf = composition_vector(pm, ps, pf)
    for k, v in {"comp_mafic": cm, "comp_silica": cs, "comp_fines": cf}.items():
        v = v.astype(np.float32)
        v[~bare] = np.nan  # TIR mineral signal is not interpretable under vegetation
        feats[k] = v
    feats["mineral_bare_mask"] = bare

    for k, v in feats.items():
        np.save(OUT / f"{k}.npy", np.asarray(v))
    np.savez_compressed(ROOT / "data" / "processed" / "earth_quantiles.npz", **q)
    print(f"valid land cells: {valid.sum():,} ({(valid * w).sum() / (land * w).sum():.1%} of land area)")
    for k in ["albedo", "ati", "pct_ti", "diurnal_range_k", "lst_min_k", "lst_max_k", "ndvi", "log_relief", "slope_deg",
              "quartz_index", "mafic_index", "tir_contrast", "comp_mafic", "comp_silica", "comp_fines"]:
        v = feats[k][valid]
        print(f"  {k:16s} p5={np.nanpercentile(v, 5):.4g} p50={np.nanmedian(v):.4g} p95={np.nanpercentile(v, 95):.4g}")


if __name__ == "__main__":
    main()
