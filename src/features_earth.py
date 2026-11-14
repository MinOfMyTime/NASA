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


def apparent_thermal_inertia(albedo, lst_day, lst_night, min_dt=10.0, min_day_k=273.15):
    """ATI = (1 - albedo) / (LST_day - LST_night) for one month.

    A month is usable only when dT >= min_dt (default 10 K) and LST_day > min_day_k: ATI diverges as dT -> 0, which
    happens in polar day/night and with near-coincident overpasses, and would otherwise mimic bedrock-like thermal
    inertia at high latitudes (the artifact documented in REPORT.md).
    """
    dt = lst_day - lst_night
    with np.errstate(invalid="ignore", divide="ignore"):
        ati = (1.0 - albedo) / dt
    ati[~((dt >= min_dt) & (lst_day > min_day_k))] = np.nan
    return ati.astype(np.float32)


def annual_ati(albedo, day, night, min_dt=10.0, min_day_k=273.15, min_months=4):
    """Median monthly ATI over usable months; NaN where fewer than min_months months are usable."""
    ati_m = np.stack([apparent_thermal_inertia(albedo, day[i], night[i], min_dt, min_day_k) for i in range(day.shape[0])])
    n_ok = np.isfinite(ati_m).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ati = np.nanmedian(ati_m, axis=0).astype(np.float32)
    ati[n_ok < min_months] = np.nan
    return ati, n_ok.astype(np.uint8)


def read_mod10cm(path):
    """MOD10CM Snow_Cover_Monthly_CMG: 0-100 = % snow; >100 are flags (107 lake ice, 111 night, 237 inland water,
    239 ocean, 250 cloud, 253 no decision, 254 saturated, 255 fill) -> NaN."""
    raw, _ = read_sds(path, "Snow_Cover_Monthly_CMG")
    s = raw.astype(np.float32)
    s[raw > 100] = np.nan
    return s


def permanent_snow(stack, min_frac=50.0, min_months=3):
    """True where every usable month has >= min_frac % snow and at least min_months months are usable."""
    n_ok = np.isfinite(stack).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mn = np.nanmin(stack, axis=0)
    return (n_ok >= min_months) & (mn >= min_frac)


def snow_features(cfg):
    files = sorted((RAW / "MOD10CM").glob("MOD10CM.A*.hdf"))
    assert len(files) == 12, f"need 12 MOD10CM files, found {len(files)}"
    stack = np.stack([read_mod10cm(f) for f in files])
    s = cfg["masks"]["snow"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out = {"snow_min_pct": np.nanmin(stack, axis=0).astype(np.float32),
               "snow_mean_pct": np.nanmean(stack, axis=0).astype(np.float32)}
    out["permanent_snow"] = permanent_snow(stack, s["min_frac_pct"], s["min_months"])
    print(f"  MOD10CM: permanent snow/ice cells = {int(out['permanent_snow'].sum()):,}")
    return out


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
        a = cfg["masks"]["ati"]
        ati, n_ok = annual_ati(albedo, day, night, a["min_dt_k"], a["min_day_k"], a["min_months"])
        feats = {
            "albedo": albedo,
            "ati": ati,
            "ati_n_months": n_ok,
            "diurnal_range_k": np.nanmean(day - night, axis=0).astype(np.float32),
            "lst_min_k": np.nanmin(night, axis=0).astype(np.float32),
            "lst_max_k": np.nanmax(day, axis=0).astype(np.float32),
            "lst_day_mean_k": np.nanmean(day, axis=0).astype(np.float32),
            "lst_night_mean_k": np.nanmean(night, axis=0).astype(np.float32),
        }
    feats.update(snow_features(cfg))
    return feats


# ----------------------------------------------------------------------------------------------------------- ASTER GED
def read_ag_tile(path):
    import h5py

    with h5py.File(path, "r") as h:
        glat = h["Geolocation/Latitude"][()]
        glon = h["Geolocation/Longitude"][()]
        # Geolocation arrays span the tile edges inclusively (e.g. -12.0..-11.0 over 100 px); only the extent is used.
        top, left = float(np.round(glat.max())), float(np.round(glon.min()))
        em = h["Emissivity/Mean"][()].astype(np.float32)
        ndvi = h["NDVI/Mean"][()].astype(np.float32)
        dem = h["ASTER GDEM/ASTGDEM"][()].astype(np.float32)
        lw = h["Land Water Map/LWmap"][()]
        nobs = h["Observations/NumObs"][()]
    em[em == -9999] = np.nan
    em *= 0.001
    em[(em < 0.5) | (em > 1.0)] = np.nan  # physically impossible emissivities occur in ~12% of tiles
    ndvi[ndvi == -9999] = np.nan
    ndvi *= 0.01
    dem[dem == -9999] = np.nan
    return top, left, em, ndvi, dem, lw, nobs


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

    def add_tile(self, top, left, name, val):
        """Area-average a 1x1 deg tile (ny x nx pixels, north-up) into this grid. Only touches the tile footprint."""
        ny, nx = val.shape
        ppd = self.g.ppd
        m = int(round(ppd))
        lr = np.minimum((np.arange(ny) + 0.5) / ny * m, m - 1).astype(int)  # local row of each pixel
        lc = np.minimum((np.arange(nx) + 0.5) / nx * m, m - 1).astype(int)
        loc = (lr[:, None] * m + lc[None, :]).ravel()
        v = val.ravel()
        ok = np.isfinite(v)
        s = np.bincount(loc[ok], weights=v[ok], minlength=m * m).reshape(m, m)
        c = np.bincount(loc[ok], minlength=m * m).reshape(m, m)
        r0 = int(round((90.0 - top) * ppd))
        c0 = int(round((left - self.g.lon0) * ppd)) % self.g.ncol
        S = self.sum[name].reshape(self.g.shape)
        C = self.cnt[name].reshape(self.g.shape)
        S[r0:r0 + m, c0:c0 + m] += s
        C[r0:r0 + m, c0:c0 + m] += c

    def mean(self, name):
        with np.errstate(invalid="ignore", divide="ignore"):
            return (self.sum[name] / self.cnt[name]).reshape(self.g.shape).astype(np.float32)


def aster_features(cfg):
    files = sorted((RAW / "AG1KM").glob("*.h5"))
    print(f"  AG1KM tiles: {len(files)}")
    acc5 = Accum(G05, ["b10", "b11", "b12", "b13", "b14", "ndvi", "land", "nobs"])
    for i, f in enumerate(files):
        try:
            top, left, em, ndvi, dem, lw, nobs = read_ag_tile(f)
        except Exception as e:  # noqa: BLE001
            print(f"  WARN unreadable tile {f.name}: {e}")
            continue
        land = lw_is_land(lw)
        for b in range(5):
            e = np.where(land, em[b], np.nan)
            acc5.add_tile(top, left, f"b1{b}", e)
        acc5.add_tile(top, left, "ndvi", np.where(land, ndvi, np.nan))
        acc5.add_tile(top, left, "land", land.astype(np.float32))
        acc5.add_tile(top, left, "nobs", np.where(land, nobs.astype(np.float32), np.nan))
        if i % 2000 == 0:
            print(f"   {i}/{len(files)}", flush=True)
    em = {f"emis_b1{b}": acc5.mean(f"b1{b}") for b in range(5)}
    land_frac = acc5.mean("land")
    return em, acc5.mean("ndvi"), land_frac, acc5.mean("nobs")


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
    cache = ROOT / "data" / "processed" / "aster_agg.npz"
    if cache.exists() and cache.stat().st_mtime > Path(__file__).stat().st_mtime:
        z = np.load(cache)
        em = {f"emis_b1{b}": z[f"emis_b1{b}"] for b in range(5)}
        ndvi, land_frac, nobs = z["ndvi"], z["land_frac"], z["nobs"]
    else:
        em, ndvi, land_frac, nobs = aster_features(cfg)
        np.savez(cache, ndvi=ndvi, land_frac=land_frac, nobs=nobs, **em)
    # Topography DEM: Copernicus GLO-90 averaged to 1/30 deg (src/download_dem.py); no tile = ocean = sea level
    dem30 = np.load(ROOT / "data" / "processed" / "earth_dem_30ppd.npy")
    feats.update(em)
    feats["ndvi"] = ndvi
    feats["land_frac"] = land_frac
    feats["aster_nobs"] = nobs
    r_, _ = G30.rowcol(G05.lats, np.zeros(G05.nrow))
    _, c_ = G30.rowcol(np.zeros(G05.ncol), G05.lons)
    feats["elevation_m"] = dem30[np.ix_(r_, c_)]

    print("Topography at matched scale (1/30 deg = 3.706 km) ...")
    t = cfg["topography"]
    z = np.where(np.isfinite(dem30), dem30, 0.0).astype(np.float32)  # no Copernicus tile = open water -> sea level
    topo30 = topo_features(z, G30.lats, t["pixel_km"], t["window_px"], t["tpi_window_px"])
    r, _ = G30.rowcol(G05.lats, np.zeros(G05.nrow))
    _, c = G30.rowcol(np.zeros(G05.ncol), G05.lons)
    for k, v in topo30.items():
        feats[k] = v[np.ix_(r, c)]

    # masks
    lat2d = np.broadcast_to(G05.lats[:, None], G05.shape)
    land = land_frac >= 0.5
    # warmest month below 0 C, or MOD10CM permanent snow/ice
    ice = ~(feats["lst_max_k"] > cfg["masks"]["ice_lst_max_k"]) | feats["permanent_snow"]
    has = np.isfinite(feats["emis_b10"]) & np.isfinite(feats["lst_max_k"]) & np.isfinite(feats["ati"])
    valid = land & ~ice & has & (lat2d >= -60) & (lat2d <= 84)
    feats["land_mask"] = land
    feats["ice_mask"] = ice & land
    feats["valid_mask"] = valid

    w = area_weights(G05.lats, G05.shape)
    bare = valid & (ndvi < cfg["masks"]["ndvi_max_for_mineral"])
    wb = np.where(bare, w, 0.0)
    # Mars TI percentiles describe an all-bare planet, so the Earth reference population is bare ground too.
    q = {"ati": quantile_table(feats["ati"], wb)}
    feats["pct_ti"] = to_percentile(feats["ati"], q["ati"])

    qi, mi, contrast = mineral_indices(*(feats[f"emis_b1{b}"] for b in range(5)))
    for k, v in {"quartz_index": qi, "mafic_index": mi, "tir_contrast": contrast}.items():
        feats[k] = v
        q[k] = quantile_table(v, wb)
    pm = to_percentile(mi, q["mafic_index"])
    ps = to_percentile(qi, q["quartz_index"])
    pf = 1.0 - to_percentile(contrast, q["tir_contrast"])
    cm, cs, cf = composition_vector(pm, ps, pf)
    m = cfg["masks"]
    # Near-flat TIR spectra (snow/frost/ice-contaminated pixels) and poorly observed cells carry no usable mineral signal.
    mineral_ok = bare & (contrast >= m["tir_contrast_min"]) & (nobs >= m["aster_nobs_min"])
    for k, v in {"comp_mafic": cm, "comp_silica": cs, "comp_fines": cf}.items():
        v = v.astype(np.float32)
        v[~mineral_ok] = np.nan
        feats[k] = v
    feats["mineral_bare_mask"] = bare
    feats["mineral_ok_mask"] = mineral_ok
    feats["mineral_lowsignal_mask"] = bare & ~mineral_ok

    # Relative thermal-extremes score: how Mars-like a cell is compared with the rest of Earth's land (large diurnal
    # range, cold nights). Target-independent; the absolute comparison is kept separately as the "Mars gap".
    wv = np.where(valid, w, 0.0)
    q["diurnal_range_k"] = quantile_table(feats["diurnal_range_k"], wv)
    q["lst_min_k"] = quantile_table(feats["lst_min_k"], wv)
    feats["extremes_rel"] = (0.5 * (to_percentile(feats["diurnal_range_k"], q["diurnal_range_k"])
                                    + (1.0 - to_percentile(feats["lst_min_k"], q["lst_min_k"])))).astype(np.float32)

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
