"""Readers for the raw Mars products. Every reader returns (array float32 with NaN fill, Grid) on a
north-up, east-positive grid. See data/SOURCES.md for provenance."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from geo import Grid

RAW = Path(__file__).resolve().parents[1] / "data" / "raw" / "mars"


def read_mola_megdr(path=RAW / "megt90n000eb.img"):
    """MOLA MEGDR MEGT90N000EB (16 ppd). PDS label: 2880x5760 MSB_INTEGER metres, line 1 = 90N,
    sample 1 western edge = 0E, east-positive. No fill value (gaps are interpolated by the producer)."""
    a = np.fromfile(path, dtype=">i2").reshape(2880, 5760).astype(np.float32)
    return a, Grid(16, lon0=0.0)


def read_tes_albedo(path=RAW / "Mars_MGS_TES_Albedo_mosaic_global_7410m.tif"):
    """USGS TES bolometric albedo mosaic, 2880x1440 (8 ppd) simple cylindrical, -180..180E, north up."""
    import rasterio

    with rasterio.open(path) as ds:
        a = ds.read(1).astype(np.float32)
        nod = ds.nodata
        t = ds.transform
    assert abs(t.c + 10670400.0) < 1 and a.shape == (1440, 2880), "unexpected albedo geometry"
    if nod is not None:
        a[a <= nod * 0.999] = np.nan  # nodata is ~-3.4e38
    a[(a < 0) | (a > 1)] = np.nan
    return a, Grid(8, lon0=-180.0)


def read_tes_ti(path=RAW / "nmap2003.bin", mask_path=RAW / "nmsk2003.bin", drop_interpolated=True):
    """Putzig & Mellon (2005) nightside TES thermal inertia, 20 ppd.

    File: 7200-byte text header + 3600 records of 7200 int16 (big-endian, verified by value range 24..800).
    Header: 'Latitudes: -90 to 90' (first record = southernmost row) and 'Longitude: 0 to 360 West'.
    0 = no data. Companion mask: 1 = real TES inertia, 0 = no data/interpolated.
    Output: north-up, east-positive 0..360 grid.
    """
    raw = np.fromfile(path, dtype=">i2", offset=7200).reshape(3600, 7200).astype(np.float32)
    raw[raw == 0] = np.nan
    if drop_interpolated:
        m = np.fromfile(mask_path, dtype=np.uint8, offset=7200).reshape(3600, 7200)
        raw[m != 1] = np.nan
    a = flip_ti_orientation(raw)
    return a, Grid(20, lon0=0.0)


def flip_ti_orientation(raw):
    """South-up, west-positive (col j covers W lon [j, j+1)/20) -> north-up, east-positive 0..360.
    East-positive column k covers E lon [k, k+1)/20 = W lon [7200-k-1, 7200-k)/20, i.e. a reversed column order."""
    return raw[::-1, ::-1].copy()


TES_MINERALS = [
    "plagioclase",
    "high_ca_pyroxene",
    "low_ca_pyroxene",
    "olivine",
    "quartz",
    "sheet_silicates_glass",
    "surface_dust",
    "hematite",
]


def read_vicar(path):
    """Minimal VICAR reader for the single-band REAL images supplied (LBLSIZE header, BSQ, little-endian)."""
    with open(path, "rb") as f:
        head = f.read(1024).decode("ascii", "replace")
    kv = {}
    for tok in head.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            kv[k] = v.strip("'")
    lbl = int(kv["LBLSIZE"])
    nl, ns = int(kv["NL"]), int(kv["NS"])
    assert kv["FORMAT"] == "REAL" and kv.get("REALFMT", "RIEEE") == "RIEEE"
    a = np.fromfile(path, dtype="<f4", offset=lbl, count=nl * ns).reshape(nl, ns)
    return a, kv


def read_tes_mineral(name, path=None, orientation=None):
    """TES mineral map (4 ppd, 1440x720 float32). The VICAR label carries no map geometry; the orientation
    (row order and longitude origin) is established empirically in scripts/inspect_mars.py by correlating the
    surface-dust map against the TES albedo mosaic, and stored in config.yaml (mars.tes_mineral_orientation)."""
    path = path or RAW / "tes_minerals" / f"tes_{name}.vic"
    a, _ = read_vicar(path)
    a = a.astype(np.float32)
    a[~np.isfinite(a)] = np.nan
    a[(a < 0) | (a > 1e30)] = np.nan  # -1 = no data in the supplied VICAR maps
    orientation = orientation or {"north_up": True, "lon0": 0.0}
    if not orientation["north_up"]:
        a = a[::-1]
    return a, Grid(4, lon0=float(orientation["lon0"]))
