"""Earth DEM for topography: Copernicus GLO-90 (public AWS bucket, no login), averaged to the 1/30 deg grid
(3.706 km N-S, matching the MOLA MEGDR pixel). Reads only the COG's 1/4 overview (~333 m) for each 1x1 deg tile and
area-averages it to 30x30 cells per degree. Replaces the ASTER GED GDEM layer, which is all-fill in 4,944 of 24,873
AG1KM tiles (scripts/gap_stats.py), including Hawaii.

Output: data/processed/earth_dem_30ppd.npy (5400 x 10800 float32, north-up, -180..180; NaN = no tile = ocean).
"""
from __future__ import annotations

import os
import re
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "5")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "2")
import rasterio  # noqa: E402
from rasterio.enums import Resampling  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BUCKET = "https://copernicus-dem-90m.s3.amazonaws.com"
OUT = ROOT / "data" / "processed" / "earth_dem_30ppd.npy"
PPD = 30


def tile_origin(name):
    """Copernicus tile names give the SOUTH-WEST corner: ..._N19_00_W156_00_DEM covers 19..20N, 156..155W."""
    m = re.search(r"_([NS])(\d+)_00_([EW])(\d+)_00_DEM", name)
    lat = int(m.group(2)) * (1 if m.group(1) == "N" else -1)
    lon = int(m.group(4)) * (1 if m.group(3) == "E" else -1)
    return lat, lon


def read_tile(name):
    url = f"/vsicurl/{BUCKET}/{name}/{name}.tif"
    for _ in range(3):
        try:
            with rasterio.open(url) as d:
                a = d.read(1, out_shape=(PPD, PPD), resampling=Resampling.average).astype(np.float32)
            return name, a
        except Exception as e:  # noqa: BLE001
            err = e
    return name, err


def main():
    names = urllib.request.urlopen(f"{BUCKET}/tileList.txt", timeout=60).read().decode().split()
    print("Copernicus GLO-90 tiles:", len(names))
    dem = np.full((180 * PPD, 360 * PPD), np.nan, np.float32)
    failed = []
    with ThreadPoolExecutor(48) as ex:
        for i, (name, a) in enumerate(ex.map(read_tile, names), 1):
            if isinstance(a, Exception):
                failed.append((name, str(a)))
                continue
            lat, lon = tile_origin(name)
            r0 = (90 - (lat + 1)) * PPD
            c0 = (lon + 180) * PPD
            dem[r0:r0 + PPD, c0:c0 + PPD] = a
            if i % 2000 == 0:
                print(f"  {i}/{len(names)}", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.save(OUT, dem)
    print(f"saved {OUT.name}; failed tiles: {len(failed)}")
    for f in failed[:20]:
        print("  FAILED", f)
    if failed:
        raise SystemExit("some tiles failed; rerun")


if __name__ == "__main__":
    main()
