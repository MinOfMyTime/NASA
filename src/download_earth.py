"""Download the Earth products with earthaccess (credentials from ~/.netrc or EARTHDATA_USERNAME/PASSWORD).

  python src/download_earth.py modis          # MYD11C3 (12 monthly) + MCD43C3 (days in config.mcd43c3_days)
  python src/download_earth.py ag1km          # ASTER GED 1 km, all tiles (~25k x 0.46 MB)
  python src/download_earth.py ag100 LAT LON  # ASTER GED 100 m tile(s) around a point (case studies)
  python src/download_earth.py emit LAT LON   # EMIT L2B mineral granules intersecting a 0.5 deg box
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import earthaccess
import yaml

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "earth"
CFG = yaml.safe_load(open(ROOT / "config.yaml"))


def login():
    strategy = "environment" if os.environ.get("EARTHDATA_USERNAME") else "netrc"
    auth = earthaccess.login(strategy=strategy)
    if not auth or not auth.authenticated:
        raise SystemExit(f"Earthdata login failed (strategy={strategy}). Check credentials / account status.")
    return auth


def fetch(results, sub, threads=8):
    out = RAW / sub
    out.mkdir(parents=True, exist_ok=True)
    have = {p.name for p in out.iterdir()}
    todo = [g for g in results if Path(g.data_links()[0]).name not in have]
    print(f"{sub}: {len(results)} granules, {len(todo)} to download")
    if todo:
        earthaccess.download(todo, str(out), threads=threads)
    return out


def modis():
    y = CFG["year"]
    r = earthaccess.search_data(short_name="MYD11C3", version="061", temporal=(f"{y}-01-01", f"{y}-12-31"))
    assert len(r) == 12, f"expected 12 MYD11C3 granules for {y}, got {len(r)}"
    fetch(r, "MYD11C3")
    picks = []
    for m in range(1, 13):
        for d in CFG["mcd43c3_days"]:
            day = f"{y}-{m:02d}-{d:02d}"
            g = earthaccess.search_data(short_name="MCD43C3", version="061", temporal=(day, day))
            # temporal search can return the neighbouring day as well; keep the exact day-of-year
            import datetime as dt

            doy = dt.date(y, m, d).timetuple().tm_yday
            g = [x for x in g if f".A{y}{doy:03d}." in x.data_links()[0]]
            assert len(g) == 1, f"MCD43C3 {day}: found {len(g)}"
            picks += g
    fetch(picks, "MCD43C3")


def ag1km():
    r = earthaccess.search_data(short_name="AG1KM", version="003", count=-1)
    print("AG1KM tiles:", len(r))
    fetch(r, "AG1KM", threads=16)


def ag100(lat, lon, half=0.5):
    r = earthaccess.search_data(short_name="AG100", version="003",
                                bounding_box=(lon - half, lat - half, lon + half, lat + half))
    fetch(r, "AG100")


def emit(lat, lon, half=0.25):
    r = earthaccess.search_data(short_name="EMITL2BMIN", version="001",
                                bounding_box=(lon - half, lat - half, lon + half, lat + half), count=20)
    print(f"EMIT L2B MIN granules near ({lat},{lon}): {len(r)}")
    fetch(r[:4], "EMITL2BMIN")
    return r


if __name__ == "__main__":
    login()
    what = sys.argv[1]
    if what == "modis":
        modis()
    elif what == "ag1km":
        ag1km()
    elif what == "ag100":
        ag100(float(sys.argv[2]), float(sys.argv[3]))
    elif what == "emit":
        emit(float(sys.argv[2]), float(sys.argv[3]))
    else:
        raise SystemExit(__doc__)
