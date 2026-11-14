"""Candidate extraction: top 0.5% of valid land cells -> 8-connected regions -> ranked list.

Regions smaller than candidates.min_cells are dropped as speckle. Regions are ranked by mean overall score; area is
reported as a column. Regions within `exclude_known_km` of a known analog are flagged `near_known`;
the remainder are the 'under-studied' discovery candidates.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo import EARTH_RADIUS_KM, Grid, haversine_km, km_per_deg  # noqa: E402
from similarity import DIMS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
NE_PATH = ROOT / "data" / "raw" / "earth" / "ne_50m_admin_0_countries.zip"
NE_URL = "https://naciscdn.org/naturalearth/50m/cultural/ne_50m_admin_0_countries.zip"


def countries():
    import geopandas as gpd

    if not NE_PATH.exists():
        import urllib.request

        NE_PATH.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(NE_URL, headers={"User-Agent": "analogfinder"})
        NE_PATH.write_bytes(urllib.request.urlopen(req, timeout=120).read())
    return gpd.read_file(NE_PATH)[["NAME", "geometry"]]


def country_of(gdf, lat, lon):
    from shapely.geometry import Point

    p = Point(lon, lat)
    hit = gdf[gdf.contains(p)]
    if len(hit):
        return hit.iloc[0]["NAME"]
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        d = gdf.distance(p)  # degrees; only used to label offshore/coastal centroids
    return gdf.iloc[int(d.idxmin())]["NAME"] + " (nearest)"


def extract(t, cfg, gdf=None):
    g = Grid(cfg["grids"]["earth_ppd"], lon0=-180.0)
    ov = np.load(PROC / "similarity" / f"{t}_overall.npy")
    dims = {d: np.load(PROC / "similarity" / f"{t}_{d}.npy", mmap_mode="r") for d in DIMS}
    thr = np.nanpercentile(ov, 100 * (1 - cfg["candidates"]["top_fraction"]))
    lab, n = ndimage.label(np.nan_to_num(ov, nan=-1) >= thr, structure=np.ones((3, 3)))
    lats, lons = g.lats, g.lons
    cell_km2 = (km_per_deg(EARTH_RADIUS_KM) * g.res) ** 2 * np.cos(np.radians(lats))
    objs = ndimage.find_objects(lab)
    rows = []
    for i, sl in enumerate(objs, start=1):
        m = lab[sl] == i
        rr, cc = np.nonzero(m)
        rr = rr + sl[0].start
        cc = cc + sl[1].start
        s = ov[rr, cc]
        w = s / s.sum()
        rows.append({
            "region_id": i, "n_cells": int(m.sum()),
            "area_km2": float(cell_km2[rr].sum()),
            "lat": float((lats[rr] * w).sum()), "lon": float((lons[cc] * w).sum()),
            "mean_score": float(s.mean()), "max_score": float(s.max()),
            **{f"s_{d}": float(np.nanmean(np.asarray(dims[d])[rr, cc])) for d in DIMS},
        })
    df = pd.DataFrame(rows)
    df = df[df.n_cells >= cfg["candidates"]["min_cells"]]
    # Ranked by mean overall score (area is shown as a column, not used for ranking)
    df["rank_metric"] = df["mean_score"]
    df = df.sort_values("rank_metric", ascending=False).reset_index(drop=True)
    known = cfg["known_analogs"]
    near, dist = [], []
    for _, r in df.iterrows():
        ds = {k: float(haversine_km(r.lat, r.lon, v["lat"], v["lon"])) for k, v in known.items()}
        k = min(ds, key=ds.get)
        dist.append(round(ds[k], 1))
        near.append(known[k]["name"] if ds[k] <= cfg["candidates"]["exclude_known_km"] else "")
    df["nearest_known_km"] = dist
    df["near_known"] = near
    top = df.head(cfg["candidates"]["n_candidates"]).copy()
    disc = df[df.near_known == ""].head(cfg["candidates"]["n_candidates"]).copy()
    out = pd.concat([top.assign(list="top"), disc.assign(list="discovery")])  # a region may be in both lists
    if gdf is not None:
        out["country"] = [country_of(gdf, a, b) for a, b in zip(out.lat, out.lon)]
    out.insert(0, "rank", out.groupby("list").cumcount() + 1)
    out["threshold"] = thr
    return out, n


def main():
    cfg = yaml.safe_load(open(ROOT / "config.yaml"))
    fps = json.load(open(PROC / "mars_fingerprints.json"))
    gdf = countries()
    for t in fps:
        df, n = extract(t, cfg, gdf)
        p = ROOT / "results" / f"candidates_{t}.csv"
        df.round(4).to_csv(p, index=False)
        print(f"== {t}: {n} regions above threshold; top 5:")
        print(df[df.list == "top"].head(5)[["rank", "lat", "lon", "country", "area_km2", "mean_score", "near_known"]]
              .to_string(index=False))
        print("   top 5 discovery (>150 km from known analogs):")
        print(df[df.list == "discovery"].head(5)[["rank", "lat", "lon", "country", "area_km2", "mean_score"]]
              .to_string(index=False))


if __name__ == "__main__":
    main()
