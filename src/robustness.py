"""Robustness of the candidate list to the weights.

Samples N weight vectors from a Dirichlet distribution centred on the default weights (concentration = 20 x defaults),
recomputes the overall score on the 0.25 deg web grid (the same per-dimension grids the app uses), re-extracts the top-0.5%
regions and records, for every baseline candidate, the share of samples in which its centroid lies inside one of the top-25
regions. Writes results/robustness_<target>.csv and adds a `robust_pct` column to results/candidates_<target>.csv.
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
from geo import block_reduce_nan  # noqa: E402
from similarity import DIMS, overall_score  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
N_SAMPLES = 500
CONCENTRATION = 20.0
SEED = 0


def web_grids(t, cfg):
    ppd = cfg["grids"]["earth_ppd"]
    f = ppd // cfg["grids"]["web_ppd"]
    lat_lo, lat_hi = cfg["grids"]["web_lat_range"]
    r0, r1 = int((90 - lat_hi) * ppd), int((90 - lat_lo) * ppd)
    return {d: block_reduce_nan(np.load(PROC / "similarity" / f"{t}_{d}.npy", mmap_mode="r")[r0:r1].astype(np.float32), f)
            for d in DIMS}


def top_regions(ov, top_fraction, min_cells, n_top):
    thr = np.nanpercentile(ov, 100 * (1 - top_fraction))
    lab, n = ndimage.label(np.nan_to_num(ov, nan=-1) >= thr, structure=np.ones((3, 3)))
    if n == 0:
        return np.zeros_like(lab, dtype=bool)
    sums = ndimage.sum_labels(np.nan_to_num(ov), lab, index=np.arange(1, n + 1))
    counts = ndimage.sum_labels(np.ones_like(ov), lab, index=np.arange(1, n + 1))
    means = np.where(counts >= min_cells, sums / np.maximum(counts, 1), -1)
    keep = np.argsort(-means)[:n_top] + 1
    return np.isin(lab, keep[means[keep - 1] > 0])


def sample_weights(defaults, rng, n):
    keys = list(defaults)
    alpha = np.array([defaults[k] for k in keys]) * CONCENTRATION
    w = rng.dirichlet(alpha, size=n)
    return [dict(zip(keys, row)) for row in w]


def main(targets=None):
    cfg = yaml.safe_load(open(ROOT / "config.yaml"))
    fps = json.load(open(PROC / "mars_fingerprints.json"))
    rng = np.random.default_rng(SEED)
    ws = sample_weights(cfg["weights"], rng, N_SAMPLES)
    ppd = cfg["grids"]["web_ppd"]
    lat_hi = cfg["grids"]["web_lat_range"][1]
    cc = cfg["candidates"]
    min_cells_web = max(1, cc["min_cells"] // (cfg["grids"]["earth_ppd"] // ppd) ** 2)
    for t in targets or list(fps):
        grids = web_grids(t, cfg)
        df = pd.read_csv(ROOT / "results" / f"candidates_{t}.csv")
        rows = np.clip(((lat_hi - df.lat) * ppd).astype(int), 0, grids["topo"].shape[0] - 1)
        cols = np.clip(((df.lon + 180) * ppd).astype(int), 0, grids["topo"].shape[1] - 1)
        hits = np.zeros(len(df))
        for w in ws:
            ov = overall_score(grids, w)
            top = top_regions(ov, cc["top_fraction"], min_cells_web, cc["n_candidates"])
            hits += top[rows, cols]
        df["robust_pct"] = np.round(100 * hits / N_SAMPLES, 1)
        df.to_csv(ROOT / "results" / f"candidates_{t}.csv", index=False)
        df[["list", "rank", "region_id", "lat", "lon", "country", "mean_score", "robust_pct"]].to_csv(
            ROOT / "results" / f"robustness_{t}.csv", index=False)
        top = df[df.list == "top"]
        print(f"{t}: robust_pct of top-25 median={top.robust_pct.median():.0f}% "
              f"(>=80%: {(top.robust_pct >= 80).sum()}, <50%: {(top.robust_pct < 50).sum()})")


if __name__ == "__main__":
    main(sys.argv[1:] or None)
