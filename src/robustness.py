"""Robustness of the candidate list to the weights.

Samples N weight vectors from a Dirichlet distribution centred on the default weights (concentration = 20 x defaults)
and recomputes the overall score on the full 0.05 deg grid. For every region that is above the top-0.5% threshold at
the default weights (the candidate pool, ~3000 regions per target) we recompute its mean score from its own cells and
record two shares over the samples:
  robust_pool_pct  - the region still scores above the top-0.5% threshold (it would still be a candidate)
  robust_top25_pct - the region is still among the 25 best-scoring regions of the pool
Regions that only emerge under other weights are not re-labelled (approximation, stated in REPORT.md).
Writes results/robustness_<target>.csv and adds both columns to results/candidates_<target>.csv.
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
from similarity import DIMS, overall_score  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
N_SAMPLES = 500
CONCENTRATION = 20.0
SEED = 0
N_PCT_SAMPLE = 2_000_000  # cells used to estimate the top-0.5% threshold per sample


def sample_weights(defaults, rng, n):
    keys = list(defaults)
    alpha = np.array([defaults[k] for k in keys]) * CONCENTRATION
    w = rng.dirichlet(alpha, size=n)
    return [dict(zip(keys, row)) for row in w]


def top_regions(ov, top_fraction, min_cells, n_top):
    """Boolean mask of the n_top best regions (by mean score) above the top_fraction threshold. Used by tests/demo."""
    thr = np.nanpercentile(ov, 100 * (1 - top_fraction))
    lab, n = ndimage.label(np.nan_to_num(ov, nan=-1) >= thr, structure=np.ones((3, 3)))
    if n == 0:
        return np.zeros_like(lab, dtype=bool)
    sums = ndimage.sum_labels(np.nan_to_num(ov), lab, index=np.arange(1, n + 1))
    counts = ndimage.sum_labels(np.ones_like(ov), lab, index=np.arange(1, n + 1))
    means = np.where(counts >= min_cells, sums / np.maximum(counts, 1), -1)
    keep = np.argsort(-means)[:n_top] + 1
    return np.isin(lab, keep[means[keep - 1] > 0])


def pool_regions(ov, cc):
    thr = np.nanpercentile(ov, 100 * (1 - cc["top_fraction"]))
    lab, n = ndimage.label(np.nan_to_num(ov, nan=-1) >= thr, structure=np.ones((3, 3)))
    counts = ndimage.sum_labels(np.ones_like(ov), lab, index=np.arange(1, n + 1))
    ids = np.flatnonzero(counts >= cc["min_cells"]) + 1
    lab[~np.isin(lab, ids)] = 0
    return lab, ids


def main(targets=None):
    cfg = yaml.safe_load(open(ROOT / "config.yaml"))
    fps = json.load(open(PROC / "mars_fingerprints.json"))
    cc = cfg["candidates"]
    rng = np.random.default_rng(SEED)
    ws = sample_weights(cfg["weights"], rng, N_SAMPLES)
    for t in targets or list(fps):
        dims = {d: np.load(PROC / "similarity" / f"{t}_{d}.npy") for d in DIMS}
        ov0 = overall_score(dims, cfg["weights"])
        lab, ids = pool_regions(ov0, cc)
        flat = lab.ravel()
        valid = np.flatnonzero(np.isfinite(ov0.ravel()))
        sub = rng.choice(valid, min(N_PCT_SAMPLE, valid.size), replace=False)
        in_pool = np.zeros(len(ids))
        in_top = np.zeros(len(ids))
        for w in ws:
            ov = overall_score(dims, w).ravel()
            thr = np.nanpercentile(ov[sub], 100 * (1 - cc["top_fraction"]))
            means = ndimage.mean(np.nan_to_num(ov), flat, index=ids)
            in_pool += means >= thr
            order = np.argsort(-means)[: cc["n_candidates"]]
            in_top[order] += 1
        pool = pd.DataFrame({"region_id": ids, "robust_pool_pct": np.round(100 * in_pool / N_SAMPLES, 1),
                             "robust_top25_pct": np.round(100 * in_top / N_SAMPLES, 1)})
        df = pd.read_csv(ROOT / "results" / f"candidates_{t}.csv")
        df = df.drop(columns=[c for c in df.columns if c.startswith("robust")]).merge(pool, on="region_id", how="left")
        df.to_csv(ROOT / "results" / f"candidates_{t}.csv", index=False)
        pool.to_csv(ROOT / "results" / f"robustness_{t}.csv", index=False)
        top = df[df.list == "top"]
        print(f"{t}: top-25 regions keep candidate status in median {top.robust_pool_pct.median():.0f}% of weight samples "
              f"and stay in the top 25 in median {top.robust_top25_pct.median():.0f}%")


if __name__ == "__main__":
    main(sys.argv[1:] or None)
