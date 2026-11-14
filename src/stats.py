"""Area-weighted percentile tables and robust summary statistics."""
from __future__ import annotations

import numpy as np

N_Q = 1001


def area_weights(lats: np.ndarray, shape) -> np.ndarray:
    return np.broadcast_to(np.cos(np.radians(lats))[:, None], shape)


def quantile_table(values: np.ndarray, weights: np.ndarray | None = None, n=N_Q) -> np.ndarray:
    """Values at n evenly spaced cumulative-area fractions (0..1). Used as the percentile transform for a planet."""
    v = values.ravel()
    w = np.ones_like(v) if weights is None else np.asarray(weights).ravel()
    m = np.isfinite(v) & (w > 0)
    v, w = v[m], w[m]
    o = np.argsort(v, kind="stable")
    v, w = v[o], w[o]
    cw = np.cumsum(w)
    cw = (cw - 0.5 * w) / cw[-1]
    return np.interp(np.linspace(0, 1, n), cw, v).astype(np.float64)


def to_percentile(x: np.ndarray, table: np.ndarray) -> np.ndarray:
    """Map values to 0..1 cumulative-area percentile using a quantile table (ties -> mid-rank)."""
    q = np.linspace(0, 1, table.size)
    lo = np.interp(x, table, q, left=0.0, right=1.0)
    # mid-rank for flat stretches (e.g. many exact zeros in TES mineral abundances)
    left = np.searchsorted(table, x, side="left") / (table.size - 1)
    right = np.searchsorted(table, x, side="right") / (table.size - 1)
    tie = right - left > 1.0 / (table.size - 1)
    out = np.where(tie, 0.5 * (left + right), lo)
    out = np.clip(out, 0, 1).astype(np.float32)
    out[~np.isfinite(x)] = np.nan
    return out


def robust_summary(v: np.ndarray) -> dict:
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {"n": 0, "median": None, "q25": None, "q75": None, "iqr": None}
    q25, med, q75 = np.percentile(v, [25, 50, 75])
    return {"n": int(v.size), "median": float(med), "q25": float(q25), "q75": float(q75), "iqr": float(q75 - q25)}
