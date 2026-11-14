"""Per-dimension similarity of every Earth cell to a Mars target fingerprint, plus the weighted overall score.

  s_f = robust distance of feature f: (x - target_median) / max(target_IQR, min_scale_f)
  d^2 = mean_f s_f^2 over the features of a dimension
  s_d = exp(-d^2 / (2 sigma^2))  in [0, 1]
  overall = sum_d w_d s_d / sum_d w_d   over dimensions with data (weights renormalised per cell)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from features_mars import COMPARE  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
DIMS = list(COMPARE)

SCALE_KEY = {
    "log_relief": "log_relief", "slope_deg": "slope_deg", "log_rough": "log_rough", "tpi_m": "tpi_m",
    "pct_ti": "pct", "albedo": "albedo", "comp_mafic": "pct", "comp_silica": "pct", "comp_fines": "pct",
    "ndvi": "ndvi", "diurnal_range_k": "diurnal_range_k", "lst_min_k": "lst_min_k", "lst_max_k": "lst_max_k",
}


def feature_distance(x, summary, min_scale):
    scale = max(summary["iqr"] or 0.0, min_scale)
    return (x - summary["median"]) / scale


def ndvi_distance(ndvi, min_scale):
    """Mars NDVI = 0. Negative NDVI (water/snow) is masked upstream; clip so bare ground (NDVI<=0) scores 1."""
    return np.clip(ndvi, 0, None) / min_scale


def dimension_score(feats: dict, fp: dict, dim: str, cfg) -> np.ndarray:
    sig = cfg["similarity"]["sigma"]
    ms = cfg["similarity"]["min_scale"]
    d2 = None
    n = 0
    for f in COMPARE[dim]:
        if f not in feats:
            continue
        if f == "ndvi":
            z = ndvi_distance(feats[f], ms["ndvi"])
        else:
            z = feature_distance(feats[f], fp["features"][f], ms[SCALE_KEY[f]])
        d2 = z * z if d2 is None else d2 + z * z
        n += 1
    if n == 0:
        raise KeyError(f"no features for dimension {dim}")
    return np.exp(-(d2 / n) / (2 * sig * sig)).astype(np.float32)


def overall_score(scores: dict, weights: dict, require_all: bool = True) -> np.ndarray:
    """Weighted mean of dimension scores. With require_all (default) a cell whose weighted dimension is n/a gets n/a:
    otherwise a missing dimension (e.g. mineral signal unavailable) would silently hand its weight to the others and
    such cells would win the ranking (this happened in round 1, see REPORT.md)."""
    num = 0.0
    den = 0.0
    missing = None
    for d, s in scores.items():
        w = float(weights.get(d, 0.0))
        ok = np.isfinite(s)
        num = num + np.where(ok, s, 0.0) * w
        den = den + ok * w
        if w > 0:
            missing = ~ok if missing is None else (missing | ~ok)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (num / den).astype(np.float32)
    out[den == 0] = np.nan
    if require_all and missing is not None:
        out[missing] = np.nan
    return out


def load_earth():
    d = PROC / "earth"
    feats = {p.stem: np.load(p, mmap_mode="r") for p in d.glob("*.npy")}
    return feats


def main(targets=None):
    cfg = yaml.safe_load(open(ROOT / "config.yaml"))
    fps = json.load(open(PROC / "mars_fingerprints.json"))
    feats = load_earth()
    valid = feats["valid_mask"].astype(bool)
    out = PROC / "similarity"
    out.mkdir(exist_ok=True)
    for t in targets or fps:
        scores = {}
        for dim in DIMS:
            if dim == "extremes":
                # Weighted score uses the RELATIVE extremes score (Mars-likeness vs the rest of Earth's land).
                # The absolute comparison is exported separately as the "Mars gap" and is not weighted.
                gap = dimension_score(feats, fps[t], dim, cfg)
                gap[~valid] = np.nan
                np.save(out / f"{t}_gap.npy", gap)
                s = np.array(feats["extremes_rel"], dtype=np.float32).copy()
            else:
                s = dimension_score(feats, fps[t], dim, cfg)
            if dim == "mineral":
                # vegetated land (NDVI >= ndvi_max_for_mineral) has no exposed mineral surface to match -> 0;
                # bare cells with no usable TIR signal (flat spectrum / too few ASTER looks) stay NaN = n/a
                low = np.asarray(feats["mineral_lowsignal_mask"]).astype(bool)
                s[valid & ~np.isfinite(s) & ~low] = 0.0
            s[~valid] = np.nan
            scores[dim] = s
            np.save(out / f"{t}_{dim}.npy", s)
        ov = overall_score(scores, cfg["weights"])
        np.save(out / f"{t}_overall.npy", ov)
        print(f"{t}: overall p50={np.nanmedian(ov):.3f} p99={np.nanpercentile(ov, 99):.3f} max={np.nanmax(ov):.3f}")


if __name__ == "__main__":
    main(sys.argv[1:] or None)
