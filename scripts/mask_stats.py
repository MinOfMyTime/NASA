"""Debug: TIR contrast / NumObs distributions among bare land, by latitude band, to set the mineral masks."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from geo import Grid  # noqa: E402

E = {k: np.load(ROOT / f"data/processed/earth/{k}.npy") for k in ["mineral_bare_mask", "tir_contrast", "aster_nobs", "valid_mask"]}
lat = np.broadcast_to(Grid(20, -180.0).lats[:, None], E["tir_contrast"].shape)
bare = E["mineral_bare_mask"]
for lo, hi in [(-60, -30), (-30, 0), (0, 30), (30, 60), (60, 70), (70, 84)]:
    m = bare & (lat >= lo) & (lat < hi)
    c, n = E["tir_contrast"][m], E["aster_nobs"][m]
    print(f"lat {lo:4d}..{hi:3d}: bare cells {m.sum():8,}  contrast p10/p50/p90 = {np.nanpercentile(c, 10):.4f}/{np.nanmedian(c):.4f}/{np.nanpercentile(c, 90):.4f}"
          f"  frac<0.02 = {np.nanmean(c < 0.02):.2f}  frac<0.01 = {np.nanmean(c < 0.01):.2f}  nobs p10/p50 = {np.nanpercentile(n, 10):.1f}/{np.nanmedian(n):.1f}  frac nobs<3 = {np.nanmean(n < 3):.2f}")
