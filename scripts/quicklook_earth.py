"""MVP-0 (Earth): MODIS-derived layers -> quicklook PNGs + range printout (runs before the ASTER mosaic is ready)."""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from features_earth import modis_features  # noqa: E402

OUT = ROOT / "results" / "quicklooks"
cfg = yaml.safe_load(open(ROOT / "config.yaml"))
f = modis_features(cfg)
spec = {"albedo": (0, 0.6), "ati": (0, 0.06), "diurnal_range_k": (0, 40), "lst_min_k": (220, 300), "lst_max_k": (250, 340)}
for k, (lo, hi) in spec.items():
    a = f[k]
    print(f"{k:16s} valid={np.isfinite(a).mean():.3f} p1={np.nanpercentile(a, 1):.4g} p50={np.nanmedian(a):.4g} p99={np.nanpercentile(a, 99):.4g}")
    fig, ax = plt.subplots(figsize=(10, 5))
    m = ax.imshow(a[::2, ::2], extent=[-180, 180, -90, 90], cmap="viridis", vmin=lo, vmax=hi)
    ax.set_title(f"Earth {k} ({cfg['year']})")
    fig.colorbar(m, ax=ax, shrink=0.7)
    fig.savefig(OUT / f"earth_{k}.png", dpi=80, bbox_inches="tight")
    plt.close(fig)
