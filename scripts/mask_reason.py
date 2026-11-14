"""Debug: which mask/feature removes a cell from the valid set."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from geo import Grid  # noqa: E402

g = Grid(20, -180.0)
keys = ["land_frac", "emis_b10", "emis_b13", "lst_max_k", "ati", "albedo", "ndvi", "valid_mask", "ice_mask"]
E = {k: np.load(ROOT / f"data/processed/earth/{k}.npy", mmap_mode="r") for k in keys}
for lat, lon in [(float(a), float(b)) for a, b in zip(sys.argv[1::2], sys.argv[2::2])]:
    r, c = g.rowcol(lat, lon)
    print(f"({lat},{lon}): " + "  ".join(f"{k}={float(E[k][r, c]):.3f}" for k in keys))
