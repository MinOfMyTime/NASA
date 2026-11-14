"""Print raw Earth feature values at the known analog sites next to each Mars target fingerprint (debug aid)."""
import json
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from geo import Grid  # noqa: E402

cfg = yaml.safe_load(open(ROOT / "config.yaml"))
fps = json.load(open(ROOT / "data/processed/mars_fingerprints.json"))
E = {p.stem: np.load(p, mmap_mode="r") for p in (ROOT / "data/processed/earth").glob("*.npy")}
g = Grid(20, -180.0)
keys = sys.argv[1:] or ["ati", "pct_ti", "albedo", "ndvi", "comp_mafic", "comp_silica", "comp_fines", "log_relief", "slope_deg"]
print(f"{'site':22s}" + "".join(f"{k:>12s}" for k in keys))
for t, fp in fps.items():
    print(f"{'MARS ' + t:22s}" + "".join(f"{fp['features'][k]['median']:12.3f}" if k in fp["features"] else f"{'':12s}" for k in keys))
for k_, s in cfg["known_analogs"].items():
    r, c = g.rowcol(s["lat"], s["lon"])
    print(f"{s['name'][:22]:22s}" + "".join(f"{float(E[k][r, c]):12.3f}" for k in keys))
