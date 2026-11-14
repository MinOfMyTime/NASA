"""Debug: elevation around a point from the 0.05 deg mosaic and from the raw AG1KM tile."""
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from geo import Grid  # noqa: E402

lat, lon = float(sys.argv[1]), float(sys.argv[2])
g = Grid(20, -180.0)
e = np.load(ROOT / "data/processed/earth/elevation_m.npy", mmap_mode="r")
r, c = g.rowcol(lat, lon)
print("mosaic 0.05 elevation 5x5:\n", np.asarray(e[r - 2:r + 3, c - 2:c + 3]).round(0))
for k in ["log_relief", "slope_deg", "land_frac", "valid_mask"]:
    print(k, np.asarray(np.load(ROOT / f"data/processed/earth/{k}.npy", mmap_mode="r")[r - 1:r + 2, c - 1:c + 2]).round(3).tolist())
top, left = int(np.ceil(lat)), int(np.floor(lon))
for f in (ROOT / "data/raw/earth/AG1KM").glob(f"AG1km.v003.{top}.{left:+04d}*.h5".replace("+", "")):
    with h5py.File(f) as h:
        d = h["ASTER GDEM/ASTGDEM"][()]
        print(f.name, "lat", h["Geolocation/Latitude"][()].min(), h["Geolocation/Latitude"][()].max(), "lon", h["Geolocation/Longitude"][()].min(), "dem max", d.max())
