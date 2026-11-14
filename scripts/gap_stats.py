"""Count AG1KM tiles whose GDEM layer is all fill, and land cells with invalid emissivity (>1) in the mosaic."""
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
files = sorted((ROOT / "data/raw/earth/AG1KM").glob("*.h5"))
n_fill, n_land_fill, n_em_hi = 0, 0, 0
names = []
for f in files:
    with h5py.File(f) as h:
        d = h["ASTER GDEM/ASTGDEM"][()]
        lw = h["Land Water Map/LWmap"][()]
        em = h["Emissivity/Mean"][()]
    land = lw == 0
    if land.any() and (d[land] == -9999).all():
        n_land_fill += 1
        names.append(f.name)
    n_em_hi += int(((em > 1000) & land[None]).any())
print(f"tiles: {len(files)}; tiles with land but all-fill GDEM: {n_land_fill}; tiles with emissivity > 1 on land: {n_em_hi}")
print(names[:40])
