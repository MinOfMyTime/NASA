"""Latitude histogram of the top-50 candidate regions per target (before/after artifact check)."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BINS = [-90, -60, -30, 0, 30, 60, 70, 90]


def main(tag):
    out = {}
    for t in ["jezero", "gale", "meridiani"]:
        df = pd.read_csv(ROOT / f"results/candidates_{t}.csv")
        df = df[df.list == "top"]
        h, _ = np.histogram(df.lat, bins=BINS)
        out[t] = {"bins": BINS, "counts": h.tolist(), "n": int(len(df)), "n_above_70": int((df.lat > 70).sum())}
        print(t, dict(zip([f"{a}..{b}" for a, b in zip(BINS[:-1], BINS[1:])], h.tolist())), "above 70N:", out[t]["n_above_70"])
    p = ROOT / "results" / f"lat_hist_{tag}.json"
    json.dump(out, open(p, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "after")
