"""MVP-0 (Earth): print structure, attributes, ranges and fill values of the Earth products; verify LWmap codes."""
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "earth"


def h5tree(path):
    with h5py.File(path, "r") as h:
        def show(name, obj):
            if isinstance(obj, h5py.Dataset):
                a = obj[()]
                att = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in obj.attrs.items()}
                print(f"  {name:32s} {a.shape} {a.dtype} min={a.min()} max={a.max()} attrs={att}")
        h.visititems(show)
        lw = h["Land Water Map/LWmap"][()]
        dem = h["ASTER GDEM/ASTGDEM"][()]
        for code in np.unique(lw):
            m = lw == code
            print(f"  LWmap code {code}: n={m.sum()} dem mean={dem[m].mean():.1f} dem==0 frac={(dem[m] == 0).mean():.2f}")


def hdf4(path):
    from pyhdf.SD import SD, SDC
    f = SD(str(path), SDC.READ)
    for name in sorted(f.datasets()):
        s = f.select(name)
        a = s.get()
        att = {k: v for k, v in s.attributes().items() if k in ("scale_factor", "add_offset", "_FillValue", "valid_range", "units", "long_name")}
        print(f"  {name:28s} {a.shape} {a.dtype} min={a.min()} max={a.max()} {att}")
    f.end()


if __name__ == "__main__":
    for p in sys.argv[1:]:
        print(p)
        (h5tree if p.endswith(".h5") else hdf4)(p)
