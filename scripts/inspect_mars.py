"""MVP-0 (Mars): print shapes/ranges/fill values of every Mars product, verify orientation empirically and write
quicklook PNGs to results/quicklooks/."""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from geo import Grid, resample_nearest  # noqa: E402
import io_mars as im  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results" / "quicklooks"
OUT.mkdir(parents=True, exist_ok=True)
G4E = Grid(4, lon0=0.0)  # comparison grid, 0..360E


def stats(name, a):
    f = np.isfinite(a)
    print(f"{name:28s} shape={a.shape} finite={f.mean():.3f} min={np.nanmin(a):.4g} "
          f"p50={np.nanmedian(a):.4g} max={np.nanmax(a):.4g}")


def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[m], b[m])[0, 1])


def quick(name, a, grid, cmap="viridis", vmin=None, vmax=None):
    fig, ax = plt.subplots(figsize=(10, 5))
    ext = [grid.lon0, grid.lon0 + 360, -90, 90]
    im_ = ax.imshow(a, extent=ext, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_title(name)
    ax.set_xlabel("longitude (deg E)")
    ax.set_ylabel("latitude")
    fig.colorbar(im_, ax=ax, shrink=0.7)
    fig.savefig(OUT / f"mars_{name}.png", dpi=80, bbox_inches="tight")
    plt.close(fig)


# --- endianness check for the TI map (header says range 24..800) ---
for dt in (">i2", "<i2"):
    raw = np.fromfile(im.RAW / "nmap2003.bin", dtype=dt, offset=7200)
    nz = raw[raw != 0]
    print(f"TI as {dt}: min={nz.min()} max={nz.max()} frac_in_24_800={np.mean((nz >= 24) & (nz <= 800)):.4f}")

mola, gm = im.read_mola_megdr()
alb, ga = im.read_tes_albedo()
ti, gt = im.read_tes_ti()
ti_all, _ = im.read_tes_ti(drop_interpolated=False)
stats("MOLA topography (m)", mola)
stats("TES albedo", alb)
stats("TES TI (TES-only)", ti)
stats("TES TI (incl. interpolated)", ti_all)

alb4 = resample_nearest(alb, ga, G4E)
mola4 = resample_nearest(mola, gm, G4E)

# --- TI orientation check: TI vs albedo is strongly anti-correlated on Mars (dusty = bright + low TI) ---
raw = np.fromfile(im.RAW / "nmap2003.bin", dtype=">i2", offset=7200).reshape(3600, 7200).astype(np.float32)
raw[raw == 0] = np.nan
for label, arr in {
    "header (south-first, west-positive)": raw[::-1, ::-1],
    "alt: north-first, west-positive": raw[:, ::-1],
    "alt: south-first, east-positive": raw[::-1, :],
}.items():
    print(f"TI orientation {label:40s} corr(TI, albedo) = {corr(resample_nearest(arr, Grid(20, 0.0), G4E), alb4):+.3f}")

# --- point checks against known features ---
for nm, lat, lon in [("Olympus Mons summit", 18.65, 226.2), ("Hellas floor", -42.4, 70.5), ("Arsia Mons", -8.26, 239.0)]:
    r, c = gm.rowcol(lat, lon)
    print(f"MOLA at {nm:22s} ({lat}, {lon}E): {mola[r, c]:.0f} m")

# --- TES mineral maps: orientation unknown -> test hypotheses using dust vs albedo ---
dust, kv = im.read_vicar(im.RAW / "tes_minerals" / "tes_surface_dust.vic")
dust = dust.astype(np.float32)
print("VICAR label keys:", {k: kv[k] for k in ("NL", "NS", "FORMAT", "ORG")})
for nm in im.TES_MINERALS:
    a, _ = im.read_vicar(im.RAW / "tes_minerals" / f"tes_{nm}.vic")
    vals, counts = np.unique(a[~np.isfinite(a) | (np.abs(a) > 1e20)], return_counts=True) if False else (None, None)
    u = np.unique(a)
    print(f"  {nm:24s} min={a.min():.4g} max={a.max():.4g} n_unique={u.size} most common={np.bincount(np.searchsorted(u, a.ravel())).argmax() and u[np.bincount(np.searchsorted(u, a.ravel())).argmax()]}")
results = {}
for north_up in (True, False):
    for lon0 in (0.0, -180.0):
        for west in (False, True):
            d = dust if north_up else dust[::-1]
            d = d[:, ::-1] if west else d
            d = np.where((d >= 0) & (d < 1e20), d, np.nan)
            d4 = resample_nearest(d, Grid(4, lon0=lon0), G4E)
            results[(north_up, lon0, west)] = corr(d4, alb4)
for k, v in sorted(results.items(), key=lambda kv: -kv[1]):
    print(f"dust orientation north_up={k[0]!s:5s} lon0={k[1]:6.0f} west_positive={k[2]!s:5s}: corr(dust, albedo)={v:+.3f}")

quick("mola_topography_m", mola, gm, cmap="cividis")
quick("tes_albedo", alb, ga, vmin=0.08, vmax=0.32)
quick("tes_thermal_inertia", ti, gt, vmin=24, vmax=600)
best = max(results, key=results.get)
print("best mineral orientation:", best)
for nm in im.TES_MINERALS:
    a, g = im.read_tes_mineral(nm, orientation={"north_up": best[0], "lon0": best[1]})
    quick(f"tes_{nm}", a, g)

# --- independent check: hematite peaks at Meridiani Planum (~2S, 354.5E; Christensen et al. 2000) ---
for north_up, lon0 in [(True, 0.0), (True, -180.0), (False, 0.0)]:
    h, g = im.read_tes_mineral("hematite", orientation={"north_up": north_up, "lon0": lon0})
    rows, cols = g.window(-1.9, 354.5, 2.0)
    win = h[np.ix_(rows, cols)]
    print(f"hematite north_up={north_up} lon0={lon0}: Meridiani window median={np.nanmedian(win):.3f} "
          f"global p99={np.nanpercentile(h, 99):.3f} global median={np.nanmedian(h):.3f}")
    for nm in im.TES_MINERALS:
        a, _ = im.read_tes_mineral(nm, orientation={"north_up": north_up, "lon0": lon0})
        if north_up and lon0 == 0.0:
            print(f"   {nm:24s} coverage={np.isfinite(a).mean():.3f} zeros={np.mean(a == 0):.3f}")
