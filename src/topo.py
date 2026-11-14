"""Topographic metrics computed IDENTICALLY on Earth and Mars on a regular lat/lon DEM whose N-S pixel size is
`pixel_km` (3.705 km on both planets, see config.yaml). E-W pixel size = pixel_km * cos(lat)."""
from __future__ import annotations

import numpy as np
from scipy import ndimage


def slope_deg(z: np.ndarray, lats: np.ndarray, pixel_km: float) -> np.ndarray:
    """Mean-gradient slope (degrees) from central differences; dx shrinks with cos(lat)."""
    dy = pixel_km * 1000.0
    dx = dy * np.clip(np.cos(np.radians(lats)), 0.05, None)[:, None]
    gz_y = (np.roll(z, -1, 0) - np.roll(z, 1, 0)) / (2 * dy)
    gz_x = (np.roll(z, -1, 1) - np.roll(z, 1, 1)) / (2 * dx)
    s = np.degrees(np.arctan(np.hypot(gz_x, gz_y)))
    s[0], s[-1] = np.nan, np.nan
    return s.astype(np.float32)


def local_relief(z: np.ndarray, window: int) -> np.ndarray:
    """max - min elevation inside a window x window box (m). NaNs propagate as NaN."""
    zf = np.where(np.isfinite(z), z, np.nan)
    hi = ndimage.maximum_filter(np.nan_to_num(zf, nan=-1e9), size=window, mode="wrap")
    lo = ndimage.minimum_filter(np.nan_to_num(zf, nan=1e9), size=window, mode="wrap")
    r = (hi - lo).astype(np.float32)
    bad = ndimage.maximum_filter((~np.isfinite(zf)).astype(np.uint8), size=window, mode="wrap") > 0
    r[bad] = np.nan
    return r


def roughness_detrended(z: np.ndarray, window: int) -> np.ndarray:
    """Std of elevation residuals after removing a least-squares plane in each window x window box (m).

    For a symmetric window with centred pixel offsets x, y (orthogonal), the plane fit z ~ c + a x + b y has
    c = E[z], a = E[zx]/E[x^2], b = E[zy]/E[y^2], and residual variance E[z^2] - c^2 - a^2 E[x^2] - b^2 E[y^2].
    Offsets are in pixel units; the residual is in metres and invariant to that choice.
    """
    h = window // 2
    off = np.arange(-h, h + 1, dtype=np.float64)
    kx = np.tile(off, (window, 1))  # columns
    ky = kx.T
    n = window * window
    zz = np.where(np.isfinite(z), z, 0.0).astype(np.float64)
    ez = ndimage.uniform_filter(zz, window, mode="wrap")
    ez2 = ndimage.uniform_filter(zz * zz, window, mode="wrap")
    # correlate with offset kernels (flip because correlate uses the kernel as-is around the centre)
    ezx = ndimage.correlate(zz, kx / n, mode="wrap")
    ezy = ndimage.correlate(zz, ky / n, mode="wrap")
    ex2 = np.mean(kx**2)
    var = ez2 - ez**2 - ezx**2 / ex2 - ezy**2 / ex2
    r = np.sqrt(np.clip(var, 0, None)).astype(np.float32)
    bad = ndimage.maximum_filter((~np.isfinite(z)).astype(np.uint8), size=window, mode="wrap") > 0
    r[bad] = np.nan
    return r


def tpi(z: np.ndarray, window: int) -> np.ndarray:
    """Topographic position index: elevation minus the window mean (m). Negative = basin/crater floor."""
    zz = np.where(np.isfinite(z), z, 0.0).astype(np.float64)
    w = ndimage.uniform_filter(np.isfinite(z).astype(np.float64), window, mode="wrap")
    m = ndimage.uniform_filter(zz, window, mode="wrap") / np.where(w > 0, w, np.nan)
    out = (z - m).astype(np.float32)
    return out


def topo_features(z: np.ndarray, lats: np.ndarray, pixel_km: float, window: int, tpi_window: int) -> dict:
    relief = local_relief(z, window)
    rough = roughness_detrended(z, window)
    slope = slope_deg(z, lats, pixel_km)
    slope_mean = ndimage.uniform_filter(np.nan_to_num(slope), window, mode="wrap").astype(np.float32)
    slope_mean[~np.isfinite(relief)] = np.nan
    h = max(window, tpi_window)
    out = {
        # log10(1 + x) compresses the heavy right tail so the robust distance is not dominated by mountains
        "log_relief": np.log10(1.0 + relief),
        "slope_deg": slope_mean,
        "log_rough": np.log10(1.0 + rough),
        "tpi_m": tpi(z, tpi_window),
    }
    for v in out.values():
        v[:h] = np.nan
        v[-h:] = np.nan  # no wrap across the poles
    return out
