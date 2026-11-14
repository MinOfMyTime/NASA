"""Shared geographic helpers: longitude conventions, regular lat/lon grids, planet constants."""
from __future__ import annotations

import numpy as np

MARS_RADIUS_KM = 3396.0  # IAU2000 sphere used by MOLA MEGDR / USGS mosaics
EARTH_RADIUS_KM = 6371.0


def km_per_deg(radius_km: float) -> float:
    return 2 * np.pi * radius_km / 360.0


def lon_e360_to_e180(lon):
    """East-positive 0..360 (Mars convention) -> east-positive -180..180."""
    lon = np.asarray(lon, dtype=float)
    return ((lon + 180.0) % 360.0) - 180.0


def lon_e180_to_e360(lon):
    """East-positive -180..180 -> east-positive 0..360."""
    return np.asarray(lon, dtype=float) % 360.0


def lon_w360_to_e360(lon_w):
    """West-positive 0..360 (e.g. Putzig & Mellon TI map) -> east-positive 0..360."""
    return (360.0 - np.asarray(lon_w, dtype=float)) % 360.0


class Grid:
    """Regular global lat/lon grid with pixel-centre coordinates.

    Rows run north -> south. `lon0` is the western edge of column 0 (0 for 0..360, -180 for -180..180).
    """

    def __init__(self, ppd: float, lon0: float = -180.0):
        self.ppd = ppd
        self.res = 1.0 / ppd
        self.lon0 = lon0
        self.nrow = int(round(180 * ppd))
        self.ncol = int(round(360 * ppd))

    @property
    def shape(self):
        return (self.nrow, self.ncol)

    @property
    def lats(self):
        return 90.0 - (np.arange(self.nrow) + 0.5) * self.res

    @property
    def lons(self):
        return self.lon0 + (np.arange(self.ncol) + 0.5) * self.res

    def rowcol(self, lat, lon):
        """Indices of the cell containing (lat, lon); lon may be in any convention."""
        lat = np.asarray(lat, dtype=float)
        lon = (np.asarray(lon, dtype=float) - self.lon0) % 360.0
        r = np.clip(np.floor((90.0 - lat) * self.ppd).astype(int), 0, self.nrow - 1)
        c = np.clip(np.floor(lon * self.ppd).astype(int), 0, self.ncol - 1)
        return r, c

    def window(self, lat, lon, size_deg):
        """Row/col slices for a size_deg x size_deg box centred on (lat, lon). Returns (rows, cols) index arrays
        (cols wrap around the dateline)."""
        half = size_deg / 2.0
        r0, _ = self.rowcol(lat + half - 1e-9, lon)
        r1, _ = self.rowcol(lat - half + 1e-9, lon)
        n = int(round(size_deg * self.ppd))
        _, c0 = self.rowcol(lat, lon - half + 1e-9)
        cols = (c0 + np.arange(n)) % self.ncol
        rows = np.arange(r0, r1 + 1)
        return rows, cols


def resample_nearest(src: np.ndarray, src_grid: Grid, dst_grid: Grid) -> np.ndarray:
    """Nearest-neighbour resample between global lat/lon grids (handles different lon0 / resolution)."""
    lat = dst_grid.lats
    lon = dst_grid.lons
    r, _ = src_grid.rowcol(lat, np.zeros_like(lat))
    _, c = src_grid.rowcol(np.zeros_like(lon), lon)
    return src[np.ix_(r, c)]


def block_reduce_nan(a: np.ndarray, f: int, func=np.nanmean) -> np.ndarray:
    """Reduce a 2-D array by an integer factor f with a NaN-aware function."""
    ny, nx = a.shape[0] // f, a.shape[1] // f
    v = a[: ny * f, : nx * f].reshape(ny, f, nx, f)
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return func(v, axis=(1, 3))


def haversine_km(lat1, lon1, lat2, lon2, radius_km=EARTH_RADIUS_KM):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = p2 - p1
    dl = np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * radius_km * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
