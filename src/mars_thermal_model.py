"""Minimal 1-D surface thermal model for Mars (no atmosphere column; constant downwelling IR term).

Used ONLY for the 'thermal extremes' dimension because no global orbital surface-temperature map was supplied.
Inputs are the TES thermal inertia and albedo of the target window (data/processed/mars_fingerprints.json) and latitude.
Sanity check against Spanovich et al. (2006, Icarus 180, 314): Opportunity max daily surface temperature 295 K at Ls=339.
"""
from __future__ import annotations

import numpy as np

SIGMA = 5.670374e-8
SOL_S = 88775.244
S0_1AU = 1361.0
A_AU, ECC, LS_PERI = 1.52368, 0.0934, 251.0
OBLIQ = 25.19
RHO_C = 1.2e6  # volumetric heat capacity J m^-3 K^-1 (typical Mars regolith value used to convert TI -> k)
EMIS = 1.0
ATM_IR_FRAC = 0.04  # downwelling IR as a fraction of noon absorbed insolation (Kieffer et al. 1977 style approx.)


def sun_distance_au(ls_deg):
    return A_AU * (1 - ECC**2) / (1 + ECC * np.cos(np.radians(ls_deg - LS_PERI)))


def declination(ls_deg):
    return np.arcsin(np.sin(np.radians(OBLIQ)) * np.sin(np.radians(ls_deg)))


def diurnal_cycle(ti, albedo, lat_deg, ls_deg, n_sols=20, n_layers=30):
    """Return (Tmin, Tmax) of the converged surface temperature for one sol at season ls_deg."""
    k = ti**2 / RHO_C
    kappa = k / RHO_C
    skin = np.sqrt(kappa * SOL_S / np.pi)
    dz = skin / 8 * 1.15 ** np.arange(n_layers)  # geometric layer thickness, ~10 skin depths total
    z_c = np.cumsum(dz) - dz / 2
    dt = 0.4 * dz[0] ** 2 / kappa
    nstep = int(np.ceil(SOL_S / dt))
    dt = SOL_S / nstep
    flux = S0_1AU / sun_distance_au(ls_deg) ** 2
    dec = declination(ls_deg)
    lat = np.radians(lat_deg)
    cos_noon = max(np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec), 0.0)
    atm = ATM_IR_FRAC * (1 - albedo) * flux * cos_noon
    # start at radiative-equilibrium mean temperature
    h = np.linspace(0, 2 * np.pi, 200, endpoint=False)
    mean_ins = np.mean(np.clip(np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(h - np.pi), 0, None))
    T = np.full(n_layers, ((((1 - albedo) * flux * mean_ins) + atm) / (EMIS * SIGMA)) ** 0.25)
    dzi = np.diff(z_c)
    tmin, tmax = np.inf, -np.inf
    for sol in range(n_sols):
        last = sol == n_sols - 1
        for i in range(nstep):
            hour_angle = 2 * np.pi * (i * dt / SOL_S) - np.pi  # local noon at mid-sol
            cz = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(hour_angle)
            q_in = (1 - albedo) * flux * max(cz, 0.0) + atm
            fl = k * (T[1:] - T[:-1]) / dzi  # downward-positive conductive flux between layers (into layer above)
            dT = np.empty_like(T)
            dT[0] = (q_in - EMIS * SIGMA * T[0] ** 4 + fl[0]) / (RHO_C * dz[0])
            dT[1:-1] = (fl[1:] - fl[:-1]) / (RHO_C * dz[1:-1])
            dT[-1] = -fl[-1] / (RHO_C * dz[-1])  # insulated bottom
            T = T + dT * dt
            if last:
                tmin, tmax = min(tmin, T[0]), max(tmax, T[0])
    return float(tmin), float(tmax)


def extremes(ti, albedo, lat_deg, ls_list=tuple(range(0, 360, 30))):
    """Annual min / max surface temperature and mean diurnal range over the given seasons."""
    res = [diurnal_cycle(ti, albedo, lat_deg, ls) for ls in ls_list]
    mins, maxs = np.array(res).T
    return {"lst_min_k": float(mins.min()), "lst_max_k": float(maxs.max()),
            "diurnal_range_k": float(np.mean(maxs - mins))}
