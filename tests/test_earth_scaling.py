import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from features_earth import apparent_thermal_inertia, lst_qc_ok, mineral_indices, scale_albedo, scale_lst  # noqa: E402


def test_lst_scale_and_fill():
    raw = np.array([0, 15000, 16000], dtype=np.uint16)
    out = scale_lst(raw, {"scale_factor": 0.02, "add_offset": 0.0, "_FillValue": 0})
    assert np.isnan(out[0])
    assert out[1] == pytest.approx(300.0) and out[2] == pytest.approx(320.0)


def test_lst_qc_bits():
    # mandatory QA in bits 0-1, LST error in bits 6-7
    good = 0b00000000          # good quality, error <= 1 K
    other_2k = 0b01000001      # other quality, error <= 2 K
    err3k = 0b10000000         # error <= 3 K -> rejected
    cloud = 0b00000010         # not produced -> rejected
    ok = lst_qc_ok(np.array([good, other_2k, err3k, cloud], dtype=np.uint8))
    assert ok.tolist() == [True, True, False, False]


def test_albedo_scale_and_fill():
    raw = np.array([32767, 250, 1200], dtype=np.int16)
    a = scale_albedo(raw, {"scale_factor": 0.001, "_FillValue": 32767})
    assert np.isnan(a[0]) and a[1] == pytest.approx(0.25)
    assert np.isnan(a[2])  # > 1 is physically invalid


def test_ati():
    ati = apparent_thermal_inertia(np.array([0.2, 0.2]), np.array([320.0, 300.0]), np.array([280.0, 299.5]))
    assert ati[0] == pytest.approx(0.8 / 40.0)
    assert np.isnan(ati[1])  # dT below 1 K is undefined


def test_mineral_indices_flat_spectrum():
    one = np.full(3, 0.95, dtype=np.float32)
    qi, mi, c = mineral_indices(one, one, one, one, one)
    assert qi == pytest.approx(1.0) and mi == pytest.approx(1.0) and c == pytest.approx(0.0)
    # quartz reststrahlen: low emissivity in B10 & B12 relative to B11 -> QI > 1
    qi, _, _ = mineral_indices(one * 0.85, one * 0.92, one * 0.84, one, one)
    assert (qi > 1).all()


def test_accum_tile_placement():
    from features_earth import Accum
    from geo import Grid

    acc = Accum(Grid(20, lon0=-180.0), ["x"])
    tile = np.ones((100, 100), np.float32)
    tile[:5, :5] = 3.0  # NW 5x5 px -> first 0.05 deg cell
    acc.add_tile(-11.0, 154.0, "x", tile)
    m = acc.mean("x")
    r0, c0 = (90 + 11) * 20, (154 + 180) * 20
    assert m[r0, c0] == 3.0 and m[r0 + 1, c0 + 1] == 1.0
    assert np.isfinite(m[r0:r0 + 20, c0:c0 + 20]).all() and np.isnan(m[r0 - 1, c0])
    acc30 = Accum(Grid(30, lon0=-180.0), ["x"])
    acc30.add_tile(-11.0, 154.0, "x", np.ones((100, 100), np.float32))
    assert np.isfinite(acc30.mean("x")[(90 + 11) * 30:(90 + 12) * 30, (154 + 180) * 30:(155 + 180) * 30]).all()


def test_ati_monthly_filter():
    from features_earth import annual_ati

    alb = np.full((2, 2), 0.2, np.float32)
    # cell (0,0): 12 good months; (0,1): dT only 3 K every month (polar-style) -> NaN;
    # (1,0): only 3 usable months -> NaN; (1,1): frozen days (LST_day < 273) -> NaN
    day = np.zeros((12, 2, 2), np.float32)
    night = np.zeros((12, 2, 2), np.float32)
    day[:, 0, 0], night[:, 0, 0] = 300.0, 280.0
    day[:, 0, 1], night[:, 0, 1] = 280.0, 277.0
    day[:, 1, 0], night[:, 1, 0] = 280.0, 277.0
    day[:3, 1, 0], night[:3, 1, 0] = 300.0, 280.0
    day[:, 1, 1], night[:, 1, 1] = 270.0, 250.0
    ati, n = annual_ati(alb, day, night, min_dt=10.0, min_day_k=273.15, min_months=4)
    assert ati[0, 0] == pytest.approx(0.8 / 20.0) and n[0, 0] == 12
    assert np.isnan(ati[0, 1]) and n[0, 1] == 0
    assert np.isnan(ati[1, 0]) and n[1, 0] == 3
    assert np.isnan(ati[1, 1])


def test_permanent_snow_rule():
    from features_earth import permanent_snow

    s = np.full((12, 3), np.nan, np.float32)
    s[:, 0] = 90.0                    # snow all year -> permanent
    s[:, 1] = [90] * 6 + [0] * 6      # seasonal -> not permanent
    s[:2, 2] = 100.0                  # only 2 usable months (polar night elsewhere) -> not enough evidence
    assert permanent_snow(s, 50, 3).tolist() == [True, False, False]


def test_tir_contrast_and_nobs_mask_logic():
    bare = np.array([True, True, True, False])
    contrast = np.array([0.05, 0.01, 0.05, 0.05])
    nobs = np.array([5.0, 5.0, 1.0, 5.0])
    ok = bare & (contrast >= 0.02) & (nobs >= 3)
    assert ok.tolist() == [True, False, False, False]


def test_robustness_sampler():
    from robustness import sample_weights, top_regions

    rng = np.random.default_rng(1)
    ws = sample_weights({"a": 0.5, "b": 0.3, "c": 0.2}, rng, 200)
    arr = np.array([[w["a"], w["b"], w["c"]] for w in ws])
    assert np.allclose(arr.sum(1), 1.0) and (arr >= 0).all()
    assert abs(arr.mean(0)[0] - 0.5) < 0.05 and abs(arr.mean(0)[2] - 0.2) < 0.05
    ov = np.random.default_rng(2).random((40, 80)).astype(np.float32)
    ov[10:14, 10:14] = 1.5  # one solid block above every random value
    top = top_regions(ov, 0.01, 2, 1)
    assert top[10:14, 10:14].all() and top.sum() == 16
