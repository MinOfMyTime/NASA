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
