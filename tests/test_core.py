import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from geo import Grid, lon_e180_to_e360, lon_e360_to_e180, lon_w360_to_e360, resample_nearest  # noqa: E402
from io_mars import flip_ti_orientation  # noqa: E402
from similarity import dimension_score, feature_distance, overall_score  # noqa: E402
from stats import quantile_table, to_percentile  # noqa: E402
from topo import roughness_detrended  # noqa: E402


def test_lon_conversions():
    assert lon_e360_to_e180(354.5) == pytest.approx(-5.5)
    assert lon_e360_to_e180(77.6) == pytest.approx(77.6)
    assert lon_e360_to_e180(180.0) == pytest.approx(-180.0)
    assert lon_e180_to_e360(-69.25) == pytest.approx(290.75)
    assert lon_w360_to_e360(5.5) == pytest.approx(354.5)
    assert lon_w360_to_e360(0.0) == pytest.approx(0.0)


def test_grid_rowcol_and_wrap():
    g = Grid(4, lon0=0.0)
    r, c = g.rowcol(89.9, 0.1)
    assert (r, c) == (0, 0)
    r, c = g.rowcol(-89.9, -0.1)  # -0.1 E == 359.9 E -> last column
    assert (r, c) == (719, 1439)
    rows, cols = g.window(0.0, 0.0, 1.0)
    assert len(rows) == 4 and len(cols) == 4 and cols[0] == 1438  # wraps across 0E


def test_resample_between_lon_conventions():
    src = Grid(1, lon0=-180.0)
    a = np.tile(np.arange(360, dtype=float), (180, 1))  # value = column index in -180..180 grid
    dst = Grid(1, lon0=0.0)
    b = resample_nearest(a, src, dst)
    # dst column 0 is 0..1E -> src column 180
    assert b[0, 0] == 180 and b[0, 359] == 179


def test_ti_flip_orientation():
    raw = np.zeros((3600, 7200), np.int16)
    raw[0, 0] = 1  # south-most row, W lon 0..0.05 -> E lon 359.95..360 (last col), north-up -> last row
    out = flip_ti_orientation(raw)
    assert out[-1, -1] == 1


def test_percentile_transform():
    t = quantile_table(np.arange(1000, dtype=float))
    p = to_percentile(np.array([0.0, 499.5, 999.0, np.nan]), t)
    assert p[0] == pytest.approx(0, abs=1e-3) and p[1] == pytest.approx(0.5, abs=2e-3) and p[2] == pytest.approx(1)
    assert np.isnan(p[3])
    ties = quantile_table(np.r_[np.zeros(500), np.arange(1, 501.0)])
    assert to_percentile(np.array([0.0]), ties)[0] == pytest.approx(0.25, abs=0.01)  # mid-rank for ties


def test_roughness_removes_plane():
    y, x = np.mgrid[0:50, 0:50].astype(float)
    plane = 3 * x + 7 * y
    r = roughness_detrended(plane, 5)
    assert np.nanmax(r[5:-5, 5:-5]) < 1e-6


def test_similarity_function():
    fp = {"features": {"pct_ti": {"median": 0.8, "iqr": 0.2}, "albedo": {"median": 0.2, "iqr": 0.0}}}
    cfg = {"similarity": {"sigma": 1.0, "min_scale": {"pct": 0.05, "albedo": 0.02}}}
    assert feature_distance(1.0, fp["features"]["pct_ti"], 0.05) == pytest.approx(1.0)
    feats = {"pct_ti": np.array([0.8, 1.0]), "albedo": np.array([0.2, 0.2])}
    s = dimension_score(feats, fp, "thermo", cfg)
    assert s[0] == pytest.approx(1.0)
    assert s[1] == pytest.approx(np.exp(-0.5 / 2))  # d^2 = (1^2 + 0^2)/2
    assert ((s >= 0) & (s <= 1)).all()


def test_overall_requires_all_weighted_dimensions():
    scores = {"a": np.array([1.0, 1.0]), "b": np.array([0.0, np.nan])}
    o = overall_score(scores, {"a": 0.5, "b": 0.5})
    assert o[0] == pytest.approx(0.5) and np.isnan(o[1])  # missing dimension -> n/a, not a free pass
    o = overall_score(scores, {"a": 0.5, "b": 0.5}, require_all=False)
    assert o[1] == pytest.approx(1.0)
    o = overall_score(scores, {"a": 1.0, "b": 0.0})  # zero-weight dimension may be missing
    assert o[1] == pytest.approx(1.0)
