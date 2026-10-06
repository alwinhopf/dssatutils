# Tests for shared solar radiation calculations and Bristow-Campbell estimation.
import numpy as np
import pandas as pd
import pytest
from unittest.mock import patch
from pathlib import Path

from dssatutils.weather_solar import (
    extraterrestrial_radiation,
    estimate_srad_bristow_campbell,
)
from dssatutils.weather_prism import process_weather_prism


def test_extraterrestrial_radiation_fao56_example8():
    """FAO-56 Example 8: Latitude 20 deg S (-20.0), DOY 246 (3 September).

    FAO-56 worked calculation gives Ra = 32.2 MJ/m2/day.
    """
    ra, ws = extraterrestrial_radiation(-20.0, 246)
    assert np.isclose(ra, 32.2, atol=0.1)
    assert np.isclose(ra, 32.19, atol=0.01)
    assert 0.0 < ws < np.pi


def test_extraterrestrial_radiation_polar_limits():
    """Sunset hour angle argument must clip to [-1, 1] without domain error."""
    # North Pole at winter solstice (DOY 355): 24h darkness, Ra = 0
    ra_np_win, ws_np_win = extraterrestrial_radiation(90.0, 355)
    assert ra_np_win == 0.0
    assert ws_np_win == 0.0

    # North Pole at summer solstice (DOY 172): 24h sunlight, ws = pi, Ra > 0
    ra_np_sum, ws_np_sum = extraterrestrial_radiation(90.0, 172)
    assert ra_np_sum > 0.0
    assert np.isclose(ws_np_sum, np.pi)

    # South Pole at winter solstice (DOY 172): 24h darkness, Ra = 0
    ra_sp_win, ws_sp_win = extraterrestrial_radiation(-90.0, 172)
    assert ra_sp_win == 0.0
    assert ws_sp_win == 0.0

    # South Pole at summer solstice (DOY 355): 24h sunlight, ws = pi, Ra > 0
    ra_sp_sum, ws_sp_sum = extraterrestrial_radiation(-90.0, 355)
    assert ra_sp_sum > 0.0
    assert np.isclose(ws_sp_sum, np.pi)


def test_bristow_campbell_physical_bounds():
    """Rs must always satisfy 0 <= Rs <= 0.80 * Ra."""
    dates = pd.date_range("2020-01-01", "2020-12-31", freq="D")
    n = len(dates)
    # Realistic seasonal temperatures with diurnal range 4 to 20 C
    tmax = 20.0 + 10.0 * np.sin(2.0 * np.pi * dates.dayofyear / 365.0)
    tmin = tmax - np.linspace(4.0, 20.0, n)

    for lat in [-45.0, -15.0, 0.0, 30.0, 50.0]:
        rs = estimate_srad_bristow_campbell(dates, tmax, tmin, lat)
        ra, _ = extraterrestrial_radiation(lat, dates.dayofyear.values)

        assert np.all(np.isfinite(rs))
        assert np.all(rs >= 0.0)
        assert np.all(rs <= 0.80 * ra + 1e-6)


def test_bristow_campbell_diurnal_range_monotonicity():
    """Higher diurnal range yields higher solar radiation up to clear-sky limit."""
    dates = pd.date_range("2020-07-01", periods=10, freq="D")
    lat = 35.0

    # Diurnal range 2 C vs 18 C
    rs_low = estimate_srad_bristow_campbell(dates, np.full(10, 22.0), np.full(10, 20.0), lat)
    rs_high = estimate_srad_bristow_campbell(dates, np.full(10, 35.0), np.full(10, 17.0), lat)

    assert np.all(rs_high > rs_low)


def test_bristow_campbell_gap_and_missing_handling():
    """Check handling of date gaps, terminal day, and NaN temperatures."""
    dates = pd.to_datetime(["2020-01-01", "2020-01-02", "2020-01-10", "2020-01-11"])
    tmax = [25.0, np.nan, 22.0, 24.0]
    tmin = [10.0, 12.0, 10.0, 11.0]

    rs = estimate_srad_bristow_campbell(dates, tmax, tmin, 30.0)
    assert len(rs) == 4
    assert np.isfinite(rs[0]) and rs[0] > 0
    assert np.isnan(rs[1])  # tmax is NaN
    assert np.isfinite(rs[2]) and rs[2] > 0  # after a gap, still works
    assert np.isfinite(rs[3]) and rs[3] > 0  # terminal day works


def test_bristow_campbell_cross_language_numerical_parity():
    """Verify numerical parity with R implementation on synthetic sequence."""
    dates = pd.date_range("2020-01-01", periods=5, freq="D")
    tmax = [15.0, 18.0, np.nan, 20.0, 22.0]
    tmin = [5.0, 7.0, 6.0, np.nan, 10.0]
    rs = estimate_srad_bristow_campbell(dates, tmax, tmin, 30.0)

    assert np.isclose(rs[0], 11.00180729, atol=1e-5)
    assert np.isclose(rs[1], 13.64421333, atol=1e-5)
    assert np.isnan(rs[2])
    assert np.isnan(rs[3])
    assert np.isclose(rs[4], 14.08366375, atol=1e-5)


def test_process_weather_prism_srad_method(tmp_path):
    """Test process_weather_prism with bristow_campbell vs none vs invalid."""
    import geopandas as gpd
    from shapely.geometry import Point

    points_gdf = gpd.GeoDataFrame({
        "ID": ["P1"],
        "LAT": [32.5],
        "LONG": [-85.0],
        "geometry": [Point(-85.0, 32.5)],
    }, crs="EPSG:4326")

    vals = {"ppt": 0.0, "tmax": 28.0, "tmin": 14.0, "tdmean": 12.0}

    with patch("dssatutils.weather_prism._download_grid", side_effect=lambda var, day, cache: var), \
            patch("dssatutils.weather_prism._sample_raster",
                  side_effect=lambda path, lats, lons: np.full(len(lats), vals[path])):

        # 1. Default (bristow_campbell)
        out_bc = tmp_path / "wth_bc"
        process_weather_prism(
            points_gdf, 2020, 2020, str(out_bc), "ID", "LAT", "LONG", 1,
            str(tmp_path / "prism_bc.log"), str(tmp_path / "cache"),
            srad_method="bristow_campbell",
        )
        wth_bc_text = (out_bc / "P1.WTH").read_text()
        assert "SRAD estimated: Bristow-Campbell 1984" in wth_bc_text
        # SRAD should not be -99
        daily_lines = [ln for ln in wth_bc_text.splitlines() if ln.strip() and ln.strip()[:4] == "2020"]
        assert len(daily_lines) > 0
        srad_val = float(daily_lines[0].split()[1])
        assert 10.0 < srad_val < 35.0

        # 2. srad_method='none'
        out_none = tmp_path / "wth_none"
        process_weather_prism(
            points_gdf, 2020, 2020, str(out_none), "ID", "LAT", "LONG", 1,
            str(tmp_path / "prism_none.log"), str(tmp_path / "cache"),
            srad_method="none",
        )
        wth_none_text = (out_none / "P1.WTH").read_text()
        assert "SRAD estimated" not in wth_none_text
        daily_lines_none = [ln for ln in wth_none_text.splitlines() if ln.strip() and ln.strip()[:4] == "2020"]
        assert len(daily_lines_none) > 0
        srad_none_val = float(daily_lines_none[0].split()[1])
        assert srad_none_val == -99.0

        # 3. Invalid srad_method raises ValueError
        with pytest.raises(ValueError, match="Unknown srad_method"):
            process_weather_prism(
                points_gdf, 2020, 2020, str(tmp_path / "wth_inv"), "ID", "LAT", "LONG", 1,
                str(tmp_path / "prism_inv.log"), str(tmp_path / "cache"),
                srad_method="invalid_option",
            )
