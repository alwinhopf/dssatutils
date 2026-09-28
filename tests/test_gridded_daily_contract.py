"""Do not silently aggregate ambiguous NetCDF time intervals."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from dssatutils.weather_gridded_common import convert_units, extract_netcdf_series


def dataset(times, values, units="mm", lat=40.0, lon=10.0):
    ds = xr.Dataset({"rain": (("time", "lat", "lon"), np.array(values).reshape(-1, 1, 1))},
                    coords={"time": pd.to_datetime(times), "lat": [lat], "lon": [lon]})
    ds.rain.attrs["units"] = units
    return ds


@pytest.mark.parametrize("times,values,units,error", [
    (["2001-01-01", "2001-01-01"], [10., 10.], "mm", "Duplicate weather timestamps"),
    (pd.date_range("2001-01-01", periods=24, freq="h"), [1 / 86400] * 24,
     "kg m-2 s-1", "Subdaily weather data"),
    (["2001-01-01 00:00", "2001-01-01 12:00"], [1., 2.], "mm", "Subdaily weather data"),
])
def test_ambiguous_time_records_fail(monkeypatch, times, values, units, error):
    monkeypatch.setattr(xr, "open_dataset", lambda p: dataset(times, values, units))
    with pytest.raises(ValueError, match=error):
        extract_netcdf_series("unused", ["rain"], ["p"], [40], [10], 2001, 2001, "rain")


def test_overlapping_files_fail(monkeypatch):
    monkeypatch.setattr(xr, "open_dataset", lambda p: dataset(["2001-01-01"], [10.]))
    with pytest.raises(ValueError, match="Duplicate weather timestamps"):
        extract_netcdf_series(["a", "b"], ["rain"], ["p"], [40], [10], 2001, 2001, "rain")


def test_daily_rate_conversion_unchanged(monkeypatch):
    monkeypatch.setattr(xr, "open_dataset", lambda p: dataset(
        ["2001-01-01", "2001-01-02"], [1 / 86400, 2 / 86400], "kg m-2 s-1"))
    result = extract_netcdf_series("unused", ["rain"], ["p"], [40], [10], 2001, 2001, "rain")
    assert result["p"] == pytest.approx({"2001001": 1., "2001002": 2.})


def test_bidirectional_longitude_wrapping(monkeypatch):
    # Case 1: 0..360 grid with negative query longitude (-80 mapped to 280)
    monkeypatch.setattr(xr, "open_dataset", lambda p: dataset(
        ["2001-01-01"], [5.0], "mm", lat=25.0, lon=280.0))
    res = extract_netcdf_series("unused", ["rain"], ["p1"], [25.0], [-80.0], 2001, 2001, "rain")
    assert res["p1"] == pytest.approx({"2001001": 5.0})

    # Case 2: -180..180 grid with query longitude > 180 (280 mapped to -80)
    monkeypatch.setattr(xr, "open_dataset", lambda p: dataset(
        ["2001-01-01"], [7.0], "mm", lat=25.0, lon=-80.0))
    res2 = extract_netcdf_series("unused", ["rain"], ["p2"], [25.0], [280.0], 2001, 2001, "rain")
    assert res2["p2"] == pytest.approx({"2001001": 7.0})


def test_convert_units_edge_cases():
    # A magnitude guess cannot establish a unit contract.
    with pytest.raises(ValueError, match="units"):
        convert_units([25.0], "unknown", "temp")

    # Temp: 'degK', 'kelvin', 'K' should trigger subtraction
    assert convert_units([298.15], "degK", "temp")[0] == pytest.approx(25.0)
    assert convert_units([298.15], "kelvin", "temp")[0] == pytest.approx(25.0)

    with pytest.raises(ValueError, match="units"):
        convert_units([15_000_000.0], "", "srad")

    # SRAD: 'MJ' units preserved
    assert convert_units([18.5], "MJ m-2 day-1", "srad")[0] == pytest.approx(18.5)
