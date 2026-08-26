#!/usr/bin/env python3
"""Offline smoke tests for public weather/soil source entry points.

These tests intentionally use tiny synthetic NetCDF/GeoTIFF/CSV fixtures and
mock live services. The goal is not to validate remote provider content; it is
to make sure every public source function can run its real entry point and emit
a DSSAT-shaped file without requiring API keys or large local datasets.
"""

from __future__ import annotations

import os
import sqlite3
import sys
import tempfile
from concurrent.futures import Future
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_REPO, "python"))


POINTS = pd.DataFrame({"ID": ["SRC1"], "LAT": [0.0], "LONG": [0.0]})
YEAR = 2001


class InlineExecutor:
    """Small ProcessPoolExecutor stand-in that executes submitted work inline."""

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def submit(self, fn, *args, **kwargs):
        fut = Future()
        try:
            fut.set_result(fn(*args, **kwargs))
        except Exception as exc:  # noqa: BLE001
            fut.set_exception(exc)
        return fut


def _dates(year: int = YEAR):
    return pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")


def _assert_wth(path: str | Path):
    path = Path(path)
    assert path.exists(), f"missing WTH: {path}"
    lines = [ln.rstrip("\n") for ln in path.read_text().splitlines() if ln.strip()]
    assert lines[0].startswith("$WEATHER DATA"), lines[:2]
    assert any(ln.lstrip().startswith("@  DATE") for ln in lines)
    assert any(ln[:4].isdigit() for ln in lines), "no daily weather rows"
    assert "nan" not in path.read_text().lower()


def _assert_sol(path: str | Path):
    path = Path(path)
    assert path.exists(), f"missing SOL: {path}"
    text = path.read_text(errors="replace")
    assert text.lstrip().startswith("*")
    assert "@  SLB" in text or "@SLB" in text
    assert "nan" not in text.lower()


def _native_nasa_frame(year: int = YEAR):
    dates = _dates(year)
    return pd.DataFrame({
        "YEAR": dates.year,
        "MM": dates.month,
        "DOY": dates.dayofyear,
        "DATE": [f"{d.year}{d.dayofyear:03d}" for d in dates],
        "T2M_MAX": np.full(len(dates), 28.0),
        "T2M_MIN": np.full(len(dates), 15.0),
        "ALLSKY_SFC_SW_DWN": np.full(len(dates), 18.0),
        "PRECTOTCORR": np.full(len(dates), 1.0),
        "T2MDEW": np.full(len(dates), 12.0),
        "RH2M": np.full(len(dates), 70.0),
        "WS2M": np.full(len(dates), 2.5),
    })


def _openmeteo_frame(year: int = YEAR):
    dates = _dates(year)
    return pd.DataFrame({
        "time": dates,
        "temperature_2m_max": np.full(len(dates), 28.0),
        "temperature_2m_min": np.full(len(dates), 15.0),
        "precipitation_sum": np.full(len(dates), 1.0),
        "shortwave_radiation_sum": np.full(len(dates), 18.0),
        "wind_speed_10m_mean": np.full(len(dates), 3.0),
        "dew_point_2m_mean": np.full(len(dates), 12.0),
        "relative_humidity_2m_mean": np.full(len(dates), 70.0),
        "YEAR": dates.year,
        "MM": dates.month,
        "DOY": dates.dayofyear,
        "DATE": [f"{d.year}{d.dayofyear:03d}" for d in dates],
    })


def _daymet_frame(year: int = YEAR):
    dates = _dates(year)
    return pd.DataFrame({
        "year": dates.year,
        "yday": dates.dayofyear,
        "tmax": np.full(len(dates), 28.0),
        "tmin": np.full(len(dates), 15.0),
        "prcp": np.full(len(dates), 1.0),
        "srad": np.full(len(dates), 250.0),
        "dayl": np.full(len(dates), 40000.0),
        "vp": np.full(len(dates), 1200.0),
    })


def _write_nc(path: Path, var: str, value: float, units: str = "", year: int = YEAR):
    xr = pytest.importorskip("xarray")
    dates = _dates(year)
    data = np.full((len(dates), 1, 1), value, dtype=float)
    da = xr.DataArray(
        data,
        coords={"time": dates, "lat": np.array([0.0]), "lon": np.array([0.0])},
        dims=["time", "lat", "lon"],
    )
    da.attrs["units"] = units
    ds = da.to_dataset(name=var)
    ds.to_netcdf(path)
    ds.close()


def _write_weather_set(nc_dir: Path, names: dict[str, tuple[str, float, str]], year: int = YEAR):
    nc_dir.mkdir(parents=True, exist_ok=True)
    for fname, (var, value, units) in names.items():
        _write_nc(nc_dir / fname, var, value, units, year)


def _write_raster(path: Path, value: float):
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=1,
        width=1,
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=from_origin(-0.5, 0.5, 1.0, 1.0),
        nodata=-9999,
    ) as dst:
        dst.write(np.array([[value]], dtype="float32"), 1)


def _write_soil_rasters(root: Path, depth_token: str = "one"):
    vals = {"sand": 40.0, "clay": 22.0, "silt": 38.0, "bulk": 130.0, "soc": 12.0}
    for prop, value in vals.items():
        _write_raster(root / f"{prop}_{depth_token}.tif", value)


def _profile_csv(path: Path):
    pd.DataFrame({
        "profile_id": ["P1", "P1"],
        "latitude": [0.0, 0.0],
        "longitude": [0.0, 0.0],
        "depth_bottom": [20, 60],
        "sand": [40.0, 38.0],
        "clay": [22.0, 24.0],
        "silt": [38.0, 38.0],
        "bdod": [1.30, 1.35],
        "soc_pct": [1.2, 0.8],
    }).to_csv(path, index=False)


LOCAL_NETCDF_WEATHER = [
    (
        "process_weather_chelsa_w5e5",
        "chelsa_nc_dir",
        {
            "tasmax_test.nc": ("tasmax", 300.0, "K"),
            "tasmin_test.nc": ("tasmin", 285.0, "K"),
            "pr_test.nc": ("pr", 2.0, "mm"),
            "rsds_test.nc": ("rsds", 200.0, "W m-2"),
        },
    ),
    (
        "process_weather_agmerra",
        "agmerra_nc_dir",
        {
            "tmax_test.nc": ("tmax", 300.0, "K"),
            "tmin_test.nc": ("tmin", 285.0, "K"),
            "prate_test.nc": ("prate", 2.0, "mm"),
            "srad_test.nc": ("srad", 200.0, "W m-2"),
        },
    ),
    (
        "process_weather_agcfsr",
        "agcfsr_nc_dir",
        {
            "tmax_test.nc": ("tmax", 300.0, "K"),
            "tmin_test.nc": ("tmin", 285.0, "K"),
            "prate_test.nc": ("prate", 2.0, "mm"),
            "srad_test.nc": ("srad", 200.0, "W m-2"),
        },
    ),
    (
        "process_weather_silo",
        "silo_nc_dir",
        {
            "tmax_test.nc": ("tmax", 27.0, "degC"),
            "tmin_test.nc": ("tmin", 14.0, "degC"),
            "rain_test.nc": ("rain", 2.0, "mm"),
            "srad_test.nc": ("srad", 18.0, "MJ m-2 day-1"),
            "vp_test.nc": ("vp", 12.27, "hPa"),
        },
    ),
    (
        "process_weather_mswx",
        "mswx_nc_dir",
        {
            "tmax_test.nc": ("tmax", 300.0, "K"),
            "tmin_test.nc": ("tmin", 285.0, "K"),
            "precip_test.nc": ("precip", 2.0, "mm"),
            "srad_test.nc": ("srad", 200.0, "W m-2"),
        },
    ),
    (
        "process_weather_crujra",
        "crujra_nc_dir",
        {
            "tmax_test.nc": ("tmax", 300.0, "K"),
            "tmin_test.nc": ("tmin", 285.0, "K"),
            "pr_test.nc": ("pr", 2.0, "mm"),
            "rsds_test.nc": ("rsds", 200.0, "W m-2"),
        },
    ),
    (
        "process_weather_pgf",
        "pgf_nc_dir",
        {
            "tmax_test.nc": ("tmax", 300.0, "K"),
            "tmin_test.nc": ("tmin", 285.0, "K"),
            "prcp_test.nc": ("prcp", 2.0, "mm"),
            "dswrf_test.nc": ("dswrf", 200.0, "W m-2"),
        },
    ),
    (
        "process_weather_merra2",
        "merra2_nc_dir",
        {
            "t2mmax_test.nc": ("t2mmax", 300.0, "K"),
            "t2mmin_test.nc": ("t2mmin", 285.0, "K"),
            "prectot_test.nc": ("prectot", 2.0, "mm"),
            "swgdn_test.nc": ("swgdn", 200.0, "W m-2"),
        },
    ),
]


@pytest.mark.parametrize("func_name,dir_arg,files", LOCAL_NETCDF_WEATHER)
def test_local_netcdf_weather_entrypoints(func_name, dir_arg, files, tmp_path):
    import dssatutils

    nc_dir = tmp_path / "nc"
    out_dir = tmp_path / "wth"
    _write_weather_set(nc_dir, files)
    getattr(dssatutils, func_name)(
        POINTS,
        YEAR,
        YEAR,
        str(out_dir),
        "ID",
        "LAT",
        "LONG",
        1,
        str(tmp_path / "weather.log"),
        **{dir_arg: str(nc_dir)},
    )
    _assert_wth(out_dir / "SRC1.WTH")


def test_anusplin_refuses_incomplete_standalone_forcing(tmp_path):
    """Core ANUSPLIN lacks solar radiation and must not emit a runnable WTH."""
    import dssatutils

    nc_dir = tmp_path / "nc"
    _write_weather_set(nc_dir, {
        "maxt_test.nc": ("maxt", 25.0, "degC"),
        "mint_test.nc": ("mint", 12.0, "degC"),
        "pcp_test.nc": ("pcp", 2.0, "mm"),
    })
    with pytest.raises(FileNotFoundError, match="requires SRAD"):
        dssatutils.process_weather_anusplin(
            POINTS, YEAR, YEAR, str(tmp_path / "wth"), "ID", "LAT", "LONG", 1,
            str(tmp_path / "weather.log"), anusplin_nc_dir=str(nc_dir),
        )


def test_process_weather_eobs_entrypoint(tmp_path):
    from dssatutils import process_weather_eobs

    nc_dir = tmp_path / "eobs"
    _write_weather_set(nc_dir, {
        "tx_test.nc": ("tx", 25.0, "degC"),
        "tn_test.nc": ("tn", 12.0, "degC"),
        "rr_test.nc": ("rr", 2.0, "mm"),
        "qq_test.nc": ("qq", 200.0, "W m-2"),
        "tg_test.nc": ("tg", 18.0, "degC"),
        "hu_test.nc": ("hu", 70.0, "%"),
        "fg_test.nc": ("fg", 3.0, "m s-1"),
    })
    out_dir = tmp_path / "wth"
    process_weather_eobs(
        POINTS, YEAR, YEAR, str(out_dir), "ID", "LAT", "LONG", 1,
        str(tmp_path / "eobs.log"), eobs_nc_dir=str(nc_dir),
    )
    _assert_wth(out_dir / "SRC1.WTH")


def test_process_weather_dwd_entrypoint(tmp_path):
    from dssatutils import process_weather_dwd

    dates = _dates()
    daily = pd.DataFrame({
        "DATE": dates,
        "TXK": 25.0,
        "TNK": 12.0,
        "RSK": 2.0,
        "UPM": 70.0,
        "FM": 2.0,
        "SDK": 8.0,
        "VPM": 12.27,
    })
    stations = pd.DataFrame({
        "station_id": ["00001"],
        "von": [20000101],
        "bis": [20011231],
        "elev": [100.0],
        "lat": [0.0],
        "lon": [0.0],
    })
    with patch("dssatutils.weather_dwd._dwd_stations", return_value=stations), \
            patch("dssatutils.weather_dwd._historical_index", return_value={"00001": "dummy.zip"}), \
            patch("dssatutils.weather_dwd._fetch_station", return_value=daily):
        out_dir = tmp_path / "wth"
        process_weather_dwd(
            POINTS, YEAR, YEAR, str(out_dir), "ID", "LAT", "LONG", 1,
            str(tmp_path / "dwd.log"), str(tmp_path / "dwd_cache"),
        )
    _assert_wth(out_dir / "SRC1.WTH")


def test_process_weather_ghcn_entrypoint(tmp_path):
    from dssatutils import process_weather_ghcn

    dates = _dates()
    frame = pd.DataFrame({
        "DATE": [f"{d.year}{d.dayofyear:03d}" for d in dates],
        "YEAR": dates.year,
        "MM": dates.month,
        "SRAD": -99.0,
        "TMAX": 25.0,
        "TMIN": 12.0,
        "RAIN": 2.0,
        "TDEW": -99.0,
        "RH2M": -99.0,
        "WIND": -99.0,
    })
    stations = pd.DataFrame({"sid": ["GHCN0001"], "lat": [0.0], "lon": [0.0]})
    with patch("dssatutils.weather_ghcn._load_stations", return_value=stations), \
            patch("dssatutils.weather_ghcn._fetch_station", return_value=frame):
        out_dir = tmp_path / "wth"
        process_weather_ghcn(
            POINTS, YEAR, YEAR, str(out_dir), "ID", "LAT", "LONG", 1,
            str(tmp_path / "ghcn.log"), str(tmp_path / "ghcn_cache"),
        )
    _assert_wth(out_dir / "SRC1.WTH")


def test_process_weather_prism_entrypoint(tmp_path):
    from dssatutils import process_weather_prism

    vals = {"ppt": 2.0, "tmax": 25.0, "tmin": 12.0, "tdmean": 10.0}
    with patch("dssatutils.weather_prism._download_grid", side_effect=lambda var, day, cache: var), \
            patch("dssatutils.weather_prism._sample_raster",
                  side_effect=lambda path, lats, lons: np.full(len(lats), vals[path])):
        out_dir = tmp_path / "wth"
        process_weather_prism(
            POINTS, YEAR, YEAR, str(out_dir), "ID", "LAT", "LONG", 1,
            str(tmp_path / "prism.log"), str(tmp_path / "prism_cache"),
        )
    _assert_wth(out_dir / "SRC1.WTH")


HYBRID_RAIN_WEATHER = [
    ("dssatutils.weather_aphrodite", "process_weather_aphrodite", "aphrodite_nc_dir", "aphrodite_pr.nc", "precip"),
    ("dssatutils.weather_tamsat", "process_weather_tamsat", "tamsat_nc_dir", "tamsat_rfe.nc", "rfe"),
    ("dssatutils.weather_mswep", "process_weather_mswep", "mswep_nc_dir", "mswep_precip.nc", "precip"),
]


@pytest.mark.parametrize("module_name,func_name,dir_arg,fname,var", HYBRID_RAIN_WEATHER)
def test_hybrid_rain_weather_entrypoints(module_name, func_name, dir_arg, fname, var, tmp_path):
    module = __import__(module_name, fromlist=[func_name])
    nc_dir = tmp_path / "nc"
    _write_weather_set(nc_dir, {fname: (var, 7.0, "mm")})
    out_dir = tmp_path / "wth"
    with patch(f"{module_name}._fetch_nasa_power", return_value=_native_nasa_frame()), \
            patch(f"{module_name}.ProcessPoolExecutor", InlineExecutor):
        getattr(module, func_name)(
            POINTS, YEAR, YEAR, str(out_dir), "ID", "LAT", "LONG", 1,
            str(tmp_path / "hybrid.log"), **{dir_arg: str(nc_dir)},
        )
    _assert_wth(out_dir / "SRC1.WTH")
    assert "   7.0" in (out_dir / "SRC1.WTH").read_text()


def test_process_weather_openmeteo_entrypoint(tmp_path):
    from dssatutils import process_weather_openmeteo

    with patch("dssatutils.weather_openmeteo._fetch_open_meteo", return_value=_openmeteo_frame()), \
            patch("dssatutils.weather_openmeteo.ProcessPoolExecutor", InlineExecutor):
        out_dir = tmp_path / "wth"
        process_weather_openmeteo(
            POINTS, YEAR, YEAR, str(out_dir), "ID", "LAT", "LONG", 1,
            str(tmp_path / "openmeteo.log"),
        )
    _assert_wth(out_dir / "SRC1.WTH")
    assert "  12.0  70.0" in (out_dir / "SRC1.WTH").read_text()


def test_process_weather_nasapower_entrypoint(tmp_path):
    from dssatutils import process_weather_nasapower

    with patch("dssatutils.weather_nasapower._fetch_nasa_power", return_value=_native_nasa_frame()), \
            patch("dssatutils.weather_nasapower.ProcessPoolExecutor", InlineExecutor):
        out_dir = tmp_path / "wth"
        process_weather_nasapower(
            POINTS, YEAR, YEAR, str(out_dir), "ID", "LAT", "LONG", 1,
            str(tmp_path / "nasa.log"),
        )
    _assert_wth(out_dir / "SRC1.WTH")


def test_process_weather_daymet_entrypoint(tmp_path):
    from dssatutils import process_weather_daymet

    with patch("dssatutils.weather_daymet._download_daymet", return_value=_daymet_frame()), \
            patch("dssatutils.weather_daymet.ProcessPoolExecutor", InlineExecutor):
        out_dir = tmp_path / "wth"
        process_weather_daymet(
            POINTS, YEAR, YEAR, str(out_dir), "ID", "LAT", "LONG", 1,
            str(tmp_path / "daymet.log"),
        )
    _assert_wth(out_dir / "SRC1.WTH")


RASTER_SOIL_SOURCES = [
    ("process_soils_slga", "slga_raster_dir"),
    ("process_soils_wise30sec", "wise30sec_raster_dir"),
    ("process_soils_gsde", "gsde_raster_dir"),
    ("process_soils_china", "china_raster_dir"),
    ("process_soils_slc", "slc_raster_dir"),
    ("process_soils_esdb", "esdb_raster_dir"),
]


@pytest.mark.parametrize("func_name,dir_arg", RASTER_SOIL_SOURCES)
def test_raster_soil_entrypoints(func_name, dir_arg, tmp_path):
    import dssatutils

    raster_dir = tmp_path / "rasters"
    _write_soil_rasters(raster_dir, "one")
    out_dir = tmp_path / "sol"
    getattr(dssatutils, func_name)(
        POINTS,
        str(raster_dir),
        str(tmp_path / "soil.csv"),
        str(out_dir),
        id_col="ID",
        lat_col="LAT",
        long_col="LONG",
        depth_specs=[("one", 30, 15.0)],
    )
    _assert_sol(out_dir / "SRC1.SOL")


def test_process_soils_hihydrosoil_entrypoint(tmp_path):
    from dssatutils import process_soils_hihydrosoil
    from dssatutils.soil_raster_common import DEPTHS

    raster_dir = tmp_path / "hihydro"
    for depth, _, _ in DEPTHS:
        for prop, value in {
            "pf42": 0.12,
            "pf25": 0.28,
            "thetas": 0.46,
            "ksat": 24.0,
            "organic": 2.0,
            "sand": 40.0,
            "clay": 22.0,
            "silt": 38.0,
        }.items():
            _write_raster(raster_dir / f"{prop}_{depth}.tif", value)
    out_dir = tmp_path / "sol"
    process_soils_hihydrosoil(
        POINTS,
        str(raster_dir),
        str(tmp_path / "soil.csv"),
        str(out_dir),
        id_col="ID",
        lat_col="LAT",
        long_col="LONG",
        integer_scale=1.0,
    )
    _assert_sol(out_dir / "SRC1.SOL")


PROFILE_SOIL_SOURCES = [
    ("process_soils_wosis", "wosis_profile_csv"),
    ("process_soils_febr", "febr_profile_csv"),
]


@pytest.mark.parametrize("func_name,csv_arg", PROFILE_SOIL_SOURCES)
def test_profile_csv_soil_entrypoints(func_name, csv_arg, tmp_path):
    import dssatutils

    src = tmp_path / "profiles.csv"
    _profile_csv(src)
    out_dir = tmp_path / "sol"
    getattr(dssatutils, func_name)(
        POINTS,
        str(src),
        str(tmp_path / "soil.csv"),
        str(out_dir),
        id_col="ID",
        lat_col="LAT",
        long_col="LONG",
    )
    _assert_sol(out_dir / "SRC1.SOL")


def test_process_soils_lucas_entrypoint(tmp_path):
    from dssatutils import process_soils_lucas

    lucas = tmp_path / "lucas.csv"
    pd.DataFrame({
        "POINTID": ["L1"],
        "TH_LAT": [0.0],
        "TH_LONG": [0.0],
        "clay": [22.0],
        "sand": [40.0],
        "silt": [38.0],
        "OC": [12.0],
    }).to_csv(lucas, index=False)
    out_dir = tmp_path / "sol"
    assert process_soils_lucas(
        POINTS,
        str(tmp_path / "soil.csv"),
        str(out_dir),
        1,
        "ID",
        "LAT",
        "LONG",
        lucas_csv=str(lucas),
        max_dist_km=1.0,
    )
    _assert_sol(out_dir / "SRC1.SOL")


def test_process_soils_isdasoil_entrypoint(tmp_path):
    from dssatutils import process_soils_isdasoil

    props = {
        "clay_content": np.array([[22.0, 24.0]]),
        "sand_content": np.array([[40.0, 38.0]]),
        "carbon_organic": np.array([[12.0, 8.0]]),
        "bulk_density": np.array([[1.30, 1.35]]),
    }
    with patch("dssatutils.soil_isdasoil._sample_property",
               side_effect=lambda prop, xs, ys: props[prop]):
        out_dir = tmp_path / "sol"
        assert process_soils_isdasoil(
            POINTS,
            str(tmp_path / "soil.csv"),
            str(out_dir),
            1,
            "ID",
            "LAT",
            "LONG",
        )
    _assert_sol(out_dir / "SRC1.SOL")


def test_process_soils_openlandmap_entrypoint(tmp_path):
    from dssatutils import process_soils_openlandmap

    values = iter([22.0, 40.0, 38.0, 130.0, 60.0] * 3)

    def fake_sample(url, lats, lons):
        return np.full(len(lats), next(values), dtype=float)

    assets = [(15.0, "asset")]
    with patch("dssatutils.soil_openlandmap._resolve_assets", return_value=assets), \
            patch("dssatutils.soil_openlandmap._nearest_asset",
                  side_effect=lambda resolved, target_mid: resolved[0][1]), \
            patch("dssatutils.soil_openlandmap._sample_cog", side_effect=fake_sample):
        out_dir = tmp_path / "sol"
        process_soils_openlandmap(
            POINTS,
            str(tmp_path / "soil.csv"),
            str(out_dir),
            id_col="ID",
            lat_col="LAT",
            long_col="LONG",
        )
    _assert_sol(out_dir / "SRC1.SOL")


def test_process_soils_agmip_entrypoint(tmp_path):
    from dssatutils import process_soils_agmip

    master = """*AGMIP001 AgMIP test profile
@SITE        COUNTRY          LAT     LONG SCS FAMILY
 SITEA       Test           0.0000   0.0000 -99
@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI
    15   -99 0.100 0.250 0.430 1.00  -99 1.30 1.20 22.0 38.0
"""
    source = tmp_path / "agmip.sol"
    source.write_text(master)
    out_dir = tmp_path / "sol"
    process_soils_agmip(POINTS, str(source), str(tmp_path / "soil.csv"), str(out_dir))
    _assert_sol(out_dir / "00000001.SOL")


def test_process_soils_hwsd_entrypoint(tmp_path):
    from dssatutils import process_soils_hwsd

    raster = tmp_path / "hwsd.tif"
    raster.write_text("placeholder")
    db = tmp_path / "hwsd.sqlite"
    conn = sqlite3.connect(db)
    conn.execute(
        "CREATE TABLE HWSD2_LAYERS ("
        "HWSD2_SMU_ID INTEGER, SEQUENCE INTEGER, SHARE REAL, "
        "TOPDEP REAL, BOTDEP REAL, SAND REAL, CLAY REAL, SILT REAL, "
        "BULK REAL, ORG_CARBON REAL, COARSE REAL)"
    )
    conn.execute(
        "INSERT INTO HWSD2_LAYERS VALUES "
        "(1, 1, 100, 0, 30, 40.0, 22.0, 38.0, 1.30, 1.2, 0.0)"
    )
    conn.commit()
    conn.close()
    with patch("dssatutils.soil_hwsd._sample_smu_ids", return_value=np.array([1])):
        out_dir = tmp_path / "sol"
        process_soils_hwsd(
            POINTS,
            str(raster),
            str(db),
            str(tmp_path / "soil.csv"),
            str(out_dir),
            id_col="ID",
            lat_col="LAT",
            long_col="LONG",
        )
    _assert_sol(out_dir / "SRC1.SOL")


def test_all_public_source_entrypoints_are_accounted_for():
    import dssatutils

    accounted = {
        # Original/live-or-mocked sources covered here or in test_comprehensive.py.
        "process_weather_daymet",
        "process_weather_gridmet",
        "process_weather_nasapower",
        "process_weather_openmeteo",
        "process_weather_agera5",
        "process_weather_nasapower_chirps",
        "process_weather_nasapower_chirps_v3",
        "process_weather_era5_land",
        "process_soils_soilgrids",
        "process_soils_soilgrids_online",
        "process_soils_ssurgo",
        "process_soils_ssurgo_alderman",
        "process_soils_polaris",
        # Newer/local or regional source entry points.
        "process_weather_cmfd",
        "process_weather_dwd",
        "process_weather_eobs",
        "process_weather_xavier",
        "process_weather_chelsa_w5e5",
        "process_weather_agmerra",
        "process_weather_agcfsr",
        "process_weather_silo",
        "process_weather_prism",
        "process_weather_mswx",
        "process_weather_mswep",
        "process_weather_crujra",
        "process_weather_terraclimate",
        "process_weather_aphrodite",
        "process_weather_anusplin",
        "process_weather_tamsat",
        "process_weather_ghcn",
        "process_weather_pgf",
        "process_weather_merra2",
        "process_soils_gnatsgo",
        "process_soils_isdasoil",
        "process_soils_lucas",
        "process_soils_hwsd",
        "process_soils_agmip",
        "process_soils_hihydrosoil",
        "process_soils_slga",
        "process_soils_wise30sec",
        "process_soils_wosis",
        "process_soils_gsde",
        "process_soils_china",
        "process_soils_febr",
        "process_soils_slc",
        "process_soils_esdb",
        "process_soils_openlandmap",
    }
    public_sources = {n for n in dssatutils.__all__ if n.startswith(("process_weather_", "process_soils_"))}
    assert public_sources == accounted
