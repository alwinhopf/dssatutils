#!/usr/bin/env python3
# File: tests/test_new_sources.py
# ---------------------------------------------------------------------------
# OFFLINE tests for the gNATSGO (soil), DWD and E-OBS (weather) sources.
# No network and no DSSAT install required: the network-bound lookups are
# isolated from the physics/formatting helpers, which are exercised here with
# deterministic synthetic inputs. A cross-language section asserts the R twins
# carry the same API + algorithms so the two implementations can't drift.
#
# Run:  python -m pytest tests/test_new_sources.py     (or python tests/test_new_sources.py)
# ---------------------------------------------------------------------------

import math
import os
import sys
import tempfile

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
_WORKSPACE = os.path.dirname(_REPO)
sys.path.insert(0, os.path.join(_REPO, "python"))


# --- imports ----------------------------------------------------------------

def test_new_modules_import_and_export():
    import dssatutils
    for name in ("process_soils_gnatsgo", "process_weather_dwd", "process_weather_eobs",
                 "process_soils_isdasoil", "process_soils_lucas",
                 "process_weather_xavier", "process_weather_cmfd",
                 "process_weather_era5_land", "process_soils_agmip",
                 "process_weather_chelsa_w5e5", "process_weather_agmerra",
                 "process_weather_agcfsr", "process_weather_silo",
                 "process_weather_prism", "process_soils_hihydrosoil",
                 "process_soils_slga", "process_weather_mswx",
                 "process_weather_mswep", "process_weather_crujra",
                 "process_weather_terraclimate", "process_soils_wise30sec",
                 "process_soils_wosis",
                 "process_weather_aphrodite", "process_weather_anusplin",
                 "process_weather_tamsat", "process_weather_ghcn",
                 "process_weather_pgf", "process_weather_merra2",
                 "extract_chirps_v3_rainfall", "process_weather_nasapower_chirps_v3",
                 "merge_rainfall_into_weather", "setup_cds_credentials",
                 "era5land_set_cds_key",
                 "process_soils_gsde", "process_soils_china", "process_soils_febr",
                 "process_soils_slc", "process_soils_esdb", "process_soils_openlandmap"):
        assert hasattr(dssatutils, name), f"dssatutils missing public {name}"


def test_setup_cds_credentials_writes_temp_cdsapirc():
    from dssatutils import setup_cds_credentials
    old_env = {k: os.environ.get(k) for k in ("CDSAPI_KEY", "CDSAPI_URL", "CDSAPI_RC")}
    try:
        for key in old_env:
            os.environ.pop(key, None)
        with tempfile.TemporaryDirectory() as work:
            rc = os.path.join(work, ".cdsapirc")
            meta = setup_cds_credentials(token="dummy-token", rc_path=rc, overwrite=True, prompt=False)
            text = open(rc, encoding="utf-8").read()
            assert "url: https://cds.climate.copernicus.eu/api" in text
            assert "key: dummy-token" in text
            assert meta["path"] == rc
            assert os.environ["CDSAPI_KEY"] == "dummy-token"
    finally:
        for key, value in old_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def test_rainfall_merge_helper_replaces_matching_dates_only():
    from dssatutils.weather_rainfall_merge import merge_rainfall_into_weather
    df = pd.DataFrame({"DATE": ["2010001", "2010002", "2010003"], "RAIN": [1, 1, 1]})
    n = merge_rainfall_into_weather(df, {"2010002": 9.5, "2010004": 4.0})
    assert n == 1
    assert df["RAIN"].tolist() == [1.0, 9.5, 1.0]
    assert str(df["RAIN"].dtype) == "float64"


def test_chirps_v3_path_builder_for_rnl_sat_and_prelim():
    from dssatutils.weather_chirps_v3 import _chirps_v3_file_info, _months_for_range
    fname, url = _chirps_v3_file_info(2010, 3, product="rnl", stream="final",
                                      fetch_mode="monthly_netcdf")
    assert fname == "chirps-v3.0.2010.03.days_p05.nc"
    assert "/daily/final/rnl/netcdf/byMonth/" in url
    fname, url = _chirps_v3_file_info(1998, product="sat", stream="final",
                                      fetch_mode="yearly_netcdf")
    assert fname == "chirps-v3.0.sat.1998.days_p05.nc"
    assert "/daily/final/sat/netcdf/byYear/" in url
    fname, url = _chirps_v3_file_info(2026, product="sat", stream="prelim",
                                      fetch_mode="yearly_netcdf")
    assert fname == "chirps-v3.0.sat.2026.days_p05.nc"
    assert "/daily/prelim/sat/netcdf/byYear/" in url
    assert _months_for_range(2010, 2010, months=[3]) == [(2010, 3)]


def test_chirps_v3_cog_month_filter_and_complete_cache(tmp_path, monkeypatch):
    import dssatutils.weather_chirps_v3 as chirps

    calls = []

    def fake_extract(pid, lat, lon, dates, product, stream):
        calls.append(list(dates))
        values = {
            f"{d.year}{d.timetuple().tm_yday:03d}": float(d.day)
            for d in dates
        }
        series = pd.Series(values, dtype=float)
        series.attrs["fetch_failures"] = 0
        return series

    monkeypatch.setattr(chirps, "_extract_point_cog", fake_extract)
    result = chirps._extract_chirps_v3_rain_remote_cog(
        ["P1"], np.array([0.0]), np.array([30.0]), 2010, 2010,
        "rnl", "final", str(tmp_path), months=[3],
    )
    assert len(calls) == 1 and len(calls[0]) == 31
    assert len(result["P1"]) == 31
    cache = list(tmp_path.glob("*_m03.csv"))
    assert len(cache) == 1

    def should_not_fetch(*args, **kwargs):
        raise AssertionError("a complete COG cache should be reused")

    monkeypatch.setattr(chirps, "_extract_point_cog", should_not_fetch)
    cached = chirps._extract_chirps_v3_rain_remote_cog(
        ["P1"], np.array([0.0]), np.array([30.0]), 2010, 2010,
        "rnl", "final", str(tmp_path), months=[3],
    )
    assert cached["P1"].equals(result["P1"])


def test_chirps_v3_mixed_coverage_does_not_extract_out_of_band_points(tmp_path, monkeypatch):
    import dssatutils.weather_chirps_v3 as chirps

    seen = {}

    def fake_remote(ids, lats, lons, *args, **kwargs):
        seen["ids"] = list(ids)
        return {pid: pd.Series({"2010001": 1.0}) for pid in ids}

    monkeypatch.setattr(chirps, "_extract_chirps_v3_rain_remote_cog", fake_remote)
    points = pd.DataFrame({
        "ID": ["IN", "OUT"],
        "LAT": [0.0, 65.0],
        "LONG": [30.0, 30.0],
    })
    result = chirps.extract_chirps_v3_rainfall(
        points, 2010, 2010, "ID", "LAT", "LONG", str(tmp_path),
        fetch_mode="remote_cog", months=[1],
    )
    assert seen["ids"] == ["IN"]
    assert result["IN"] == {"2010001": 1.0}
    assert result["OUT"] == {}

    try:
        chirps.extract_chirps_v3_rainfall(
            points.iloc[:1], 2010, 2010, "ID", "LAT", "LONG", str(tmp_path),
            fetch_mode="remote_cog", months=[13],
        )
    except ValueError as exc:
        assert "1..12" in str(exc)
    else:
        raise AssertionError("invalid CHIRPS month was accepted")


def test_agera5_timeseries_backend_dispatches_and_parses_csv(tmp_path, monkeypatch):
    import dssatutils.weather_agera5 as agera5

    captured = {}

    def fake_timeseries(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(agera5, "_process_weather_agera5_timeseries", fake_timeseries)
    points = pd.DataFrame({"ID": ["P1"], "LAT": [40.0], "LONG": [-90.0]})
    agera5.process_weather_agera5(
        points, 2010, 2010, str(tmp_path / "out"), "ID", "LAT", "LONG",
        1, str(tmp_path / "errors.log"), str(tmp_path / "cache"),
        agera5_backend="time-series", agera5_data_format="csv",
        agera5_timeseries_chunk_degrees=3.5,
    )
    assert captured["agera5_data_format"] == "csv"
    assert captured["agera5_timeseries_chunk_degrees"] == 3.5

    csv_path = tmp_path / "agera5.csv"
    pd.DataFrame({
        "valid_time": ["2010-01-01"],
        "latitude": [40.0],
        "longitude": [-90.0],
        "Temperature_Air_2m_Max_24h": [300.0],
        "Temperature_Air_2m_Min_24h": [280.0],
        "Solar_Radiation_Flux": [12_000_000.0],
        "Precipitation_Flux": [2.0],
        "Dew_Point_Temperature_2m_Mean_24h": [275.0],
        "Relative_Humidity_2m_15h": [60.0],
        "Wind_Speed_10m_Mean_24h": [3.0],
    }).to_csv(csv_path, index=False)
    parsed = agera5._read_agera5_timeseries_csv(str(csv_path))
    assert parsed.loc[0, "DATE"] == "2010001"
    assert abs(parsed.loc[0, "TMAX"] - 26.85) < 1e-9
    assert parsed.loc[0, "SRAD"] == 12.0


def test_agera5_timeseries_cache_chunks_are_global_grid_anchored():
    import dssatutils.weather_agera5 as agera5

    first = agera5._split_agera5_timeseries_chunks(
        np.array([33.6816]), np.array([-102.5220]), chunk_degrees=0.1
    )
    same_cell = agera5._split_agera5_timeseries_chunks(
        np.array([33.7040]), np.array([-102.4960]), chunk_degrees=0.1
    )
    adjacent = agera5._split_agera5_timeseries_chunks(
        np.array([33.7510]), np.array([-102.4490]), chunk_degrees=0.1
    )
    assert np.allclose(first[0]["area"], [33.75, -102.55, 33.65, -102.45])
    assert first[0]["area"][0] > first[0]["area"][2]
    assert first[0]["area"][3] > first[0]["area"][1]
    assert np.allclose(same_cell[0]["area"], first[0]["area"])
    assert not np.allclose(adjacent[0]["area"], first[0]["area"])


def test_agera5_rejects_unknown_backend_before_network(tmp_path):
    import dssatutils.weather_agera5 as agera5

    points = pd.DataFrame({"ID": ["P1"], "LAT": [40.0], "LONG": [-90.0]})
    try:
        agera5.process_weather_agera5(
            points, 2010, 2010, str(tmp_path / "out"), "ID", "LAT", "LONG",
            1, str(tmp_path / "errors.log"), str(tmp_path / "cache"),
            agera5_backend="not-a-backend",
        )
    except ValueError as exc:
        assert "gridded" in str(exc) and "timeseries" in str(exc)
    else:
        raise AssertionError("unknown AgERA5 backend was accepted")


def test_agera5_cache_key_includes_geographic_area(tmp_path, monkeypatch):
    import dssatutils.weather_agera5 as agera5
    import zipfile

    destinations = []
    class Client:
        def retrieve(self, dataset, request, destination):
            destinations.append(destination)
            with zipfile.ZipFile(destination, "w") as archive:
                archive.writestr("data.nc", b"valid")

    monkeypatch.setattr(agera5, "_make_cds_client", lambda _: Client())
    monkeypatch.setitem(sys.modules, "cdsapi", type("CDS", (), {})())
    first = agera5._download_agera5_var(
        "temperature", None, None, 2010, [42, -94, 40, -92], str(tmp_path))
    second = agera5._download_agera5_var(
        "temperature", None, None, 2010, [32, -102, 30, -100], str(tmp_path))
    assert first != second
    assert len(destinations) == 2


def test_gridded_filename_matching_is_exact_and_multiyear(tmp_path):
    from dssatutils.weather_gridded_common import find_nc_files
    for name in ("pr_2001.nc", "pr_2002.nc", "pressure_2001.nc"):
        (tmp_path / name).write_bytes(b"")
    found = [os.path.basename(path) for path in find_nc_files(str(tmp_path), ["pr"])]
    assert found == ["pr_2001.nc", "pr_2002.nc"]


def test_agera5_timeseries_backend_writes_weather_file(tmp_path, monkeypatch):
    import dssatutils.weather_agera5 as agera5

    dates = pd.date_range("2010-01-01", "2010-12-31", freq="D")
    csv_path = tmp_path / "agera5_timeseries.csv"
    pd.DataFrame({
        "valid_time": dates,
        "latitude": 40.0,
        "longitude": -90.0,
        "Temperature_Air_2m_Max_24h": 300.0,
        "Temperature_Air_2m_Min_24h": 280.0,
        "Solar_Radiation_Flux": 12_000_000.0,
        "Precipitation_Flux": 2.0,
        "Dew_Point_Temperature_2m_Mean_24h": 275.0,
        "Relative_Humidity_2m_15h": 60.0,
        "Wind_Speed_10m_Mean_24h": 3.0,
    }).to_csv(csv_path, index=False)
    monkeypatch.setattr(
        agera5, "_download_agera5_timeseries",
        lambda year, area, cache_dir, data_format: str(csv_path),
    )
    points = pd.DataFrame({"ID": ["P1"], "LAT": [40.0], "LONG": [-90.0]})
    out_dir = tmp_path / "out"
    agera5.process_weather_agera5(
        points, 2010, 2010, str(out_dir), "ID", "LAT", "LONG", 1,
        str(tmp_path / "errors.log"), str(tmp_path / "cache"),
        agera5_backend="timeseries",
    )
    lines = (out_dir / "P1.WTH").read_text().splitlines()
    assert len([line for line in lines if line.startswith("2010")]) == 365
    assert "  26.9" in lines[4]


def test_alderman_coordinate_aliases_and_point_geometry():
    from shapely.geometry import Point
    from dssatutils.soil_ssurgo_alderman import _alderman_coordinates

    assert _alderman_coordinates(lat=40, long=-90, required=True) == (40.0, -90.0)
    assert _alderman_coordinates(pt_geom=Point(-89, 41), required=True) == (41.0, -89.0)
    try:
        _alderman_coordinates(lat=40, lon=-90, long=-91)
    except ValueError as exc:
        assert "different values" in str(exc)
    else:
        raise AssertionError("conflicting lon/long aliases were accepted")


def test_gridded_weather_writer_and_unit_helpers():
    from dssatutils import weather_gridded_common as g
    dates = pd.date_range("2001-01-01", "2001-12-31", freq="D")
    df = pd.DataFrame({
        "DATE": [f"{d.year}{d.dayofyear:03d}" for d in dates],
        "YEAR": dates.year, "MM": dates.month,
        "SRAD": 18.0, "TMAX": 25.0, "TMIN": 12.0, "RAIN": 2.0,
        "TDEW": 9.0, "RH2M": 70.0, "WIND": 2.0,
    })
    with tempfile.TemporaryDirectory() as work:
        g.write_wth(df, "WX1", 35.0, -90.0, work, "TESTGRID", "TGRD")
        text = open(os.path.join(work, "WX1.WTH")).read()
        assert "$WEATHER DATA: TESTGRID" in text and "@  DATE" in text
    assert abs(g.convert_units(np.array([300.0]), "K", "temp")[0] - 26.85) < 1e-6
    assert abs(g.convert_units(np.array([100.0]), "W m-2", "srad")[0] - 8.64) < 1e-6


def test_convert_units_wind_and_vapour_pressure():
    from dssatutils import weather_gridded_common as g
    # 10 m -> 2 m wind (FAO-56 log profile factor 0.748).
    assert abs(g.convert_units(np.array([4.0]), "m s-1", "wind")[0] - 4.0 * 0.748) < 1e-9
    # Vapour pressure ~12.27 hPa -> dewpoint near 10 C (inverse Magnus).
    td = g.convert_units(np.array([12.27]), "hPa", "vp")[0]
    assert abs(td - 10.0) < 0.5, f"vp->dewpoint off: {td}"


def test_corrupt_netcdf_validators_reject_bogus_cache_files():
    from dssatutils.weather_gridmet import _validate_gridmet_nc
    from dssatutils.weather_nasapower_chirps import _validate_chirps_nc
    from dssatutils.weather_chirps_v3 import _validate_chirps_v3_nc
    with tempfile.TemporaryDirectory() as work:
        bogus = os.path.join(work, "bad.nc")
        with open(bogus, "wb") as fh:
            fh.write(b"not a netcdf")
        assert not _validate_gridmet_nc(bogus, "tmmn")
        assert not _validate_chirps_nc(bogus)
        assert not _validate_chirps_v3_nc(bogus)


def test_ssks_passthrough_in_sol_writer():
    """A finite source SSKS column overrides the texture-based estimate."""
    from dssatutils import soil_soilgrids_online as sg
    base = dict(ID="SK1", latitude=10.0, longitude=20.0, depth_bottom=15,
                depth_center=7.5, sand=40.0, clay=20.0, silt=40.0, bdod=1.30,
                soc_pct=1.2, cfvo=0.0, SLLL=0.10, SDUL=0.25, SSAT=0.43)
    with tempfile.TemporaryDirectory() as work:
        sg._format_dssat_sol_file(pd.DataFrame([{**base, "SSKS": 12.3}]), work,
                                  source_name="Test", source_tag="ssks")
        with_ssks = open(os.path.join(work, "SK1.SOL")).read()
        sg._format_dssat_sol_file(pd.DataFrame([base]), work,
                                  source_name="Test", source_tag="tex")
        without = open(os.path.join(work, "SK1.SOL")).read()
    assert " 12.3 " in with_ssks, "source SSKS not honored"
    assert " 12.3 " not in without, "texture fallback should not equal source SSKS"


def test_soilgrids_online_accepts_custom_id_column(tmp_path, monkeypatch):
    """The internal canonical ID must not leak into the public id_col join."""
    try:
        import geopandas as gpd
        from shapely.geometry import Point
    except ImportError:
        return

    from dssatutils import soil_soilgrids_online as sg

    rows = []
    values = {"clay": 200, "sand": 400, "silt": 400, "soc": 150,
              "bdod": 130, "cfvo": 10}
    for prop, value in values.items():
        rows.append({"prop": prop, "depth_label": "0-5cm", "depth_bottom": 5,
                     "depth_center": 2.5, "value": value})

    monkeypatch.setattr(sg, "USE_REST_API", True)
    monkeypatch.setattr(sg, "_fetch_soilgrids_rest",
                        lambda lat, lon: pd.DataFrame(rows))
    points = gpd.GeoDataFrame(
        {"point_id": ["CUSTOM_1"]}, geometry=[Point(-90.0, 40.0)], crs="EPSG:4326"
    )
    mapping = tmp_path / "soil_map.csv"
    sol_dir = tmp_path / "sol"

    sg.process_soils_soilgrids_online(points, str(mapping), str(sol_dir), "point_id")

    result = pd.read_csv(mapping)
    assert result.loc[0, "point_id"] == "CUSTOM_1"
    assert (sol_dir / "CUSTOM_1.SOL").is_file()


def test_terraclimate_monthly_disaggregated_to_daily():
    """A 12-month synthetic TerraClimate NetCDF expands to continuous daily rows."""
    xr = pytest_importorskip_xarray()
    if xr is None:
        return
    from dssatutils import process_weather_terraclimate
    lat = np.array([-10.0, -11.0, -12.0])
    lon = np.array([30.0, 31.0, 32.0])
    times = pd.date_range("2001-01-01", periods=12, freq="MS")

    def _da(value, units):
        arr = np.full((12, 3, 3), value, dtype=float)
        da = xr.DataArray(arr, coords={"time": times, "lat": lat, "lon": lon},
                          dims=["time", "lat", "lon"])
        da.attrs["units"] = units
        return da

    with tempfile.TemporaryDirectory() as work:
        nc_dir = os.path.join(work, "nc"); os.makedirs(nc_dir)
        _da(30.0, "degC").to_dataset(name="tmmx").to_netcdf(os.path.join(nc_dir, "TerraClimate_tmmx_2001.nc"))
        _da(15.0, "degC").to_dataset(name="tmmn").to_netcdf(os.path.join(nc_dir, "TerraClimate_tmmn_2001.nc"))
        _da(31.0, "mm").to_dataset(name="ppt").to_netcdf(os.path.join(nc_dir, "TerraClimate_ppt_2001.nc"))
        out_dir = os.path.join(work, "wth")
        log = os.path.join(work, "log.txt")
        points = pd.DataFrame({"ID": ["TC1"], "LAT": [-11.1], "LON": [31.1]})
        process_weather_terraclimate(points, 2001, 2001, out_dir, "ID", "LAT", "LON",
                                     1, log, nc_dir)
        lines = [l.rstrip("\n") for l in open(os.path.join(out_dir, "TC1.WTH")) if l.strip()]
        data = lines[4:]
        assert len(data) == 365, f"expected 365 continuous daily rows, got {len(data)}"
        # Continuous DOY 1..365 with no monthly gaps.
        doys = [int(r.split()[0][4:]) for r in data]
        assert doys == list(range(1, 366))
        # January precip total 31 mm spread over 31 days -> 1.0 mm/day.
        jan = data[0].split()
        assert abs(float(jan[4]) - 1.0) < 0.05, f"daily rain off: {jan[4]}"


def pytest_importorskip_xarray():
    try:
        import xarray as xr
        return xr
    except Exception:
        return None


def test_raster_soil_texture_helper():
    from dssatutils import soil_raster_common as s
    sand, silt, clay = s.texture_to_pct(9)
    assert (sand, silt, clay) == (32, 34, 34)


def test_agmip_wrapper_maps_external_sol_profiles():
    import dssatutils
    master = """*AGMIP001 AgMIP test profile A
@SITE        COUNTRY          LAT     LONG SCS FAMILY
 SITEA       Test          10.0000  20.0000 -99
@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI
    15   -99 0.100 0.250 0.430 1.00  -99 1.30 1.20 22.0 38.0
*AGMIP002 AgMIP test profile B
@SITE        COUNTRY          LAT     LONG SCS FAMILY
 SITEB       Test         -12.0000  35.0000 -99
@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI
    15   -99 0.120 0.270 0.440 1.00  -99 1.25 1.40 28.0 34.0
"""
    points = pd.DataFrame({"ID": ["P7"], "LAT": [-12.1], "LONG": [35.1]})
    with tempfile.TemporaryDirectory() as work:
        master_path = os.path.join(work, "AGMIP_TEST.SOL")
        out_dir = os.path.join(work, "sol")
        map_path = os.path.join(work, "soil.csv")
        with open(master_path, "w", encoding="latin-1") as fh:
            fh.write(master)
        mapping = dssatutils.process_soils_agmip(points, master_path, map_path, out_dir)
        assert mapping.loc[0, "SOURCE_SOIL_ID"] == "AGMIP002"
        out = os.path.join(out_dir, "00000007.SOL")
        assert os.path.exists(out)
        text = open(out, encoding="latin-1").read()
        assert text.startswith("*00000007")


def test_isdasoil_back_transform_and_writer():
    from dssatutils import soil_isdasoil as g
    raw = np.array([42, 117, 30, 255], dtype=float)  # clay, bd, oc, nodata
    assert g._back_transform("clay_content", raw)[0] == 42
    assert abs(g._back_transform("bulk_density", raw)[1] - 1.17) < 1e-9
    assert abs(g._back_transform("carbon_organic", raw)[2] - (math.exp(3.0) - 1)) < 1e-6
    assert np.isnan(g._back_transform("clay_content", raw)[3])  # 255 -> nodata
    prof = pd.DataFrame({
        "ID": ["X"] * 2, "latitude": [-0.5, -0.5], "longitude": [37.0, 37.0],
        "depth_top": [0, 20], "depth_bottom": [20, 50],
        "clay_pct": [42.0, 42.0], "sand_pct": [27.0, 27.0], "silt_pct": [31.0, 31.0],
        "om_pct": [3.29, 2.15], "bulk_density": [1.17, 1.21],
        "SLLL": [0.25, 0.25], "SDUL": [0.38, 0.38], "SSAT": [0.46, 0.46],
    })
    with tempfile.TemporaryDirectory() as work:
        g._write_sol(prof, work)
        txt = open(os.path.join(work, "X.SOL")).read()
        assert "Africa iSDAsoil" in txt and "@  SLB" in txt


def test_xavier_no_srad_conversion_and_writer():
    from dssatutils import weather_xavier as x
    # Xavier Rs is already MJ/m2/day -> the writer must NOT scale it.
    dates = pd.date_range("2015-01-01", "2015-12-31", freq="D")
    df = pd.DataFrame({
        "DATE": [f"{d.year}{d.dayofyear:03d}" for d in dates],
        "YEAR": dates.year, "MM": dates.month,
        "SRAD": 22.0, "TMAX": 32.0, "TMIN": 22.0, "RAIN": 3.0,
        "TDEW": 20.0, "RH2M": 75.0, "WIND": 2.5,
    })
    with tempfile.TemporaryDirectory() as work:
        x._write_wth(df, "BR1", -15.8, -47.9, work)
        lines = [l.rstrip("\n") for l in open(os.path.join(work, "BR1.WTH")) if l.strip()]
        assert lines[0].startswith("$WEATHER DATA: BR-DWGD/Xavier")
        assert abs(float(lines[4][7:13]) - 22.0) < 0.05   # SRAD passed through


def test_cmfd_rh_tdew_from_specific_humidity():
    from dssatutils import weather_cmfd as c
    rh, tdew = c._rh_tdew_from_shum(np.array([0.004]), np.array([10.0]), np.array([90000.0]))
    assert 30.0 < rh[0] < 70.0          # plausible RH
    assert tdew[0] < 10.0               # dew point below air temp


# --- DWD: solar radiation from sunshine duration ----------------------------

def test_dwd_srad_from_sunshine_physical():
    from dssatutils import weather_dwd as dwd
    lat = 50.0
    # Mid-winter (DOY 15) vs mid-summer (DOY 196), full sunshine.
    Ra_win, _ = dwd._extraterrestrial_radiation(lat, [15])
    Ra_sum, _ = dwd._extraterrestrial_radiation(lat, [196])
    assert Ra_sum[0] > Ra_win[0] > 0, "extraterrestrial radiation must peak in summer"

    # Clear-sky (n=N) summer radiation should land in a sane MJ/m2/day range.
    doy = np.array([196])
    _, ws = dwd._extraterrestrial_radiation(lat, doy)
    N = 24.0 / np.pi * ws
    rs = dwd._srad_from_sunshine(lat, doy, N)        # full sunshine
    assert 18.0 < rs[0] < 32.0, f"summer clear-sky SRAD out of range: {rs[0]}"

    # Zero sunshine -> only the diffuse fraction (a_s * Ra), and < clear-sky.
    rs0 = dwd._srad_from_sunshine(lat, doy, np.array([0.0]))
    assert 0.0 < rs0[0] < rs[0]
    # Missing sunshine -> NaN (written as -99 downstream).
    assert np.isnan(dwd._srad_from_sunshine(lat, doy, np.array([np.nan]))[0])


def test_dwd_tdew_from_vapour_pressure():
    from dssatutils import weather_dwd as dwd
    # Saturation VP of ~12.27 hPa corresponds to a dew point near 10 C.
    td = dwd._tdew_from_vapour_pressure(np.array([12.27]))[0]
    assert abs(td - 10.0) < 0.5, f"dew point off: {td}"


# --- E-OBS: dew point from RH + the .WTH writer -----------------------------

def test_eobs_tdew_from_rh_and_writer():
    from dssatutils import weather_eobs as eobs
    td = eobs._tdew_from_rh(np.array([23.0]), np.array([65.0]))[0]
    assert 15.0 < td < 17.0, f"dew point from RH off: {td}"

    # Writer: synthetic 1-year frame -> well-formed .WTH with the qq->MJ already applied.
    dates = pd.date_range("2019-01-01", "2019-12-31", freq="D")
    df = pd.DataFrame({
        "DATE": [f"{d.year}{d.dayofyear:03d}" for d in dates],
        "YEAR": dates.year, "MM": dates.month,
        "SRAD": 15.0, "TMAX": 20.0, "TMIN": 9.0, "RAIN": 1.0,
        "TDEW": 8.0, "RH2M": 70.0, "WIND": 3.0,
    })
    with tempfile.TemporaryDirectory() as work:
        eobs._write_wth(df, "EOBS_T", 50.0, 8.0, work)
        out = os.path.join(work, "EOBS_T.WTH")
        lines = [l.rstrip("\n") for l in open(out) if l.strip()]
        assert lines[0].startswith("$WEATHER DATA: E-OBS")
        assert lines[1].lstrip().startswith("@ INSI")
        assert lines[3].lstrip().startswith("@  DATE")
        assert len(lines[4:]) == 365
        assert lines[4][:7] == "2019001"


# --- gNATSGO: the .SOL writer (no network) ----------------------------------

def test_gnatsgo_sol_writer():
    from dssatutils import soil_gnatsgo as g
    profile = pd.DataFrame({
        "ID": ["00000099"] * 2, "latitude": [42.0, 42.0], "longitude": [-93.0, -93.0],
        "depth_top": [0, 5], "depth_bottom": [5, 20],
        "clay_pct": [25.0, 26.0], "sand_pct": [35.0, 34.0], "silt_pct": [40.0, 40.0],
        "om_pct": [3.0, 2.5], "bulk_density": [1.3, 1.35],
        "SLLL": [0.12, 0.13], "SDUL": [0.28, 0.29], "SSAT": [0.45, 0.45],
    })
    with tempfile.TemporaryDirectory() as work:
        g._write_sol(profile, work)
        out = os.path.join(work, "00000099.SOL")
        text = open(out).read()
        assert "USA gNATSGO Soil Profiles" in text
        assert "@  SLB  SLMH  SLLL  SDUL  SSAT" in text
        # Two layers + headers; depth column shows 5 and 20.
        assert "     5   -99" in text and "    20   -99" in text


# --- Cross-language parity: the R twins carry the same API + algorithms ------

def _read(rel):
    return open(os.path.join(_WORKSPACE, "dssatutils", rel), encoding="utf-8", errors="replace").read()


def test_r_python_parity_markers():
    checks = {
        "R/soil_gnatsgo.R": ("process_soils_gnatsgo", "mukey.wcs", "no-tabular", "Saxton"),
        "python/dssatutils/soil_gnatsgo.py": ("process_soils_gnatsgo", "gnatsgo", "no-tabular", "_saxton_rawls"),
        "R/weather_dwd.R": ("process_weather_dwd", "Angstrom", "Magnus" if False else "6.1094"),
        "python/dssatutils/weather_dwd.py": ("process_weather_dwd", "Angstrom", "6.1094"),
        "R/weather_eobs.R": ("process_weather_eobs", "0.0864", "tx"),
        "python/dssatutils/weather_eobs.py": ("process_weather_eobs", "0.0864", "tx"),
        "R/soil_isdasoil.R": ("process_soils_isdasoil", "isdasoil", "carbon_organic"),
        "python/dssatutils/soil_isdasoil.py": ("process_soils_isdasoil", "isdasoil", "expm1"),
        "R/soil_lucas.R": ("process_soils_lucas", "LUCAS", "esdac", "2.65"),
        "python/dssatutils/soil_lucas.py": ("process_soils_lucas", "LUCAS", "esdac", "2.65"),
        "R/weather_xavier.R": ("process_weather_xavier", "XAVR", "BR-DWGD"),
        "python/dssatutils/weather_xavier.py": ("process_weather_xavier", "XAVR", "BR-DWGD"),
        "R/weather_cmfd.R": ("process_weather_cmfd", "CMFD", "0.0864"),
        "python/dssatutils/weather_cmfd.py": ("process_weather_cmfd", "CMFD", "0.0864"),
        "R/soil_agmip.R": ("process_soils_agmip", "AgMIP", "10.7910/DVN/1PEEY0"),
        "python/dssatutils/soil_agmip.py": ("process_soils_agmip", "AgMIP", "10.7910/DVN/1PEEY0"),
        "R/weather_chelsa_w5e5.R": ("process_weather_chelsa_w5e5", "CHELSA-W5E5"),
        "python/dssatutils/weather_chelsa_w5e5.py": ("process_weather_chelsa_w5e5", "CHELSA-W5E5"),
        "R/weather_agmip.R": ("process_weather_agmerra", "process_weather_agcfsr", "AgMERRA", "AgCFSR"),
        "python/dssatutils/weather_agmip.py": ("process_weather_agmerra", "process_weather_agcfsr", "AgMERRA", "AgCFSR"),
        "R/weather_silo.R": ("process_weather_silo", "SILO"),
        "python/dssatutils/weather_silo.py": ("process_weather_silo", "SILO"),
        "R/weather_prism.R": ("process_weather_prism", "PRISM"),
        "python/dssatutils/weather_prism.py": ("process_weather_prism", "PRISM"),
        "R/soil_hihydrosoil.R": ("process_soils_hihydrosoil", "HiHydroSoil"),
        "python/dssatutils/soil_hihydrosoil.py": ("process_soils_hihydrosoil", "HiHydroSoil"),
        "R/soil_slga.R": ("process_soils_slga", "SLGA"),
        "python/dssatutils/soil_slga.py": ("process_soils_slga", "SLGA"),
        "R/weather_mswx.R": ("process_weather_mswx", "MSWX"),
        "python/dssatutils/weather_mswx.py": ("process_weather_mswx", "MSWX"),
        "R/weather_mswep.R": ("process_weather_mswep", "MSWEP"),
        "python/dssatutils/weather_mswep.py": ("process_weather_mswep", "MSWEP"),
        "R/weather_crujra.R": ("process_weather_crujra", "CRU-JRA"),
        "python/dssatutils/weather_crujra.py": ("process_weather_crujra", "CRU-JRA"),
        "R/weather_terraclimate.R": ("process_weather_terraclimate", "TerraClimate"),
        "python/dssatutils/weather_terraclimate.py": ("process_weather_terraclimate", "TerraClimate"),
        "R/soil_wise30sec.R": ("process_soils_wise30sec", "WISE30sec"),
        "python/dssatutils/soil_wise30sec.py": ("process_soils_wise30sec", "WISE30sec"),
        "R/soil_wosis.R": ("process_soils_wosis", "WoSIS"),
        "python/dssatutils/soil_wosis.py": ("process_soils_wosis", "WoSIS"),
        "R/weather_aphrodite.R": ("process_weather_aphrodite", "APHRODITE"),
        "python/dssatutils/weather_aphrodite.py": ("process_weather_aphrodite", "APHRODITE"),
        "R/weather_anusplin.R": ("process_weather_anusplin", "ANUSPLIN"),
        "python/dssatutils/weather_anusplin.py": ("process_weather_anusplin", "ANUSPLIN"),
        "R/weather_tamsat.R": ("process_weather_tamsat", "TAMSAT"),
        "python/dssatutils/weather_tamsat.py": ("process_weather_tamsat", "TAMSAT"),
        "R/weather_ghcn.R": ("process_weather_ghcn", "GHCN"),
        "python/dssatutils/weather_ghcn.py": ("process_weather_ghcn", "GHCN"),
        "R/weather_pgf.R": ("process_weather_pgf", "PGF"),
        "python/dssatutils/weather_pgf.py": ("process_weather_pgf", "PGF"),
        "R/weather_merra2.R": ("process_weather_merra2", "MERRA-2"),
        "python/dssatutils/weather_merra2.py": ("process_weather_merra2", "MERRA-2"),
        "R/weather_rainfall_merge.R": ("merge_rainfall_into_weather", "n_replaced"),
        "python/dssatutils/weather_rainfall_merge.py": ("merge_rainfall_into_weather", "date_col", "rain_col"),
        "R/weather_chirps_v3.R": ("extract_chirps_v3_rainfall", "process_weather_nasapower_chirps_v3", "rnl", "sat"),
        "python/dssatutils/weather_chirps_v3.py": ("extract_chirps_v3_rainfall", "process_weather_nasapower_chirps_v3", "rnl", "sat"),
        "R/soil_gsde.R": ("process_soils_gsde", "GSDE"),
        "python/dssatutils/soil_gsde.py": ("process_soils_gsde", "GSDE"),
        "R/soil_china.R": ("process_soils_china", "BNU"),
        "python/dssatutils/soil_china.py": ("process_soils_china", "BNU"),
        "R/soil_febr.R": ("process_soils_febr", "FEBR"),
        "python/dssatutils/soil_febr.py": ("process_soils_febr", "FEBR"),
        "R/soil_slc.R": ("process_soils_slc", "Soil Landscapes of Canada"),
        "python/dssatutils/soil_slc.py": ("process_soils_slc", "Soil Landscapes of Canada"),
        "R/soil_esdb.R": ("process_soils_esdb", "ESDB"),
        "python/dssatutils/soil_esdb.py": ("process_soils_esdb", "ESDB"),
        "R/soil_openlandmap.R": ("process_soils_openlandmap", "OpenLandMap"),
        "python/dssatutils/soil_openlandmap.py": ("process_soils_openlandmap", "OpenLandMap"),
    }
    for rel, markers in checks.items():
        src = _read(rel)
        for m in markers:
            if not m:
                continue
            assert m in src, f"{rel} missing expected marker {m!r}"


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"  [ok] {name}")
            except Exception as exc:  # noqa: BLE001
                print(f"  [FAIL] {name}: {exc}")
                fails += 1
    print("ALL NEW-SOURCE TESTS PASSED." if not fails else f"{fails} TEST(S) FAILED.")
    sys.exit(1 if fails else 0)
