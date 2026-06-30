#!/usr/bin/env python3
# File: tests/test_comprehensive.py
# ---------------------------------------------------------------------------
# Comprehensive test suite (Python) — tests all 6 weather + 4 soil sources.
# Uses mocks and synthetic data to run completely offline (no API keys, no
# network, no large rasters).
#
# Run: python -m pytest tests/test_comprehensive.py -v
# ---------------------------------------------------------------------------

import os
import sys
import tempfile
import shutil
import sqlite3
import unittest
from unittest.mock import patch, MagicMock

import numpy as np
import pandas as pd

# Insert the 'python' folder where 'dssatutils' package resides
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_REPO, "python"))


class TestComprehensive(unittest.TestCase):
    """Offline mocked tests for every weather and soil source in dssatutils."""

    def setUp(self):
        self.work = tempfile.mkdtemp(prefix="dssat_comprehensive_test_")
        self.shapefile = pd.DataFrame({
            "ID": ["TEST_1"],
            "LAT": [40.0],
            "LONG": [-90.0],
        })

    def tearDown(self):
        shutil.rmtree(self.work, ignore_errors=True)

    # ------------------------------------------------------------------
    # Helper: verify a .WTH file exists and has basic structure
    # ------------------------------------------------------------------
    def _assert_wth_valid(self, path):
        self.assertTrue(os.path.exists(path), f"WTH file not found: {path}")
        with open(path) as f:
            content = f.read()
        self.assertIn("$WEATHER DATA", content)
        self.assertIn("@  DATE", content)
        lines = content.strip().split("\n")
        self.assertGreater(len(lines), 4, "WTH file has too few lines")

    # ------------------------------------------------------------------
    # Helper: verify a .SOL file exists and has basic structure
    # ------------------------------------------------------------------
    def _assert_sol_valid(self, path):
        self.assertTrue(os.path.exists(path), f"SOL file not found: {path}")
        with open(path) as f:
            content = f.read()
        # SOL files may start with *SOILS: header or directly with the profile
        has_header = "*SOILS:" in content or content.startswith("*")
        self.assertTrue(has_header, f"SOL file missing header: {path}")
        self.assertIn("@  SLB", content)

    # ===================================================================
    # WEATHER SOURCES
    # ===================================================================

    @patch("dssatutils.weather_openmeteo._fetch_open_meteo")
    @patch("dssatutils.weather_openmeteo.ProcessPoolExecutor")
    def test_weather_openmeteo(self, mock_pool_cls, mock_fetch):
        """Open-Meteo: mock the fetch function and bypass ProcessPoolExecutor."""
        from dssatutils import process_weather_openmeteo

        # Build a synthetic DataFrame mimicking _fetch_open_meteo output
        dates = pd.date_range("2010-01-01", "2010-12-31", freq="D")
        n = len(dates)
        mock_df = pd.DataFrame({
            "time": dates,
            "temperature_2m_max": np.random.uniform(20, 30, n),
            "temperature_2m_min": np.random.uniform(5, 15, n),
            "precipitation_sum": np.random.uniform(0, 10, n),
            "shortwave_radiation_sum": np.random.uniform(10, 25, n),
            "wind_speed_10m_max": np.random.uniform(1, 8, n),
        })
        mock_df["YEAR"] = dates.year
        mock_df["MM"] = dates.month
        mock_df["DOY"] = dates.day_of_year
        mock_df["DATE"] = mock_df["YEAR"].astype(str) + mock_df["DOY"].astype(str).str.zfill(3)
        mock_fetch.return_value = mock_df

        # Make ProcessPoolExecutor run synchronously in the main process
        mock_pool = MagicMock()
        mock_pool_cls.return_value.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)

        # Instead of mocking the pool, just directly call _process_single_point
        from dssatutils.weather_openmeteo import _process_single_point
        out_dir = os.path.join(self.work, "openmeteo")
        os.makedirs(out_dir, exist_ok=True)
        _process_single_point({
            "latitude": 40.0, "longitude": -90.0, "point_id": "TEST_1",
            "output_dir": out_dir, "start_date": "2010-01-01",
            "end_date": "2010-12-31",
            "log_file": os.path.join(self.work, "error.log"),
        })

        self._assert_wth_valid(os.path.join(out_dir, "TEST_1.WTH"))

    @patch("dssatutils.weather_nasapower._fetch_nasa_power")
    @patch("dssatutils.weather_nasapower.ProcessPoolExecutor")
    def test_weather_nasapower(self, mock_pool_cls, mock_fetch):
        """NASA POWER: mock the fetch function and bypass ProcessPoolExecutor."""
        from dssatutils.weather_nasapower import _process_single_point

        dates = pd.date_range("2010-01-01", "2010-12-31", freq="D")
        n = len(dates)
        mock_df = pd.DataFrame({
            "T2M_MAX": np.random.uniform(20, 35, n),
            "T2M_MIN": np.random.uniform(5, 15, n),
            "ALLSKY_SFC_SW_DWN": np.random.uniform(10, 25, n),
            "PRECTOTCORR": np.random.uniform(0, 10, n),
            "T2MDEW": np.random.uniform(5, 15, n),
            "RH2M": np.random.uniform(50, 90, n),
            "WS2M": np.random.uniform(1, 5, n),
            "YEAR": dates.year,
            "MM": dates.month,
            "DOY": dates.day_of_year,
        })
        mock_df["DATE"] = mock_df["YEAR"].astype(str) + mock_df["DOY"].astype(str).str.zfill(3)
        mock_fetch.return_value = mock_df

        out_dir = os.path.join(self.work, "nasapower")
        os.makedirs(out_dir, exist_ok=True)
        _process_single_point({
            "latitude": 40.0, "longitude": -90.0, "point_id": "TEST_1",
            "output_dir": out_dir, "start_date": "20100101",
            "end_date": "20101231",
            "log_file": os.path.join(self.work, "error.log"),
        })

        self._assert_wth_valid(os.path.join(out_dir, "TEST_1.WTH"))

    @patch("dssatutils.weather_nasapower_chirps._download_chirps_year")
    @patch("dssatutils.weather_nasapower_chirps._extract_chirps_rain")
    @patch("dssatutils.weather_nasapower_chirps._fetch_nasa_power")
    def test_weather_nasapower_chirps(self, mock_fetch_nasa, mock_extract_chirps, mock_download_chirps):
        """NASA POWER + CHIRPS hybrid: mock both data sources."""
        from dssatutils import process_weather_nasapower_chirps

        mock_download_chirps.return_value = "dummy_chirps.nc"

        # Create 365 days of CHIRPS rain data
        dates = pd.date_range("2010-01-01", "2010-12-31", freq="D")
        n = len(dates)
        rain_keys = [f"{d.year}{d.day_of_year:03d}" for d in dates]
        rain_vals = np.random.uniform(0, 15, n)
        mock_extract_chirps.return_value = {
            "TEST_1": pd.Series(dict(zip(rain_keys, rain_vals)))
        }

        # Create NASA POWER base data
        mock_df = pd.DataFrame({
            "YEAR": dates.year,
            "MM": dates.month,
            "DOY": dates.day_of_year,
            "DATE": [f"{d.year}{d.day_of_year:03d}" for d in dates],
            "TMAX": np.random.uniform(20, 35, n),
            "TMIN": np.random.uniform(5, 15, n),
            "SRAD": np.random.uniform(10, 25, n),
            "RAIN": np.random.uniform(0, 10, n),  # will be replaced by CHIRPS
            "TDEW": np.random.uniform(5, 15, n),
            "RH2M": np.random.uniform(50, 90, n),
            "WIND": np.random.uniform(1, 5, n),
        })
        mock_fetch_nasa.return_value = mock_df

        out_dir = os.path.join(self.work, "chirps")
        process_weather_nasapower_chirps(
            shapefile=self.shapefile,
            start_year=2010,
            end_year=2010,
            output_dir=out_dir,
            id_col="ID",
            lat_col="LAT",
            lon_col="LONG",
            n_cores=1,
            log_file=os.path.join(self.work, "error.log"),
            chirps_cache_dir=os.path.join(self.work, "chirps_cache"),
        )
        self._assert_wth_valid(os.path.join(out_dir, "TEST_1.WTH"))

    @patch("dssatutils.weather_chirps_v3._download_chirps_v3_file")
    @patch("dssatutils.weather_chirps_v3._extract_chirps_v3_rain")
    @patch("dssatutils.weather_chirps_v3._fetch_nasa_power")
    def test_weather_nasapower_chirps_v3(self, mock_fetch_nasa, mock_extract_chirps, mock_download_chirps):
        """NASA POWER + CHIRPS v3 hybrid: mock rnl/sat-capable rainfall layer."""
        from dssatutils import process_weather_nasapower_chirps_v3

        mock_download_chirps.return_value = "dummy_chirps_v3.nc"

        dates = pd.date_range("2010-01-01", "2010-12-31", freq="D")
        n = len(dates)
        rain_keys = [f"{d.year}{d.day_of_year:03d}" for d in dates]
        rain_vals = np.linspace(1.0, 12.0, n)
        mock_extract_chirps.return_value = {
            "TEST_1": pd.Series(dict(zip(rain_keys, rain_vals)))
        }

        mock_df = pd.DataFrame({
            "YEAR": dates.year,
            "MM": dates.month,
            "DOY": dates.day_of_year,
            "DATE": [f"{d.year}{d.day_of_year:03d}" for d in dates],
            "TMAX": np.random.uniform(20, 35, n),
            "TMIN": np.random.uniform(5, 15, n),
            "SRAD": np.random.uniform(10, 25, n),
            "RAIN": np.zeros(n),
            "TDEW": np.random.uniform(5, 15, n),
            "RH2M": np.random.uniform(50, 90, n),
            "WIND": np.random.uniform(1, 5, n),
        })
        mock_fetch_nasa.return_value = mock_df

        out_dir = os.path.join(self.work, "chirps_v3")
        process_weather_nasapower_chirps_v3(
            shapefile=self.shapefile,
            start_year=2010,
            end_year=2010,
            output_dir=out_dir,
            id_col="ID",
            lat_col="LAT",
            lon_col="LONG",
            n_cores=1,
            log_file=os.path.join(self.work, "error.log"),
            chirps_cache_dir=os.path.join(self.work, "chirps_v3_cache"),
            chirps_product="sat",
            chirps_stream="final",
        )
        path = os.path.join(out_dir, "TEST_1.WTH")
        self._assert_wth_valid(path)
        with open(path) as fh:
            txt = fh.read()
        self.assertIn("CHIRPS-v3 final/sat", txt)
        self.assertIn("NCV3", txt)

    @patch("dssatutils.weather_daymet._download_daymet")
    @patch("dssatutils.weather_daymet.ProcessPoolExecutor")
    def test_weather_daymet(self, mock_pool_cls, mock_download):
        """Daymet: mock the CSV download and bypass ProcessPoolExecutor."""
        from dssatutils.weather_daymet import _process_single_point

        dates = pd.date_range("2010-01-01", "2010-12-31", freq="D")
        n = len(dates)
        mock_download.return_value = pd.DataFrame({
            "year": dates.year,
            "yday": dates.day_of_year,
            "tmax": np.random.uniform(20, 30, n),
            "tmin": np.random.uniform(5, 15, n),
            "prcp": np.random.uniform(0, 10, n),
            "srad": np.random.uniform(200, 400, n),  # W/m²
            "dayl": np.random.uniform(35000, 45000, n),  # seconds
            "vp": np.random.uniform(800, 1500, n),  # Pa
        })

        out_dir = os.path.join(self.work, "daymet")
        os.makedirs(out_dir, exist_ok=True)
        _process_single_point({
            "latitude": 40.0, "longitude": -90.0, "point_id": "TEST_1",
            "output_dir": out_dir, "start_year": 2010, "end_year": 2010,
            "log_file": os.path.join(self.work, "error.log"),
        })

        self._assert_wth_valid(os.path.join(out_dir, "TEST_1.WTH"))

    @patch("dssatutils.weather_gridmet._download_nc")
    def test_weather_gridmet(self, mock_download_nc):
        """GridMET: mock NetCDF downloads and create synthetic xarray datasets."""
        import xarray as xr
        from dssatutils import process_weather_gridmet

        cache_dir = os.path.join(self.work, "gridmet_cache")
        os.makedirs(cache_dir, exist_ok=True)

        # Create synthetic NetCDF files for each variable and year
        dates = pd.date_range("2010-01-01", "2010-12-31", freq="D")
        n = len(dates)
        lat = np.array([39.5, 40.0, 40.5])  # 3 lat cells
        lon = np.array([-90.5, -90.0, -89.5])  # 3 lon cells

        var_configs = {
            "tmmn": ("TMIN", np.random.uniform(273, 288, (n, 3, 3))),  # Kelvin
            "tmmx": ("TMAX", np.random.uniform(293, 308, (n, 3, 3))),  # Kelvin
            "pr":   ("RAIN", np.random.uniform(0, 15, (n, 3, 3))),     # mm/day
            "srad": ("SRAD", np.random.uniform(100, 350, (n, 3, 3))),  # W/m²
        }

        for abbrev, (_, data) in var_configs.items():
            fname = f"{abbrev}_2010.nc"
            fpath = os.path.join(cache_dir, fname)
            ds = xr.Dataset(
                {abbrev: (["day", "lat", "lon"], data)},
                coords={"day": dates, "lat": lat, "lon": lon},
            )
            ds.to_netcdf(fpath)
            ds.close()

        # Mock _download_nc to just return True (files already created)
        mock_download_nc.return_value = True

        out_dir = os.path.join(self.work, "gridmet")
        process_weather_gridmet(
            shapefile=self.shapefile,
            start_year=2010,
            end_year=2010,
            output_dir=out_dir,
            id_col="ID",
            lat_col="LAT",
            lon_col="LONG",
            n_cores=1,
            log_file=os.path.join(self.work, "error.log"),
            gridmet_cache_dir=cache_dir,
        )

        self._assert_wth_valid(os.path.join(out_dir, "TEST_1.WTH"))

    @patch("dssatutils.weather_agera5._download_agera5_var")
    @patch("dssatutils.weather_agera5._open_agera5")
    def test_weather_agera5(self, mock_open, mock_download):
        """AgERA5: mock CDS download and NetCDF parsing."""
        import xarray as xr
        from dssatutils import process_weather_agera5

        # Build synthetic xarray datasets for each variable
        dates = pd.date_range("2010-01-01", "2010-12-31", freq="D")
        n = len(dates)
        lat = np.array([40.0])
        lon = np.array([-90.0])

        var_data = {
            "TMAX": np.random.uniform(293, 308, (n, 1, 1)),  # K
            "TMIN": np.random.uniform(273, 288, (n, 1, 1)),  # K
            "SRAD": np.random.uniform(5e6, 25e6, (n, 1, 1)),  # J/m²/day
            "RAIN": np.random.uniform(0, 10, (n, 1, 1)),  # mm
            "TDEW": np.random.uniform(273, 288, (n, 1, 1)),  # K
            "RH2M": np.random.uniform(50, 90, (n, 1, 1)),  # %
            "WIND": np.random.uniform(1, 8, (n, 1, 1)),  # m/s
        }

        # Mock download to return a fake path
        mock_download.return_value = "/fake/agera5.nc"

        # Mock _open_agera5 to return appropriate xarray datasets
        call_count = [0]
        var_order = list(var_data.keys())

        def mock_open_side_effect(path):
            idx = call_count[0] % len(var_order)
            vname = var_order[idx]
            call_count[0] += 1
            ds = xr.Dataset(
                {"data_var": (["time", "lat", "lon"], var_data[vname])},
                coords={"time": dates, "lat": lat, "lon": lon},
            )
            return ds

        mock_open.side_effect = mock_open_side_effect

        cache_dir = os.path.join(self.work, "agera5_cache")
        out_dir = os.path.join(self.work, "agera5")

        process_weather_agera5(
            shapefile=self.shapefile,
            start_year=2010,
            end_year=2010,
            output_dir=out_dir,
            id_col="ID",
            lat_col="LAT",
            lon_col="LONG",
            n_cores=1,
            log_file=os.path.join(self.work, "error.log"),
            agera5_cache_dir=cache_dir,
        )

        self._assert_wth_valid(os.path.join(out_dir, "TEST_1.WTH"))

    @patch("dssatutils.weather_era5land._download_era5_land_point_csv")
    def test_weather_era5land(self, mock_download):
        """ERA5-Land: mock download to write a synthetic CSV and run parsing/aggregation."""
        from dssatutils import process_weather_era5_land

        def mock_download_side_effect(latitude, longitude, start_date_str, end_date_str, target_file, cds_user):
            # Create a synthetic hourly dataset for 2 days
            dates = pd.date_range("2010-01-01 00:00:00", "2010-01-02 23:00:00", freq="h")
            n = len(dates)
            df = pd.DataFrame({
                "time": dates,
                "2m_temperature": np.random.uniform(280, 295, n),  # K
                "2m_dewpoint_temperature": np.random.uniform(275, 285, n),  # K
                "total_precipitation": np.random.uniform(0, 0.005, n),  # m
                "surface_solar_radiation_downwards": np.random.uniform(1e5, 1e6, n),  # J/m2
                "10m_u_component_of_wind": np.random.uniform(-5, 5, n),
                "10m_v_component_of_wind": np.random.uniform(-5, 5, n)
            })
            os.makedirs(os.path.dirname(target_file), exist_ok=True)
            df.to_csv(target_file, index=False)

        mock_download.side_effect = mock_download_side_effect

        out_dir = os.path.join(self.work, "era5land")
        process_weather_era5_land(
            shapefile=self.shapefile,
            start_year=2010,
            end_year=2010,
            output_dir=out_dir,
            id_col="ID",
            lat_col="LAT",
            lon_col="LONG",
            n_cores=1,
            log_file=os.path.join(self.work, "error.log")
        )

        self._assert_wth_valid(os.path.join(out_dir, "TEST_1.WTH"))

    # ===================================================================
    # SOIL SOURCES
    # ===================================================================

    def test_soil_soilgrids(self):
        """SoilGrids 10K (offline): match nearest profile from a master .SOL."""
        from dssatutils import process_soils_soilgrids

        dummy_sol = os.path.join(self.work, "dummy_master.SOL")
        with open(dummy_sol, "w") as f:
            f.write(
                "*SOILS: Dummy Master File\n"
                "*TEST0001     TEST_SOIL       40.000   -90.000\n"
                "@SITE        COUNTRY          LAT     LONG SCS FAMILY\n"
                " TEST_SOIL   World         40.000   -90.000 \n"
                "@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE\n"
                "    BN   .13     6    .6    73     1     1 IB001 IB001 IB001\n"
                "@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC\n"
                "     5   -99 0.100 0.200 0.300  1.00  10.0  1.40  1.00  20.0  40.0   0.0   -99   -99   -99   -99   -99\n"
            )

        out_csv = os.path.join(self.work, "soil_map.csv")
        out_sol_dir = os.path.join(self.work, "individual_sol")

        process_soils_soilgrids(
            grid_points=self.shapefile,
            source_sol_file=dummy_sol,
            output_csv_path=out_csv,
            output_sol_dir=out_sol_dir,
            id_col="ID",
            numeric_only_ids=False,
        )

        self.assertTrue(os.path.exists(out_csv))
        self._assert_sol_valid(os.path.join(out_sol_dir, "TEST_1.SOL"))

    @patch("dssatutils.soil_soilgrids_online.requests.get")
    def test_soil_soilgrids_online(self, mock_get):
        """SoilGrids 2.0 REST API: mock HTTP responses with correct structure."""
        from dssatutils import process_soils_soilgrids_online
        import dssatutils.soil_soilgrids_online as sg_module

        # Force REST API mode
        sg_module.USE_REST_API = True

        # Build the mock REST API response matching the actual ISRIC format:
        # Each layer has "depths" as a list of dicts with "label" and "values"
        mock_response = MagicMock()
        mock_response.status_code = 200

        depths_template = [
            {"label": "0-5cm",     "values": {"mean": None}},
            {"label": "5-15cm",    "values": {"mean": None}},
            {"label": "15-30cm",   "values": {"mean": None}},
            {"label": "30-60cm",   "values": {"mean": None}},
            {"label": "60-100cm",  "values": {"mean": None}},
            {"label": "100-200cm", "values": {"mean": None}},
        ]

        # Property values in SoilGrids native units:
        # clay, sand, silt: g/kg (divide by 10 for %)
        # soc: dg/kg (divide by 100 for %)
        # bdod: cg/cm³ (divide by 100 for g/cm³)
        # cfvo: cm³/dm³ (divide by 10 for %)
        prop_values = {
            "clay": [200, 210, 220, 230, 240, 250],  # g/kg -> 20-25%
            "sand": [400, 390, 380, 370, 360, 350],  # g/kg -> 35-40%
            "silt": [400, 400, 400, 400, 400, 400],  # g/kg -> 40%
            "soc":  [150, 120, 100, 80, 50, 30],     # dg/kg -> 0.3-1.5%
            "bdod": [130, 135, 140, 145, 150, 155],   # cg/cm³ -> 1.3-1.55 g/cm³
            "cfvo": [10, 10, 20, 20, 30, 30],        # cm³/dm³ -> 1-3%
        }

        layers = []
        for prop, values in prop_values.items():
            import copy
            depths = copy.deepcopy(depths_template)
            for i, val in enumerate(values):
                depths[i]["values"]["mean"] = val
            layers.append({"name": prop, "depths": depths})

        mock_response.json.return_value = {
            "properties": {"layers": layers}
        }
        mock_get.return_value = mock_response

        # Create a GeoDataFrame input
        try:
            import geopandas as gpd
            from shapely.geometry import Point
            gdf = gpd.GeoDataFrame(
                self.shapefile,
                geometry=[Point(-90.0, 40.0)],
                crs="EPSG:4326",
            )
        except ImportError:
            self.skipTest("geopandas/shapely required for SoilGrids online test")

        out_csv = os.path.join(self.work, "soil_map.csv")
        out_sol_dir = os.path.join(self.work, "individual_sol")

        process_soils_soilgrids_online(
            gridfile=gdf,
            soilfile_csv_path=out_csv,
            output_sol_dir=out_sol_dir,
            id_col="ID",
        )

        self.assertTrue(os.path.exists(out_csv))
        self._assert_sol_valid(os.path.join(out_sol_dir, "TEST_1.SOL"))

    @patch("dssatutils.soil_ssurgo.requests.post")
    @patch("dssatutils.soil_ssurgo._sda_spatial_mukeys")
    def test_soil_ssurgo(self, mock_mukeys, mock_post):
        """USDA SSURGO SDA: mock spatial query and attribute queries."""
        from dssatutils import process_soils_ssurgo

        mock_mukeys.return_value = ["12345"]

        def mock_sda_post_handler(url, data=None, json=None, timeout=None, **kwargs):
            sql = ""
            if data:
                sql = data.get("query", "") if isinstance(data, dict) else str(data)
            if json:
                sql = json.get("query", "")
            mock_resp = MagicMock()
            mock_resp.status_code = 200
            if "brockdepmin" in sql.lower():
                mock_resp.json.return_value = {
                    "Table": [
                        ["mukey", "brockdepmin"],
                        ["12345", "200.0"],
                    ]
                }
            else:
                mock_resp.json.return_value = {
                    "Table": [
                        ["mukey", "cokey", "comppct_r", "hzdept_r", "hzdepb_r",
                         "claytotal_r", "sandtotal_r", "om_r", "dbthirdbar_r"],
                        ["12345", "cokey1", "100", "0", "15", "20.0", "40.0", "1.5", "1.4"],
                        ["12345", "cokey1", "100", "15", "30", "22.0", "38.0", "1.2", "1.45"],
                    ]
                }
            return mock_resp

        mock_post.side_effect = mock_sda_post_handler

        out_csv = os.path.join(self.work, "soil_map.csv")
        out_sol_dir = os.path.join(self.work, "individual_sol")

        process_soils_ssurgo(
            grid_points=self.shapefile,
            output_dir_csv=out_csv,
            output_dir_individual=out_sol_dir,
            n_cores=1,
            id_col="ID",
            lat_col="LAT",
            long_col="LONG",
        )

        self.assertTrue(os.path.exists(out_csv))

    @patch("dssatutils.soil_ssurgo_alderman.requests.post")
    @patch("dssatutils.soil_ssurgo_alderman._sda_spatial_mukeys")
    def test_soil_ssurgo_alderman(self, mock_mukeys, mock_post):
        """USDA SSURGO SDA (Alderman): mock spatial query and attribute queries."""
        from dssatutils import process_soils_ssurgo_alderman

        mock_mukeys.return_value = ["12345"]

        def mock_sda_post_handler(url, data=None, json=None, timeout=None, **kwargs):
            sql = ""
            if data:
                sql = data.get("query", "") if isinstance(data, dict) else str(data)
            if json:
                sql = json.get("query", "")
            mock_resp = MagicMock()
            mock_resp.status_code = 200
            if "brockdepmin" in sql.lower():
                mock_resp.json.return_value = {
                    "Table": [
                        ["mukey", "brockdepmin"],
                        ["12345", "200.0"],
                    ]
                }
            elif "component" in sql.lower():
                mock_resp.json.return_value = {
                    "Table": [
                        ["compname", "cokey", "mukey", "comppct_r", "hydgrp", "slope_r", "drainage", "albedodry_r"],
                        ["Miami", "cokey1", "12345", "100", "B", "2.0", "Well drained", "0.13"]
                    ]
                }
            else:
                mock_resp.json.return_value = {
                    "Table": [
                        ["hzdept_r", "hzdepb_r", "dbovendry_r", "dbtenthbar_r", "dbthirdbar_r", "dbfifteenbar_r",
                         "wsatiated_r", "wtenthbar_r", "wthirdbar_r", "partdensity", "ksat_r", "wfifteenbar_r",
                         "sandtotal_r", "claytotal_r", "silttotal_r", "om_r", "hzname", "fragvol_r", "cokey"],
                        ["0", "15", "1.45", "1.4", "1.4", "1.5", "45", "25", "20", "2.65", "15", "10", "40", "20", "40", "1.5", "Ap", None, "cokey1"]
                    ]
                }
            return mock_resp

        mock_post.side_effect = mock_sda_post_handler

        out_csv = os.path.join(self.work, "soil_map_alderman.csv")
        out_sol_dir = os.path.join(self.work, "individual_sol_alderman")

        process_soils_ssurgo_alderman(
            grid_points=self.shapefile,
            output_dir_csv=out_csv,
            output_dir_individual=out_sol_dir,
            n_cores=1,
            id_col="ID",
            lat_col="LAT",
            long_col="LONG",
        )

        self.assertTrue(os.path.exists(out_csv))

    @patch("dssatutils.soil_hwsd._sample_smu_ids")
    def test_soil_hwsd(self, mock_sample):
        """FAO HWSD v2: mock raster sampling, create real SQLite DB."""
        from dssatutils import process_soils_hwsd

        mock_sample.return_value = np.array([1])

        # Create a dummy raster file (just needs to exist for the file check)
        dummy_raster = os.path.join(self.work, "dummy_hwsd.tif")
        with open(dummy_raster, "w") as f:
            f.write("")

        # Create a real SQLite DB with HWSD2_LAYERS schema
        dummy_db = os.path.join(self.work, "dummy_hwsd.sqlite")
        conn = sqlite3.connect(dummy_db)
        conn.execute(
            "CREATE TABLE HWSD2_LAYERS ("
            "HWSD2_SMU_ID INTEGER, SEQUENCE INTEGER, SHARE REAL, "
            "TOPDEP REAL, BOTDEP REAL, SAND REAL, CLAY REAL, SILT REAL, "
            "BULK REAL, ORG_CARBON REAL, COARSE REAL)"
        )
        # Insert 3 layers for SMU 1
        conn.execute(
            "INSERT INTO HWSD2_LAYERS VALUES (1, 1, 100, 0, 15, 40.0, 20.0, 40.0, 1.40, 1.0, 0.0)"
        )
        conn.execute(
            "INSERT INTO HWSD2_LAYERS VALUES (1, 1, 100, 15, 30, 38.0, 22.0, 40.0, 1.45, 0.8, 1.0)"
        )
        conn.execute(
            "INSERT INTO HWSD2_LAYERS VALUES (1, 1, 100, 30, 60, 35.0, 25.0, 40.0, 1.50, 0.5, 2.0)"
        )
        conn.commit()
        conn.close()

        out_csv = os.path.join(self.work, "soil_map.csv")
        out_sol_dir = os.path.join(self.work, "individual_sol")

        process_soils_hwsd(
            grid_points=self.shapefile,
            hwsd_raster_file=dummy_raster,
            hwsd_db_file=dummy_db,
            output_csv_path=out_csv,
            output_sol_dir=out_sol_dir,
            id_col="ID",
            lat_col="LAT",
            long_col="LONG",
        )

        self.assertTrue(os.path.exists(out_csv))
        self._assert_sol_valid(os.path.join(out_sol_dir, "TEST_1.SOL"))


if __name__ == "__main__":
    unittest.main()
