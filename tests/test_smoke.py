#!/usr/bin/env python3
# File: tests/test_smoke.py
# ---------------------------------------------------------------------------
# OFFLINE smoke test — runs in CI with no internet and no DSSAT install.
#
# Validates:
#   1. All helper modules import cleanly from the package.
#   2. The DSSAT .WTH writer produces a well-formed file. We monkeypatch the
#      Open-Meteo network call to return a synthetic 2-year daily frame, so the
#      formatting / TAV / AMP / column-width logic is tested deterministically.
#
# Run:  python tests/test_smoke.py     (exit 0 = pass)
# ---------------------------------------------------------------------------

import os
import sys
import tempfile
import shutil
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
# Insert the 'python' folder where 'dssatutils' package resides
sys.path.insert(0, os.path.join(_REPO, "python"))

# ---------------------------------------------------------------------------
# pytest tests
# ---------------------------------------------------------------------------

def test_module_imports():
    # Verify we can import all modules from the package
    modules = [
        "dssatutils.weather_daymet",
        "dssatutils.weather_nasapower",
        "dssatutils.weather_gridmet",
        "dssatutils.weather_openmeteo",
        "dssatutils.weather_nasapower_chirps",
        "dssatutils.weather_chirps_v3",
        "dssatutils.weather_rainfall_merge",
        "dssatutils.weather_agera5",
        "dssatutils.soil_hwsd",
        "dssatutils.soil_soilgrids",
        "dssatutils.soil_soilgrids_online",
        "dssatutils.soil_ssurgo",
        "dssatutils.soil_ssurgo_alderman",
        "dssatutils.weather_era5land",
    ]
    for mod in modules:
        try:
            __import__(mod)
        except ImportError as exc:
            # Skip if optional heavy dep is missing, raise only if it's not a ModuleNotFoundError
            if "ModuleNotFoundError" in str(type(exc)):
                print(f"  [info] Optional dep missing for {mod}: {exc}")
            else:
                raise

def _make_fake_fetch():
    def _fake_fetch(lat, lon, start, end, retries=4, backoff=5.0):
        dates = pd.date_range(start, end, freq="D")
        n = len(dates)
        df = pd.DataFrame({
            "time": dates,
            "temperature_2m_max": [20.0 + (i % 10) for i in range(n)],
            "temperature_2m_min": [8.0 + (i % 5) for i in range(n)],
            "precipitation_sum": [0.0 if i % 3 else 5.0 for i in range(n)],
            "shortwave_radiation_sum": [18.0 for _ in range(n)],
            "wind_speed_10m_max": [3.0 for _ in range(n)],
        })
        df["YEAR"] = df["time"].dt.year
        df["MM"] = df["time"].dt.month
        df["DOY"] = df["time"].dt.day_of_year
        df["DATE"] = df["YEAR"].astype(str) + df["DOY"].astype(str).str.zfill(3)
        return df
    return _fake_fetch


def test_wth_writer_synthetic():
    from dssatutils import weather_openmeteo as om
    om._fetch_open_meteo = _make_fake_fetch()

    with tempfile.TemporaryDirectory() as work:
        log = os.path.join(work, "errors.log")
        om._process_single_point(dict(
            latitude=52.0, longitude=5.0, point_id="EU_TEST",
            output_dir=work, start_date="2010-01-01", end_date="2011-12-31",
            log_file=log,
        ))
        out = os.path.join(work, "EU_TEST.WTH")
        assert os.path.exists(out), "EU_TEST.WTH not written"

        lines = [ln.rstrip("\n") for ln in open(out) if ln.strip()]
        assert lines[0].startswith("$WEATHER"), "title line missing"
        assert lines[1].lstrip().startswith("@ INSI"), "INSI header missing"
        assert lines[3].lstrip().startswith("@  DATE"), "DATE column header missing"

        data = lines[4:]
        assert len(data) == 730, f"expected 730 daily rows, got {len(data)}"
        assert data[0][:7].strip() == "2010001", \
            f"first DATE token wrong: {data[0][:7].strip()!r}"
        assert "nan" not in "".join(data).lower(), "NaN found in data block"


def test_gridmet_amp_matches_dssat_definition():
    from dssatutils.weather_gridmet import _calc_amp

    dates = pd.date_range("2001-01-01", "2002-12-31", freq="D")
    monthly_mean = dates.month.to_numpy(dtype=float)
    tmax = monthly_mean + 5.0
    tmin = monthly_mean - 5.0

    # Calendar-month means span 1..12 C; DSSAT AMP is half that range.
    assert _calc_amp(tmax, tmin, dates) == 5.5


# ---------------------------------------------------------------------------
# Standalone script mode
# ---------------------------------------------------------------------------

def _run_standalone():
    _fail = 0

    def check(cond, msg):
        nonlocal _fail
        if cond:
            print(f"  [ok] {msg}")
        else:
            print(f"  [FAIL] {msg}")
            _fail += 1

    print("\n[1/2] package module imports...")
    modules = [
        "dssatutils.weather_daymet",
        "dssatutils.weather_nasapower",
        "dssatutils.weather_gridmet",
        "dssatutils.weather_openmeteo",
        "dssatutils.weather_nasapower_chirps",
        "dssatutils.weather_chirps_v3",
        "dssatutils.weather_rainfall_merge",
        "dssatutils.weather_agera5",
        "dssatutils.soil_hwsd",
        "dssatutils.soil_soilgrids",
        "dssatutils.soil_soilgrids_online",
        "dssatutils.soil_ssurgo",
        "dssatutils.soil_ssurgo_alderman",
        "dssatutils.weather_era5land",
    ]
    for mod in modules:
        try:
            __import__(mod)
            check(True, f"import {mod}")
        except ImportError as exc:
            print(f"  [skip] import {mod}: optional dep missing ({exc})")

    print("\n[2/2] Open-Meteo .WTH writer (synthetic, no network)...")
    from dssatutils import weather_openmeteo as om
    om._fetch_open_meteo = _make_fake_fetch()

    with tempfile.TemporaryDirectory() as work:
        log = os.path.join(work, "errors.log")
        om._process_single_point(dict(
            latitude=52.0, longitude=5.0, point_id="EU_TEST",
            output_dir=work, start_date="2010-01-01", end_date="2011-12-31",
            log_file=log,
        ))
        out = os.path.join(work, "EU_TEST.WTH")
        check(os.path.exists(out), "EU_TEST.WTH written")
        if os.path.exists(out):
            lines = [ln.rstrip("\n") for ln in open(out) if ln.strip()]
            check(lines[0].startswith("$WEATHER"), "title line present")
            check(lines[1].lstrip().startswith("@ INSI"), "INSI header present")
            check(lines[3].lstrip().startswith("@  DATE"), "DATE column header present")
            data = lines[4:]
            check(len(data) == 730, f"730 daily rows (got {len(data)})")
            check(data[0][:7].strip() == "2010001", f"first DATE token = {data[0][:7].strip()}")
            check("nan" not in "".join(data).lower(), "no NaN in data block")

    print("\n=== RESULT ===")
    if _fail == 0:
        print("ALL SMOKE CHECKS PASSED.")
        sys.exit(0)
    else:
        print(f"{_fail} CHECK(S) FAILED.")
        sys.exit(1)


if __name__ == "__main__":
    _run_standalone()
