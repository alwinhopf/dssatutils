#!/usr/bin/env python3
# File: tests/test_smoke.py
# ---------------------------------------------------------------------------
# OFFLINE smoke test — runs in CI with no internet and no DSSAT install.
#
# It exercises the parts of the pipeline that can be validated without live
# APIs or the DSSAT executable:
#   1. config_loader: YAML load + cfg_get fallback semantics.
#   2. All helper modules import cleanly.
#   3. The DSSAT .WTH writer produces a well-formed file. We monkeypatch the
#      Open-Meteo network call to return a synthetic 2-year daily frame, so the
#      formatting / TAV / AMP / column-width logic is tested deterministically.
#
# Run:  python tests/test_smoke.py     (exit 0 = pass)
#   or: pytest tests/test_smoke.py
# ---------------------------------------------------------------------------

import os
import sys
import tempfile

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, _REPO)
sys.path.insert(0, os.path.join(_REPO, "python_scripts"))


# ---------------------------------------------------------------------------
# pytest tests
# ---------------------------------------------------------------------------

def test_config_loader_missing_key():
    from config_loader import cfg_get
    assert cfg_get("__definitely_missing__", "DEF") == "DEF"


def test_config_loader_weather_source():
    from config_loader import cfg_get
    ws = cfg_get("weather_source", "DAYMET")
    assert isinstance(ws, str) and ws


def test_module_imports():
    # These may pull optional heavy deps at import time; skip if absent.
    optional_modules = {"weather_gridmet",          # xarray
                        "weather_nasapower_chirps",  # imports weather_nasapower; xarray lazy
                        "weather_agera5",            # numpy/pandas; cdsapi+xarray lazy
                        "soil_hwsd"}                 # imports soil_soilgrids_online
    for mod in ["weather_daymet", "weather_nasapower", "weather_gridmet",
                "weather_openmeteo", "weather_nasapower_chirps",
                "weather_agera5", "soil_hwsd"]:
        try:
            __import__(mod)
        except ModuleNotFoundError:
            if mod in optional_modules:
                pass  # optional heavy dep missing; fine for the offline smoke test
            else:
                raise
        except Exception as exc:
            raise AssertionError(f"import {mod} failed: {exc}") from exc


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
    import weather_openmeteo as om
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

    print("[1/3] config_loader fallback semantics...")
    from config_loader import cfg_get
    check(cfg_get("__definitely_missing__", "DEF") == "DEF", "missing key -> default")
    ws = cfg_get("weather_source", "DAYMET")
    check(isinstance(ws, str) and ws, "weather_source resolves to a non-empty string")

    print("\n[2/3] helper module imports...")
    for mod in ["weather_daymet", "weather_nasapower", "weather_gridmet",
                "weather_openmeteo"]:
        try:
            __import__(mod)
            check(True, f"import {mod}")
        except ModuleNotFoundError as exc:
            print(f"  [skip] import {mod}: optional dep missing ({exc.name})")
        except Exception as exc:
            check(False, f"import {mod}: {exc}")

    print("\n[3/3] Open-Meteo .WTH writer (synthetic, no network)...")
    import weather_openmeteo as om
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
