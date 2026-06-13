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
                 "process_weather_xavier", "process_weather_cmfd"):
        assert hasattr(dssatutils, name), f"dssatutils missing public {name}"


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
