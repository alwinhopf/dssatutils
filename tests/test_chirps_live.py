#!/usr/bin/env python3
"""Live CHIRPS v2/v3 integration checks.

These tests intentionally download real CHIRPS NetCDFs and call NASA POWER.
They are skipped by default; run explicitly with:

    DSSATUTILS_RUN_LIVE_CHIRPS=1 python -m pytest tests/test_chirps_live.py -m live -q -s
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from dssatutils import process_weather_nasapower_chirps, process_weather_nasapower_chirps_v3


pytestmark = pytest.mark.live


def _live_enabled() -> bool:
    return os.getenv("DSSATUTILS_RUN_LIVE_CHIRPS", "").strip().lower() in {"1", "true", "yes"}


def _read_wth(path: Path) -> tuple[str, pd.DataFrame]:
    text = path.read_text()
    rows = [line.split() for line in text.splitlines() if line[:4].isdigit()]
    df = pd.DataFrame(
        rows,
        columns=["DATE", "SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND"],
    )
    for col in df.columns[1:]:
        df[col] = pd.to_numeric(df[col])
    return text, df


@pytest.mark.skipif(not _live_enabled(), reason="set DSSATUTILS_RUN_LIVE_CHIRPS=1 to run live CHIRPS downloads")
def test_live_nasapower_chirps_v2_v3_real_data_comparison(tmp_path):
    """Generate real DSSAT weather files from CHIRPS v2 and CHIRPS v3 rainfall."""
    point = pd.DataFrame({"ID": ["UG_KLA"], "LAT": [0.3476], "LONG": [32.5825]})
    cache_root = Path(os.getenv("DSSATUTILS_LIVE_CACHE", Path.cwd() / ".live_cache"))
    cache_root.mkdir(parents=True, exist_ok=True)

    out_v2 = tmp_path / "wth_v2"
    out_v3 = tmp_path / "wth_v3"

    import dssatutils.weather_nasapower_chirps as chirps_v2

    old_res = chirps_v2.CHIRPS_RESOLUTION
    chirps_v2.CHIRPS_RESOLUTION = "p25"
    try:
        process_weather_nasapower_chirps(
            shapefile=point,
            start_year=2010,
            end_year=2010,
            output_dir=str(out_v2),
            id_col="ID",
            lat_col="LAT",
            lon_col="LONG",
            n_cores=1,
            log_file=str(tmp_path / "v2_errors.log"),
            chirps_cache_dir=str(cache_root / "chirps_v2_cache"),
        )
    finally:
        chirps_v2.CHIRPS_RESOLUTION = old_res

    process_weather_nasapower_chirps_v3(
        shapefile=point,
        start_year=2010,
        end_year=2010,
        output_dir=str(out_v3),
        id_col="ID",
        lat_col="LAT",
        lon_col="LONG",
        n_cores=1,
        log_file=str(tmp_path / "v3_errors.log"),
        chirps_cache_dir=str(cache_root / "chirps_v3_cache"),
        chirps_product="rnl",
        chirps_stream="final",
        chirps_fetch_mode="monthly_netcdf",
        chirps_months=[3],
    )

    text_v2, df_v2 = _read_wth(out_v2 / "UG_KLA.WTH")
    text_v3, df_v3 = _read_wth(out_v3 / "UG_KLA.WTH")
    assert "NASA-POWER + CHIRPS rain" in text_v2
    assert "CHIRPS(p25) where available, 365 days" in text_v2
    assert "NASA-POWER + CHIRPS-v3 rain" in text_v3
    assert "CHIRPS-v3 final/rnl where available, 31 days" in text_v3
    assert len(df_v2) == 365
    assert len(df_v3) == 365

    march = [f"2010{d:03d}" for d in range(60, 91)]
    merged = df_v2[df_v2["DATE"].isin(march)][["DATE", "RAIN"]].merge(
        df_v3[df_v3["DATE"].isin(march)][["DATE", "RAIN"]],
        on="DATE",
        suffixes=("_v2", "_v3"),
    )
    assert len(merged) == 31
    assert np.isfinite(merged[["RAIN_v2", "RAIN_v3"]].to_numpy()).all()
    assert merged["RAIN_v2"].sum() >= 0
    assert merged["RAIN_v3"].sum() >= 0
    assert (merged["RAIN_v2"] - merged["RAIN_v3"]).abs().sum() > 0

    v2_files = list((cache_root / "chirps_v2_cache").glob("chirps-v2.0.2010.days_p25.nc"))
    v3_files = list((cache_root / "chirps_v3_cache").glob("v3_final_rnl_monthly_netcdf/chirps-v3.0.2010.03.days_p05.nc"))
    assert v2_files and v2_files[0].stat().st_size > 10_000_000
    assert v3_files and v3_files[0].stat().st_size > 100_000_000


@pytest.mark.skipif(not _live_enabled(), reason="set DSSATUTILS_RUN_LIVE_CHIRPS=1 to run live CHIRPS downloads")
def test_live_remote_cog(tmp_path):
    """Test CHIRPS v3 extraction via remote_cog fetch mode."""
    point = pd.DataFrame({"ID": ["UG_KLA"], "LAT": [0.3476], "LONG": [32.5825]})
    cache_root = Path(os.getenv("DSSATUTILS_LIVE_CACHE", Path.cwd() / ".live_cache"))
    cache_root.mkdir(parents=True, exist_ok=True)

    out_nc = tmp_path / "wth_nc"
    out_cog = tmp_path / "wth_cog"

    # Fetch using monthly_netcdf (only March 2010)
    process_weather_nasapower_chirps_v3(
        shapefile=point,
        start_year=2010,
        end_year=2010,
        output_dir=str(out_nc),
        id_col="ID",
        lat_col="LAT",
        lon_col="LONG",
        n_cores=1,
        log_file=str(tmp_path / "nc_errors.log"),
        chirps_cache_dir=str(cache_root / "chirps_v3_cache_nc"),
        chirps_product="rnl",
        chirps_stream="final",
        chirps_fetch_mode="monthly_netcdf",
        chirps_months=[3],
    )

    # Fetch using remote_cog (only March 2010)
    process_weather_nasapower_chirps_v3(
        shapefile=point,
        start_year=2010,
        end_year=2010,
        output_dir=str(out_cog),
        id_col="ID",
        lat_col="LAT",
        lon_col="LONG",
        n_cores=1,
        log_file=str(tmp_path / "cog_errors.log"),
        chirps_cache_dir=str(cache_root / "chirps_v3_cache_cog"),
        chirps_product="rnl",
        chirps_stream="final",
        chirps_fetch_mode="remote_cog",
        chirps_months=[3],
    )

    text_nc, df_nc = _read_wth(out_nc / "UG_KLA.WTH")
    text_cog, df_cog = _read_wth(out_cog / "UG_KLA.WTH")

    march = [f"2010{d:03d}" for d in range(60, 91)]
    df_nc_march = df_nc[df_nc["DATE"].isin(march)].reset_index(drop=True)
    df_cog_march = df_cog[df_cog["DATE"].isin(march)].reset_index(drop=True)

    assert len(df_nc_march) == 31
    assert len(df_cog_march) == 31
    np.testing.assert_allclose(df_nc_march["RAIN"], df_cog_march["RAIN"], atol=1e-3)


@pytest.mark.skipif(not _live_enabled(), reason="set DSSATUTILS_RUN_LIVE_CHIRPS=1 to run live CHIRPS downloads")
def test_live_gee(tmp_path):
    """Test CHIRPS v3 extraction via GEE fetch mode."""
    try:
        from dssatutils.weather_chirps_v3 import _init_gee
        _init_gee()
    except Exception as e:
        pytest.skip(f"GEE not initialized/configured: {e}")

    point = pd.DataFrame({"ID": ["UG_KLA"], "LAT": [0.3476], "LONG": [32.5825]})
    cache_root = Path(os.getenv("DSSATUTILS_LIVE_CACHE", Path.cwd() / ".live_cache"))
    cache_root.mkdir(parents=True, exist_ok=True)

    out_nc = tmp_path / "wth_nc"
    out_gee = tmp_path / "wth_gee"

    # Fetch using monthly_netcdf (only March 2010)
    process_weather_nasapower_chirps_v3(
        shapefile=point,
        start_year=2010,
        end_year=2010,
        output_dir=str(out_nc),
        id_col="ID",
        lat_col="LAT",
        lon_col="LONG",
        n_cores=1,
        log_file=str(tmp_path / "nc_errors.log"),
        chirps_cache_dir=str(cache_root / "chirps_v3_cache_nc"),
        chirps_product="rnl",
        chirps_stream="final",
        chirps_fetch_mode="monthly_netcdf",
        chirps_months=[3],
    )

    # Fetch using GEE (only March 2010)
    process_weather_nasapower_chirps_v3(
        shapefile=point,
        start_year=2010,
        end_year=2010,
        output_dir=str(out_gee),
        id_col="ID",
        lat_col="LAT",
        lon_col="LONG",
        n_cores=1,
        log_file=str(tmp_path / "gee_errors.log"),
        chirps_cache_dir=str(cache_root / "chirps_v3_cache_gee"),
        chirps_product="rnl",
        chirps_stream="final",
        chirps_fetch_mode="gee",
        chirps_months=[3],
    )

    text_nc, df_nc = _read_wth(out_nc / "UG_KLA.WTH")
    text_gee, df_gee = _read_wth(out_gee / "UG_KLA.WTH")

    march = [f"2010{d:03d}" for d in range(60, 91)]
    df_nc_march = df_nc[df_nc["DATE"].isin(march)].reset_index(drop=True)
    df_gee_march = df_gee[df_gee["DATE"].isin(march)].reset_index(drop=True)

    assert len(df_nc_march) == 31
    assert len(df_gee_march) == 31
    np.testing.assert_allclose(df_nc_march["RAIN"], df_gee_march["RAIN"], atol=1e-3)

