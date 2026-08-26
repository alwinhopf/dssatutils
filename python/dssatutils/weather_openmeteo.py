# File: weather_openmeteo.py
# ---------------------------------------------------------------------------
# Open-Meteo Historical Weather (ERA5 / ERA5-Land reanalysis) -> DSSAT .WTH.
#
# WHY THIS SOURCE: gives truly GLOBAL daily coverage (Europe, Asia, Africa,
# Oceania, South America) from 1940 onward with NO API KEY and no registration,
# complementing DAYMET (North America only) and GRIDMET (US only). NASA POWER is
# also global; ERA5-Seamless preserves ERA5-Land temperature/humidity while
# supplying the complete forcing record DSSAT needs.
#
# API docs: https://open-meteo.com/en/docs/historical-weather-api
# License:  data is CC-BY 4.0 (ERA5 by Copernicus/ECMWF). Cite when publishing.
#
# Mirrors the process_weather_nasapower() signature exactly so it is a drop-in
# WEATHER_SOURCE for the pipeline.
# ---------------------------------------------------------------------------

from __future__ import annotations

import os
import math
import time
import logging
from datetime import date, timedelta
# HTTP-bound work is better served by threads; the legacy local name remains a
# stable monkeypatch seam for downstream offline tests.
from concurrent.futures import ThreadPoolExecutor as ProcessPoolExecutor, as_completed

import pandas as pd
import requests

from .config import get_config_number, get_config_value

logger = logging.getLogger(__name__)

_ARCHIVE_URL = get_config_value(
    "weather.openmeteo.archive_url",
    "https://archive-api.open-meteo.com/v1/archive",
)
_OPEN_METEO_MODEL = get_config_value(
    "weather.openmeteo.model",
    "era5_seamless",
)

# Daily variables requested from Open-Meteo (ERA5 archive).
_DAILY_VARS = [
    "temperature_2m_max",          # degC
    "temperature_2m_min",          # degC
    "precipitation_sum",           # mm
    "shortwave_radiation_sum",     # MJ/m2  (DSSAT SRAD)
    "wind_speed_10m_mean",         # m/s daily mean (windspeed_unit=ms)
    "dew_point_2m_mean",           # degC daily mean (DSSAT TDEW)
    "relative_humidity_2m_mean",   # % daily mean (DSSAT RH2M)
]

# Log-wind-profile factor to convert 10 m wind to 2 m (FAO-56):
#   u2 = u10 * 4.87 / ln(67.8*10 - 5.42)  ->  ~0.748
_WIND_10M_TO_2M = 0.748


def _calc_tav(df: pd.DataFrame) -> float:
    return float(((df["TMAX"] + df["TMIN"]) / 2.0).mean())


def _calc_amp(df: pd.DataFrame) -> float:
    df = df.copy()
    df["TAVG"] = (df["TMAX"] + df["TMIN"]) / 2.0
    monthly = df.groupby(["YEAR", "MM"])["TAVG"].mean().reset_index()
    annual = monthly.groupby("YEAR")["TAVG"].agg(lambda x: x.max() - x.min())
    return float(annual.mean())


def _fetch_open_meteo(lat: float, lon: float, start: str, end: str,
                      retries: int | None = None,
                      backoff: float | None = None) -> pd.DataFrame:
    """Fetch Open-Meteo daily archive for one point. start/end: 'YYYY-MM-DD'."""
    retries = int(retries if retries is not None else get_config_number(
        "weather.openmeteo.fetch_retries", 4
    ))
    backoff = float(backoff if backoff is not None else get_config_number(
        "weather.openmeteo.fetch_backoff_seconds", 5
    ))
    params = {
        "latitude": lat,
        "longitude": lon,
        "start_date": start,
        "end_date": end,
        "daily": ",".join(_DAILY_VARS),
        "windspeed_unit": "ms",
        "timezone": "UTC",
        # ERA5-Seamless uses ERA5-Land for temperature/humidity and ERA5 for
        # forcing variables that ERA5-Land does not expose through Open-Meteo
        # (notably radiation and wind). This produces a complete DSSAT record.
        "models": _OPEN_METEO_MODEL,
    }
    for attempt in range(retries):
        try:
            r = requests.get(_ARCHIVE_URL, params=params, timeout=180)
            # Open-Meteo returns 429 when rate-limited — back off and retry.
            if r.status_code == 429:
                raise requests.HTTPError("429 rate-limited")
            r.raise_for_status()
            daily = r.json()["daily"]
            df = pd.DataFrame(daily)
            df["time"] = pd.to_datetime(df["time"])
            df["YEAR"] = df["time"].dt.year
            df["MM"] = df["time"].dt.month
            df["DOY"] = df["time"].dt.day_of_year
            df["DATE"] = df["YEAR"].astype(str) + df["DOY"].astype(str).str.zfill(3)
            return df
        except Exception as exc:  # noqa: BLE001
            if attempt < retries - 1:
                time.sleep(backoff * (attempt + 1))
            else:
                raise RuntimeError(f"Open-Meteo fetch failed: {exc}") from exc


def _process_single_point(args: dict) -> None:
    lat = args["latitude"]
    lon = args["longitude"]
    pid = args["point_id"]
    output_dir = args["output_dir"]
    start_date = args["start_date"]
    end_date = args["end_date"]
    log_file = args["log_file"]

    out_path = os.path.join(output_dir, f"{pid}.WTH")
    if os.path.exists(out_path):
        return

    try:
        df = _fetch_open_meteo(lat, lon, start_date, end_date)
        if df.empty:
            raise ValueError("No data returned from Open-Meteo.")

        df = df.rename(columns={
            "shortwave_radiation_sum": "SRAD",
            "temperature_2m_max": "TMAX",
            "temperature_2m_min": "TMIN",
            "precipitation_sum": "RAIN",
            "wind_speed_10m_mean": "WIND",
            "dew_point_2m_mean": "TDEW",
            "relative_humidity_2m_mean": "RH2M",
        })
        required = ["SRAD", "TMAX", "TMIN", "RAIN", "WIND", "TDEW", "RH2M"]
        missing = [col for col in required if col not in df.columns]
        empty = [col for col in required if col in df.columns and df[col].isna().all()]
        if missing or empty:
            details = []
            if missing:
                details.append("missing columns: " + ", ".join(missing))
            if empty:
                details.append("all-null columns: " + ", ".join(empty))
            raise ValueError(
                "Open-Meteo did not return complete daily DSSAT weather ("
                + "; ".join(details) + ")."
            )
        # 10 m -> 2 m wind adjustment.
        df["WIND"] = df["WIND"] * _WIND_10M_TO_2M
        # Fill any gaps with the DSSAT missing sentinel.
        for col in required:
            df[col] = df[col].fillna(-99.0)

        tav = _calc_tav(df)
        amp = _calc_amp(df)

        header = (
            f"$WEATHER DATA: OPEN-METEO ERA5-SEAMLESS (Point ID: {pid})\n"
            f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
            f"  OMET {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}   2.0   2.0\n"
            f"@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND"
        )
        lines = []
        for _, row in df.iterrows():
            line = (
                f"{row['DATE']:>7s}"
                f"{row['SRAD']:6.1f}{row['TMAX']:6.1f}{row['TMIN']:6.1f}"
                f"{row['RAIN']:6.1f}{row['TDEW']:6.1f}{row['RH2M']:6.1f}{row['WIND']:6.1f}"
            )
            line = line.replace(" -99.0", "   -99")
            lines.append(line)

        with open(out_path, "w") as fh:
            fh.write(header + "\n")
            fh.write("\n".join(lines) + "\n")

    except Exception as exc:  # noqa: BLE001
        msg = (f"\n--- ERROR ---\nFailed: Point ID {pid} | "
               f"Lat {lat:.3f}, Lon {lon:.3f}\nError: {exc}\n")
        print(msg)
        with open(log_file, "a") as lf:
            lf.write(msg)


def process_weather_openmeteo(shapefile, start_year, end_year, output_dir,
                              id_col, lat_col, lon_col, n_cores, log_file):
    """Download Open-Meteo ERA5-Seamless daily weather and write DSSAT .WTH.

    ERA5-Seamless combines ERA5-Land temperature/humidity with the ERA5
    forcing variables needed for complete DSSAT weather. Drop-in replacement
    for process_weather_nasapower(); same signature.
    """
    today = date.today()
    start_date = f"{start_year}-01-01"
    if end_year >= today.year:
        archive_lag_days = int(get_config_number("weather.openmeteo.archive_lag_days", 6))
        safe_end = today - timedelta(days=archive_lag_days)
        end_date = safe_end.isoformat()
        print(f"End year is current/future. Fetching up to: {end_date}")
    else:
        end_date = f"{end_year}-12-31"

    print(f"--- Starting Open-Meteo ({_OPEN_METEO_MODEL}) Download (Years: {start_year}-{end_year}) ---")
    print(f"Registered {n_cores} cores for parallel Open-Meteo download.")
    os.makedirs(output_dir, exist_ok=True)

    tasks = [
        dict(latitude=float(row[lat_col]), longitude=float(row[lon_col]),
             point_id=str(row[id_col]), output_dir=output_dir,
             start_date=start_date, end_date=end_date, log_file=log_file)
        for _, row in shapefile.iterrows()
    ]

    with ProcessPoolExecutor(max_workers=n_cores) as pool:
        futures = {pool.submit(_process_single_point, t): t["point_id"] for t in tasks}
        for fut in as_completed(futures):
            pid = futures[fut]
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                print(f"ERROR (point {pid}): {exc}")

    print(f"\nOpen-Meteo processing complete. Check the '{output_dir}' directory.\n")
