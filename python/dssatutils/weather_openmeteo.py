# File: weather_openmeteo.py
# ---------------------------------------------------------------------------
# Open-Meteo Historical Weather (ERA5 / ERA5-Land reanalysis) -> DSSAT .WTH.
#
# WHY THIS SOURCE: gives truly GLOBAL daily coverage (Europe, Asia, Africa,
# Oceania, South America) from 1940 onward with NO API KEY and no registration,
# complementing DAYMET (North America only) and GRIDMET (US only). NASA POWER is
# also global; Open-Meteo (ERA5-Land, ~9 km) is a higher-resolution alternative
# for non-US regions.
#
# API docs: https://open-meteo.com/en/docs/historical-weather-api
# License:  data is CC-BY 4.0 (ERA5 by Copernicus/ECMWF). Cite when publishing.
#
# Mirrors the process_weather_nasapower() signature exactly so it is a drop-in
# WEATHER_SOURCE for the pipeline.
# ---------------------------------------------------------------------------

import os
import math
import time
import logging
from datetime import date, timedelta
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import requests

logger = logging.getLogger(__name__)

_ARCHIVE_URL = "https://archive-api.open-meteo.com/v1/archive"

# Daily variables requested from Open-Meteo (ERA5 archive).
_DAILY_VARS = [
    "temperature_2m_max",          # degC
    "temperature_2m_min",          # degC
    "precipitation_sum",           # mm
    "shortwave_radiation_sum",     # MJ/m2  (DSSAT SRAD)
    "wind_speed_10m_max",          # m/s (windspeed_unit=ms)
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
                      retries: int = 4, backoff: float = 5.0) -> pd.DataFrame:
    """Fetch Open-Meteo daily archive for one point. start/end: 'YYYY-MM-DD'."""
    params = {
        "latitude": lat,
        "longitude": lon,
        "start_date": start,
        "end_date": end,
        "daily": ",".join(_DAILY_VARS),
        "windspeed_unit": "ms",
        "timezone": "UTC",
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
            "wind_speed_10m_max": "WIND",
        })
        # 10 m -> 2 m wind adjustment.
        df["WIND"] = df["WIND"] * _WIND_10M_TO_2M
        # Open-Meteo does not provide daily dewpoint / RH -> mark missing (-99).
        df["TDEW"] = -99.0
        df["RH2M"] = -99.0
        # Fill any gaps with the DSSAT missing sentinel.
        for col in ["SRAD", "TMAX", "TMIN", "RAIN", "WIND"]:
            df[col] = df[col].fillna(-99.0)

        tav = _calc_tav(df)
        amp = _calc_amp(df)

        header = (
            f"$WEATHER DATA: OPEN-METEO ERA5 (Point ID: {pid})\n"
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
    """Download Open-Meteo (ERA5) daily weather for every point and write .WTH.

    Drop-in replacement for process_weather_nasapower(); same signature.
    """
    today = date.today()
    start_date = f"{start_year}-01-01"
    if end_year >= today.year:
        # ERA5 archive lags ~5 days.
        safe_end = today - timedelta(days=6)
        end_date = safe_end.isoformat()
        print(f"End year is current/future. Fetching up to: {end_date}")
    else:
        end_date = f"{end_year}-12-31"

    print(f"--- Starting Open-Meteo (ERA5) Download (Years: {start_year}-{end_year}) ---")
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
