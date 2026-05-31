# File: weather_nasapower.py
# Python port of weather_nasapower.R
#
# Fetches NASA POWER daily data for each grid point and writes
# DSSAT-formatted .WTH files.  Covers the globe; available from ~1981 onward.
#
# API docs: https://power.larc.nasa.gov/docs/services/api/temporal/daily/

import os
import time
import logging
from datetime import date
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests

logger = logging.getLogger(__name__)

# NASA POWER parameters required for DSSAT
_NASA_PARAMS = [
    "T2M_MAX",            # Max air temperature at 2 m  (°C)
    "T2M_MIN",            # Min air temperature at 2 m  (°C)
    "ALLSKY_SFC_SW_DWN",  # Solar radiation (MJ/m²/day)
    "PRECTOTCORR",        # Precipitation (mm/day)
    "T2MDEW",             # Dew-point temperature (°C)
    "RH2M",               # Relative humidity at 2 m (%)
    "WS2M",               # Wind speed at 2 m (m/s)
]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _calc_tav(df: pd.DataFrame) -> float:
    tavg = (df["TMAX"] + df["TMIN"]) / 2.0
    return float(tavg.mean())


def _calc_amp(df: pd.DataFrame) -> float:
    df = df.copy()
    df["TAVG"] = (df["TMAX"] + df["TMIN"]) / 2.0
    monthly = df.groupby(["YEAR", "MM"])["TAVG"].mean().reset_index()
    annual = monthly.groupby("YEAR")["TAVG"].agg(lambda x: x.max() - x.min())
    return float(annual.mean())


def _fetch_nasa_power(lat: float, lon: float, start: str, end: str,
                      retries: int = 3, backoff: float = 5.0) -> pd.DataFrame:
    """
    Fetch NASA POWER daily data for a single point.
    start / end: "YYYYMMDD" strings.
    Returns a DataFrame with columns matching _NASA_PARAMS plus YEAR, MM, DOY.
    """
    url = "https://power.larc.nasa.gov/api/temporal/daily/point"
    params = {
        "parameters": ",".join(_NASA_PARAMS),
        "community": "AG",
        "longitude": lon,
        "latitude": lat,
        "start": start,
        "end": end,
        "format": "JSON",
    }
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=180)
            r.raise_for_status()
            data = r.json()
            param_dict = data["properties"]["parameter"]
            # param_dict[param][YYYYMMDD] = value
            records = {}
            for param, daily in param_dict.items():
                records[param] = daily
            df = pd.DataFrame(records)
            df.index = pd.to_datetime(df.index, format="%Y%m%d")
            df["YEAR"] = df.index.year
            df["MM"] = df.index.month
            df["DOY"] = df.index.day_of_year
            df["DATE"] = df["YEAR"].astype(str) + df["DOY"].astype(str).str.zfill(3)
            df = df.reset_index(drop=True)
            # Replace NASA missing-value sentinel (-999) with -99
            df = df.replace(-999, -99)
            return df
        except Exception as exc:
            if attempt < retries - 1:
                time.sleep(backoff * (attempt + 1))
            else:
                raise RuntimeError(f"NASA POWER fetch failed: {exc}") from exc


def _process_single_point(args: dict) -> None:
    """Worker function executed in a subprocess."""
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
        df = _fetch_nasa_power(lat, lon, start_date, end_date)

        if df.empty:
            raise ValueError("No data returned from NASA POWER.")

        # Rename parameter columns to DSSAT convention
        rename = {
            "ALLSKY_SFC_SW_DWN": "SRAD",
            "T2M_MAX": "TMAX",
            "T2M_MIN": "TMIN",
            "PRECTOTCORR": "RAIN",
            "T2MDEW": "TDEW",
            "RH2M": "RH2M",
            "WS2M": "WIND",
        }
        df = df.rename(columns=rename)

        tav = _calc_tav(df)
        amp = _calc_amp(df)

        # --- Write WTH file ---
        header = (
            f"$WEATHER DATA: NASA-POWER (Point ID: {pid})\n"
            f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
            f"  NASA {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}   2.0   2.0\n"
            f"@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND"
        )

        lines = []
        for _, row in df.iterrows():
            line = (
                f"{row['DATE']:>7s}"
                f"{row['SRAD']:6.1f}"
                f"{row['TMAX']:6.1f}"
                f"{row['TMIN']:6.1f}"
                f"{row['RAIN']:6.1f}"
                f"{row['TDEW']:6.1f}"
                f"{row['RH2M']:6.1f}"
                f"{row['WIND']:6.1f}"
            )
            line = line.replace(" -99.0", "   -99")
            lines.append(line)

        with open(out_path, "w") as fh:
            fh.write(header + "\n")
            fh.write("\n".join(lines) + "\n")

    except Exception as exc:
        msg = (
            f"\n--- ERROR ---\n"
            f"Failed: Point ID {pid} | Lat {lat:.3f}, Lon {lon:.3f}\n"
            f"Error: {exc}\n"
        )
        print(msg)
        with open(log_file, "a") as lf:
            lf.write(msg)


# ---------------------------------------------------------------------------
# Public entry point (mirrors R function signature exactly)
# ---------------------------------------------------------------------------

def process_weather_nasapower(
    shapefile,           # GeoDataFrame
    start_year: int,
    end_year: int,
    output_dir: str,
    id_col: str,
    lat_col: str,
    lon_col: str,
    n_cores: int,
    log_file: str,
) -> None:
    """
    Download NASA POWER weather data for every point in *shapefile* and write
    DSSAT .WTH files to *output_dir*.  Mirrors the R ``process_weather_nasapower``
    function, including the current-year partial-data handling.
    """
    today = date.today()
    current_year = today.year
    start_date = f"{start_year}0101"
    if end_year == current_year:
        # NASA POWER data typically lags by a day or two
        safe_end = today - __import__("datetime").timedelta(days=2)
        end_date = safe_end.strftime("%Y%m%d")
        print(f"End year is current year. Fetching data up to: {safe_end.isoformat()}")
    else:
        end_date = f"{end_year}1231"

    print(f"--- Starting NASA-POWER Download (Years: {start_year}–{end_year}) ---")
    print(f"Registered {n_cores} cores for parallel NASA-POWER download.")

    os.makedirs(output_dir, exist_ok=True)

    tasks = []
    for _, row in shapefile.iterrows():
        tasks.append(
            dict(
                latitude=float(row[lat_col]),
                longitude=float(row[lon_col]),
                point_id=str(row[id_col]),
                output_dir=output_dir,
                start_date=start_date,
                end_date=end_date,
                log_file=log_file,
            )
        )

    with ProcessPoolExecutor(max_workers=n_cores) as pool:
        futures = {pool.submit(_process_single_point, t): t["point_id"] for t in tasks}
        for fut in as_completed(futures):
            pid = futures[fut]
            try:
                fut.result()
            except Exception as exc:
                print(f"ERROR (point {pid}): {exc}")

    print(f"\nNASA-POWER processing complete. Check the '{output_dir}' directory.\n")
