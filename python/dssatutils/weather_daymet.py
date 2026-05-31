# File: weather_daymet.py
# Python port of weather_daymet.R
#
# Fetches Daymet single-pixel data for each grid point and writes
# DSSAT-formatted .WTH files.  Daymet covers North America only and
# does NOT provide data for the current calendar year.
#
# API docs: https://daymet.ornl.gov/single-pixel/

import os
import io
import math
import logging
import time
from datetime import date, datetime, timedelta
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _calc_tav(df: pd.DataFrame) -> float:
    """Mean of daily mean temperature across the full period."""
    return float(df["t_avg"].mean())


def _calc_amp(df: pd.DataFrame) -> float:
    """
    Mean annual temperature amplitude:
      for each year → max(monthly mean) − min(monthly mean).
    Average those annual amplitudes.
    """
    monthly = df.groupby(["year", "month"])["t_avg"].mean().reset_index()
    annual = monthly.groupby("year")["t_avg"].agg(lambda x: x.max() - x.min())
    return float(annual.mean())


def _daymet_api_url(lat: float, lon: float, start_year: int, end_year: int) -> str:
    return (
        "https://daymet.ornl.gov/single-pixel/api/data"
        f"?lat={lat}&lon={lon}"
        "&vars=dayl,prcp,srad,tmax,tmin,vp"
        f"&start={start_year}-01-01&end={end_year}-12-31"
    )


def _download_daymet(lat: float, lon: float, start_year: int, end_year: int,
                     retries: int = 3, backoff: float = 5.0) -> pd.DataFrame:
    """
    Download Daymet data for a single point.
    Returns a DataFrame with columns:
      year, yday, dayl, prcp, srad, tmax, tmin, vp
    """
    url = _daymet_api_url(lat, lon, start_year, end_year)
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=120)
            r.raise_for_status()
            # Response is a CSV; skip the metadata header lines (non-numeric first lines)
            lines = r.text.splitlines()
            # Find the header row (contains 'year')
            header_idx = next(i for i, ln in enumerate(lines) if "year" in ln.lower())
            csv_text = "\n".join(lines[header_idx:])
            df = pd.read_csv(io.StringIO(csv_text))
            # Normalize column names: strip whitespace and units
            df.columns = [c.strip().split(" ")[0] for c in df.columns]
            return df
        except Exception as exc:
            if attempt < retries - 1:
                time.sleep(backoff * (attempt + 1))
            else:
                raise RuntimeError(f"Daymet download failed: {exc}") from exc


def _is_leap(year: int) -> bool:
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def _process_single_point(args: dict) -> None:
    """
    Worker function (run in a subprocess).
    args keys: latitude, longitude, point_id, output_dir, start_year, end_year, log_file
    """
    lat = args["latitude"]
    lon = args["longitude"]
    pid = args["point_id"]
    output_dir = args["output_dir"]
    start_year = args["start_year"]
    end_year = args["end_year"]
    log_file = args["log_file"]

    out_path = os.path.join(output_dir, f"{pid}.WTH")
    if os.path.exists(out_path):
        return

    try:
        raw = _download_daymet(lat, lon, start_year, end_year)

        # Rename to standard names
        col_map = {
            "dayl": "dayl",
            "prcp": "prcp",
            "srad": "srad",
            "tmax": "tmax",
            "tmin": "tmin",
            "vp":   "vp",
        }
        raw = raw.rename(columns={k: v for k, v in col_map.items() if k in raw.columns})

        # Mean temperature & month for TAV/AMP
        raw["t_avg"] = (raw["tmax"] + raw["tmin"]) / 2.0
        raw["month"] = pd.to_datetime(
            raw["year"].astype(str) + raw["yday"].astype(str).str.zfill(3),
            format="%Y%j"
        ).dt.month

        tav = _calc_tav(raw)
        amp = _calc_amp(raw)

        # Solar radiation: srad (W/m²) × dayl (s/day) / 1e6  → MJ/m²/day
        raw["srad_mj"] = (raw["srad"] * raw["dayl"]) / 1_000_000.0

        # Dewpoint from vapour pressure (Pa)
        vp_Pa = raw["vp"]
        raw["tdew"] = (237.3 * np.log(vp_Pa / 611.2)) / (17.27 - np.log(vp_Pa / 611.2))

        # Relative humidity
        es_Pa = 611.2 * np.exp((17.67 * raw["t_avg"]) / (raw["t_avg"] + 243.5))
        raw["rh2m"] = np.clip(100.0 * (vp_Pa / es_Pa), 0, 100)

        raw["wind"] = -99.0
        raw["DATE"] = raw["year"].astype(str) + raw["yday"].astype(str).str.zfill(3)

        # Daymet uses a 365-day year; duplicate DOY 365 as DOY 366 for leap years
        leap_years = [y for y in range(start_year, end_year + 1) if _is_leap(y)]
        extra_rows = []
        for yr in leap_years:
            mask = (raw["year"] == yr) & (raw["yday"] == 365)
            if mask.any():
                row = raw[mask].copy()
                row["yday"] = 366
                row["DATE"] = f"{yr}366"
                extra_rows.append(row)
        if extra_rows:
            raw = pd.concat([raw] + extra_rows, ignore_index=True)

        raw = raw.sort_values(["year", "yday"]).reset_index(drop=True)

        # --- Write WTH file ---
        header = (
            f"$WEATHER DATA: DayMet Data (Point ID: {pid})\n"
            f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
            f" DMET  {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}   -99   -99\n"
            f"@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND"
        )

        lines = []
        for _, row in raw.iterrows():
            line = (
                f"{row['DATE']:>7s}"
                f"{row['srad_mj']:6.1f}"
                f"{row['tmax']:6.1f}"
                f"{row['tmin']:6.1f}"
                f"{row['prcp']:6.1f}"
                f"{row['tdew']:6.1f}"
                f"{row['rh2m']:6.1f}"
                f"{row['wind']:6.1f}"
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

def process_weather_daymet(
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
    Download Daymet weather data for every point in *shapefile* and write
    DSSAT .WTH files to *output_dir*.  Mirrors the R ``process_weather_daymet``
    function, including the current-year guard and per-file resume logic.
    """
    today = date.today()
    current_year = today.year
    if end_year >= current_year:
        print(
            f"WARNING: Daymet data not available for the current year ({current_year}). "
            f"Adjusting end year to {current_year - 1}."
        )
        end_year = current_year - 1

    print(f"--- Starting DAYMET Download (Years: {start_year}–{end_year}) ---")
    print(f"Registered {n_cores} cores for parallel Daymet download.")

    os.makedirs(output_dir, exist_ok=True)

    tasks = []
    for _, row in shapefile.iterrows():
        tasks.append(
            dict(
                latitude=float(row[lat_col]),
                longitude=float(row[lon_col]),
                point_id=str(row[id_col]),
                output_dir=output_dir,
                start_year=start_year,
                end_year=end_year,
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

    print(f"\nDaymet processing complete. Check the '{output_dir}' directory.\n")
