# File: weather_era5land.py
# Python port of weather_era5land.R
#
# Downloads ERA5-Land point time series data from Copernicus Climate Data Store (CDS)
# in CSV format, aggregates it to daily statistics, and writes DSSAT-formatted .WTH files.
#
# CDS dataset: "reanalysis-era5-land-timeseries"

import os
import re
import gc
import datetime
from datetime import date, timedelta
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

from .credentials import make_cds_client

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _log_worker_message(log_file, level="INFO", point_id=None, msg=""):
    if not log_file:
        return
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    id_part = f" [ID={point_id}]" if point_id else ""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    clean_msg = " ".join(str(msg).split())
    line = f"[{timestamp}] [{level}] [WEATHER_ERA5_LAND]{id_part} {clean_msg}\n"
    with open(log_file, "a") as f:
        f.write(line)


def _fill_na_with_neighbor_mean(x, window=2, min_neighbors=2, max_iter=10):
    out = np.array(x, dtype=float)
    if not np.isnan(out).any():
        return out
    n = len(out)
    for _ in range(max_iter):
        missing_idx = np.where(np.isnan(out))[0]
        if len(missing_idx) == 0:
            break
        prev_out = out.copy()
        for idx in missing_idx:
            start_i = max(0, idx - window)
            end_i = min(n, idx + window + 1)
            neighbor_vals = [prev_out[i] for i in range(start_i, end_i) if i != idx]
            neighbor_vals = [v for v in neighbor_vals if not np.isnan(v)]
            if len(neighbor_vals) >= min_neighbors:
                out[idx] = sum(neighbor_vals) / len(neighbor_vals)
        if np.array_equal(np.isnan(out), np.isnan(prev_out)):
            break
    return out


def _find_first_matching_column(cols, patterns):
    for pattern in patterns:
        for col in cols:
            if re.search(pattern, col, re.IGNORECASE):
                return col
    return None


def _parse_hourly_datetime(df):
    cols = df.columns
    date_col = _find_first_matching_column(cols, [r"^date$", r"date", r"valid_date"])
    time_col = _find_first_matching_column(cols, [r"^time$", r"hour", r"valid_time"])
    dt_col = _find_first_matching_column(cols, [r"datetime", r"date_time", r"valid_datetime", r"^timestamp$", r"^time$"])

    if dt_col is not None:
        try:
            return pd.to_datetime(df[dt_col], utc=True)
        except Exception:
            pass

    if date_col is not None and time_col is not None:
        try:
            return pd.to_datetime(df[date_col].astype(str) + " " + df[time_col].astype(str), utc=True)
        except Exception:
            pass

    if date_col is not None:
        try:
            return pd.to_datetime(df[date_col], utc=True)
        except Exception:
            pass

    raise ValueError("Could not identify a datetime column in the ERA5-Land CSV download.")


def _calc_rh_from_temp_dew(temp_c, dew_c):
    rh = 100 * np.exp((17.625 * dew_c) / (243.04 + dew_c) - (17.625 * temp_c) / (243.04 + temp_c))
    return np.clip(rh, 0.0, 100.0)


def _cap_end_date_era5_land(end_year, lag_days=5):
    requested_end = date(end_year, 12, 31)
    safe_end = date.today() - timedelta(days=lag_days)
    return min(requested_end, safe_end)


def _download_era5_land_point_csv(latitude, longitude, start_date_str, end_date_str, target_file, cds_user="ecmwfr"):
    try:
        import cdsapi  # imported lazily
    except ImportError as exc:
        raise ImportError("ERA5-Land requires the 'cdsapi' Python package.") from exc
    c = make_cds_client(cdsapi)

    req = {
        'variable': [
            '2m_temperature',
            '2m_dewpoint_temperature',
            'total_precipitation',
            'surface_solar_radiation_downwards',
            '10m_u_component_of_wind',
            '10m_v_component_of_wind'
        ],
        'location': {
            'latitude': float(latitude),
            'longitude': float(longitude)
        },
        'date': {
            'start': start_date_str,
            'end': end_date_str
        },
        'data_format': 'csv'
    }

    os.makedirs(os.path.dirname(target_file), exist_ok=True)
    c.retrieve('reanalysis-era5-land-timeseries', req, target_file)


def _read_era5_land_hourly_csv(csv_file):
    df = pd.read_csv(csv_file)
    if df.empty:
        raise ValueError(f"Downloaded ERA5-Land file is empty: {csv_file}")

    cols = df.columns
    time_values = _parse_hourly_datetime(df)

    temp_col = _find_first_matching_column(cols, [r"^2m_temperature$", r"2m_temperature"])
    dew_col  = _find_first_matching_column(cols, [r"^2m_dewpoint_temperature$", r"2m_dewpoint_temperature"])
    tp_col   = _find_first_matching_column(cols, [r"^total_precipitation$", r"total_precipitation"])
    ssrd_col = _find_first_matching_column(cols, [r"^surface_solar_radiation_downwards$", r"surface_solar_radiation_downwards"])
    u10_col  = _find_first_matching_column(cols, [r"^10m_u_component_of_wind$", r"10m_u_component_of_wind"])
    v10_col  = _find_first_matching_column(cols, [r"^10m_v_component_of_wind$", r"10m_v_component_of_wind"])

    needed = [temp_col, dew_col, tp_col, ssrd_col, u10_col, v10_col]
    if any(c is None for c in needed):
        found = ", ".join(cols)
        raise ValueError(f"Could not find all required ERA5-Land columns in {csv_file}. Found: {found}")

    return pd.DataFrame({
        'DATETIME_UTC': time_values,
        'T2M_K': pd.to_numeric(df[temp_col], errors='coerce'),
        'DEW_K': pd.to_numeric(df[dew_col], errors='coerce'),
        'TP_M': pd.to_numeric(df[tp_col], errors='coerce'),
        'SSRD_J': pd.to_numeric(df[ssrd_col], errors='coerce'),
        'U10': pd.to_numeric(df[u10_col], errors='coerce'),
        'V10': pd.to_numeric(df[v10_col], errors='coerce'),
    })


def _aggregate_era5_land_to_daily(hourly_df, start_date_str, end_date_str, utc_offset_hours=None):
    if utc_offset_hours is not None:
        hourly_df["DATETIME_LOCAL"] = hourly_df["DATETIME_UTC"] + pd.to_timedelta(float(utc_offset_hours), unit="h")
    else:
        hourly_df["DATETIME_LOCAL"] = hourly_df["DATETIME_UTC"]

    hourly_df["DATE_obj"] = hourly_df["DATETIME_LOCAL"].dt.date
    hourly_df["T2M_C"] = hourly_df["T2M_K"] - 273.15
    hourly_df["DEW_C"] = hourly_df["DEW_K"] - 273.15
    hourly_df["RAIN_MM"] = np.maximum(hourly_df["TP_M"], 0.0) * 1000.0
    hourly_df["SRAD_MJ"] = np.maximum(hourly_df["SSRD_J"], 0.0) / 1e6
    hourly_df["WIND_MS"] = np.sqrt(hourly_df["U10"]**2 + hourly_df["V10"]**2)
    hourly_df["RH2M"] = _calc_rh_from_temp_dew(hourly_df["T2M_C"], hourly_df["DEW_C"])

    daily = hourly_df.groupby("DATE_obj").agg(
        SRAD=("SRAD_MJ", "sum"),
        TMAX=("T2M_C", "max"),
        TMIN=("T2M_C", "min"),
        RAIN=("RAIN_MM", "sum"),
        TDEW=("DEW_C", "mean"),
        RH2M=("RH2M", "mean"),
        WIND=("WIND_MS", "mean")
    ).reset_index()

    start_d = pd.to_datetime(start_date_str).date()
    end_d = pd.to_datetime(end_date_str).date()

    full_dates = pd.date_range(start_d, end_d, freq="D").date
    full_calendar = pd.DataFrame({"DATE_obj": full_dates})
    weather_data = pd.merge(full_calendar, daily, on="DATE_obj", how="left")

    vars_to_repair = ["SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND"]
    missing_before = {v: int(weather_data[v].isna().sum()) for v in vars_to_repair}

    for v in vars_to_repair:
        weather_data[v] = _fill_na_with_neighbor_mean(weather_data[v].values)

    weather_data["RAIN"] = np.maximum(weather_data["RAIN"], 0.0)
    weather_data["RH2M"] = np.clip(weather_data["RH2M"], 0.0, 100.0)

    bad_temp_idx = weather_data["TMAX"] < weather_data["TMIN"]
    if bad_temp_idx.any():
        tmp = weather_data.loc[bad_temp_idx, "TMAX"].copy()
        weather_data.loc[bad_temp_idx, "TMAX"] = weather_data.loc[bad_temp_idx, "TMIN"]
        weather_data.loc[bad_temp_idx, "TMIN"] = tmp

    missing_after = {v: int(np.isnan(weather_data[v].values).sum()) for v in vars_to_repair}

    weather_data["DATE_obj"] = pd.to_datetime(weather_data["DATE_obj"])
    weather_data["YEAR"] = weather_data["DATE_obj"].dt.year
    weather_data["MM"] = weather_data["DATE_obj"].dt.month
    weather_data["DOY"] = weather_data["DATE_obj"].dt.dayofyear

    return weather_data, missing_before, missing_after


def _write_dssat_weather_file(weather_data, latitude, longitude, output_file, point_id):
    weather_data["TAVG"] = (weather_data["TMAX"] + weather_data["TMIN"]) / 2.0
    tav = float(weather_data["TAVG"].mean())

    monthly_temps = weather_data.groupby(["YEAR", "MM"])["TAVG"].mean().reset_index()
    annual_amps = monthly_temps.groupby("YEAR")["TAVG"].agg(lambda x: x.max() - x.min()).reset_index()
    amp = float(annual_amps["TAVG"].mean())

    lines = [
        f"$WEATHER DATA: ERA5-LAND  (Point ID: {point_id})",
        "@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT",
        f"E5LD    {latitude:8.4f} {longitude:8.4f}   -99 {tav:5.1f} {amp:5.1f}   -99   -99",
        "@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND"
    ]

    for _, row in weather_data.iterrows():
        date_str = f"{int(row['YEAR'])}{int(row['DOY']):03d}"
        srad = -99.0 if pd.isna(row['SRAD']) else float(row['SRAD'])
        tmax = -99.0 if pd.isna(row['TMAX']) else float(row['TMAX'])
        tmin = -99.0 if pd.isna(row['TMIN']) else float(row['TMIN'])
        rain = -99.0 if pd.isna(row['RAIN']) else float(row['RAIN'])
        tdew = -99.0 if pd.isna(row['TDEW']) else float(row['TDEW'])
        rh2m = -99.0 if pd.isna(row['RH2M']) else float(row['RH2M'])
        wind = -99.0 if pd.isna(row['WIND']) else float(row['WIND'])

        line = (
            f"{date_str:>7s}"
            f"{srad:6.1f}"
            f"{tmax:6.1f}"
            f"{tmin:6.1f}"
            f"{rain:6.1f}"
            f"{tdew:6.1f}"
            f"{rh2m:6.1f}"
            f"{wind:6.1f}"
        )
        lines.append(line)

    with open(output_file, "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def process_weather_era5_land(
    shapefile,
    start_year: int,
    end_year: int,
    output_dir: str,
    id_col: str,
    lat_col: str,
    lon_col: str,
    n_cores: int,
    log_file: str,
    cds_user: str = "ecmwfr",
    utc_offset_hours = None,
    cache_dir: str = None,
    keep_raw_downloads: bool = False,
    availability_lag_days: int = 5
) -> None:
    """Download ERA5-Land point time series data and write DSSAT .WTH files.

    Requires cdsapi installed and configured with ~/.cdsapirc API key.
    """
    print(f"--- Starting ERA5-Land Download (Years: {start_year}–{end_year}) ---")

    start_date_str = f"{start_year}-01-01"
    end_date = _cap_end_date_era5_land(end_year, lag_days=availability_lag_days)
    end_date_str = str(end_date)

    os.makedirs(output_dir, exist_ok=True)
    if cache_dir is None:
        cache_dir = os.path.join(output_dir, "_era5_cache")
    os.makedirs(cache_dir, exist_ok=True)

    jobs = []
    for _, row in shapefile.iterrows():
        point_id = str(row[id_col])
        lat = float(row[lat_col])
        lon = float(row[lon_col])
        out_file = os.path.join(output_dir, f"{point_id}.WTH")
        raw_csv = os.path.join(cache_dir, f"{point_id}_era5land_hourly.csv")
        if not os.path.exists(out_file):
            jobs.append((point_id, lat, lon, out_file, raw_csv))

    if not jobs:
        print("All ERA5-Land weather files already exist. Processing skipped.")
        return

    def _process_one(job):
        point_id, lat, lon, out_file, raw_csv = job
        try:
            _download_era5_land_point_csv(lat, lon, start_date_str, end_date_str, raw_csv, cds_user)
            hourly = _read_era5_land_hourly_csv(raw_csv)
            daily, missing_before, missing_after = _aggregate_era5_land_to_daily(
                hourly, start_date_str, date_str := end_date_str, utc_offset_hours
            )

            if any(val > 0 for val in missing_before.values()):
                repair_msg = "; ".join(f"{v} {missing_before[v]}->{missing_after[v]}" for v in missing_before)
                _log_worker_message(log_file, "WARN", point_id,
                                    f"Missing ERA5-Land daily values detected and repaired: {repair_msg}")

            _write_dssat_weather_file(daily, lat, lon, out_file, point_id)

            if not keep_raw_downloads and os.path.exists(raw_csv):
                os.remove(raw_csv)
            _log_worker_message(log_file, "INFO", point_id, f"Successfully created ERA5-Land weather file: {os.path.basename(out_file)}")
        except Exception as e:
            _log_worker_message(log_file, "ERROR", point_id, str(e))
            print(f"Point {point_id} failed: {e}")

    requested_cores = max(1, int(n_cores))
    if requested_cores == 1 or len(jobs) == 1:
        for job in jobs:
            _process_one(job)
    else:
        with ThreadPoolExecutor(max_workers=min(requested_cores, len(jobs))) as executor:
            list(executor.map(_process_one, jobs))

    print(f"ERA5-Land weather processing complete. Output: {output_dir}")
