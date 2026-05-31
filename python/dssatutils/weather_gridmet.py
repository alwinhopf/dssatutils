# File: weather_gridmet.py
# Python port of weather_gridmet.R  (Chunked-serial version)
#
# Downloads GridMET NetCDF files (~1 file/variable/year), pre-computes
# cell indices once, then extracts time-series for each grid point and
# writes DSSAT-formatted .WTH files.  No parallel cluster is used;
# bottleneck is I/O, not CPU.
#
# GridMET docs: https://www.climatologylab.org/gridmet.html
# Data URL base: http://www.northwestknowledge.net/metdata/data/

import gc
import math
import os
from datetime import date, timedelta

import numpy as np
import pandas as pd
import requests
import xarray as xr


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_GRIDMET_VARS = {
    "TMIN": "tmmn",  # Kelvin
    "TMAX": "tmmx",  # Kelvin
    "RAIN": "pr",    # mm/day
    "SRAD": "srad",  # W/m²
}

_GRIDMET_BASE_URL = "http://www.northwestknowledge.net/metdata/data/"


def _calc_tav(tmax_arr: np.ndarray, tmin_arr: np.ndarray,
              dates: pd.DatetimeIndex) -> float:
    """Mean of daily mean temperature across the full period."""
    return float(((tmax_arr + tmin_arr) / 2.0).mean())


def _calc_amp(tmax_arr: np.ndarray, tmin_arr: np.ndarray,
              dates: pd.DatetimeIndex) -> float:
    """Mean annual temperature amplitude (monthly means)."""
    df = pd.DataFrame({"tmax": tmax_arr, "tmin": tmin_arr}, index=dates)
    df["tavg"] = (df["tmax"] + df["tmin"]) / 2.0
    monthly = df["tavg"].resample("ME").mean()
    annual = monthly.resample("YE").agg(lambda x: x.max() - x.min())
    return float(annual.mean())


def _download_nc(url: str, dest: str, timeout: int = 3600) -> bool:
    """Stream-download a NetCDF file; return True on success."""
    try:
        with requests.get(url, stream=True, timeout=timeout) as r:
            r.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
        return True
    except Exception as exc:
        print(f"Download error ({url}): {exc}")
        return False


def _find_nearest_indices(ds: xr.Dataset, lons: np.ndarray,
                           lats: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Return arrays of integer (lat_idx, lon_idx) for each (lon, lat) pair.
    Returns NaN index where the point falls outside the dataset extent.
    """
    grid_lats = ds["lat"].values
    grid_lons = ds["lon"].values

    lat_min, lat_max = grid_lats.min(), grid_lats.max()
    lon_min, lon_max = grid_lons.min(), grid_lons.max()

    lat_idxs = np.full(len(lats), -1, dtype=int)
    lon_idxs = np.full(len(lons), -1, dtype=int)

    for i, (la, lo) in enumerate(zip(lats, lons)):
        if not (lat_min <= la <= lat_max and lon_min <= lo <= lon_max):
            lat_idxs[i] = -1
            lon_idxs[i] = -1
        else:
            lat_idxs[i] = int(np.argmin(np.abs(grid_lats - la)))
            lon_idxs[i] = int(np.argmin(np.abs(grid_lons - lo)))

    return lat_idxs, lon_idxs


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def process_weather_gridmet(
    shapefile,               # GeoDataFrame
    start_year: int,
    end_year: int,
    output_dir: str,
    id_col: str,
    lat_col: str,
    lon_col: str,
    n_cores: int,            # kept for API compatibility; serial implementation
    log_file: str,
    gridmet_cache_dir: str,
    chunk_size: int = 3000,
) -> None:
    """
    Download GridMET NetCDF files and write DSSAT .WTH files for every point
    in *shapefile*.  Mirrors the R ``process_weather_gridmet`` function exactly,
    including chunked processing, unit conversions, and file-exists resume.

    Unit conversions applied:
      • Temperature (K → °C): subtract 273.15
      • Solar radiation (W/m² → MJ/m²/day): multiply by 0.0864
    """
    print(f"--- Starting GridMET Serial Download (Years: {start_year}–{end_year}) ---")

    # -----------------------------------------------------------------------
    # 1. Date logic (clip to data availability)
    # -----------------------------------------------------------------------
    start_date = pd.Timestamp(f"{start_year}-01-01")
    latest_safe = pd.Timestamp(date.today() - timedelta(days=2))
    requested_end = pd.Timestamp(f"{end_year}-12-31")
    end_date = min(requested_end, latest_safe)
    if end_date < requested_end:
        print(f"NOTICE: Adjusting end date to {end_date.date()} (data availability limit).")

    full_date_seq = pd.date_range(start_date, end_date, freq="D")
    years_needed = sorted(full_date_seq.year.unique().tolist())

    # -----------------------------------------------------------------------
    # 2. Download NetCDF files
    # -----------------------------------------------------------------------
    os.makedirs(gridmet_cache_dir, exist_ok=True)
    downloaded: dict[str, dict[int, str]] = {v: {} for v in _GRIDMET_VARS}

    for var_name, abbrev in _GRIDMET_VARS.items():
        for yr in years_needed:
            fname = f"{abbrev}_{yr}.nc"
            dest = os.path.join(gridmet_cache_dir, fname)
            url = _GRIDMET_BASE_URL + fname
            if not os.path.exists(dest):
                print(f"  Downloading {fname}...")
                _download_nc(url, dest)
            if os.path.exists(dest):
                downloaded[var_name][yr] = dest

    # -----------------------------------------------------------------------
    # 3. Pre-compute cell indices using the first TMIN file as spatial reference
    # -----------------------------------------------------------------------
    print("\n--- Calculating Grid Indices ---")
    first_yr = years_needed[0]
    ref_path = downloaded.get("TMIN", {}).get(first_yr)
    if ref_path is None:
        raise FileNotFoundError("Cannot find base TMIN NetCDF for spatial reference.")

    with xr.open_dataset(ref_path) as ref_ds:
        point_lats = shapefile[lat_col].values.astype(float)
        point_lons = shapefile[lon_col].values.astype(float)
        lat_idxs, lon_idxs = _find_nearest_indices(ref_ds, point_lons, point_lats)

    invalid = lat_idxs == -1
    if invalid.any():
        print(f"WARNING: {invalid.sum()} point(s) fall outside the GridMET coverage area. "
              "They will be skipped.")

    # -----------------------------------------------------------------------
    # 4. Chunked processing
    # -----------------------------------------------------------------------
    total_points = len(shapefile)
    num_chunks = math.ceil(total_points / chunk_size)
    print(f"\n--- Processing {total_points} points in {num_chunks} chunks "
          f"(chunk size: {chunk_size}) ---")

    os.makedirs(output_dir, exist_ok=True)

    for k in range(num_chunks):
        start_idx = k * chunk_size
        end_idx = min((k + 1) * chunk_size, total_points)
        chunk_indices = list(range(start_idx, end_idx))
        print(f"\nProcessing Chunk {k+1}/{num_chunks} "
              f"(Points {start_idx+1} to {end_idx})...")

        chunk_lat_idxs = lat_idxs[chunk_indices]
        chunk_lon_idxs = lon_idxs[chunk_indices]

        if (chunk_lat_idxs == -1).all():
            print("  Skipping chunk (all points outside GridMET grid).")
            continue

        # Build per-variable time-series matrices: shape (n_chunk_pts, n_days)
        chunk_data: dict[str, np.ndarray] = {}

        for var_name in _GRIDMET_VARS:
            year_arrays = []
            sorted_years = sorted(downloaded[var_name].keys())

            for yr in sorted_years:
                fpath = downloaded[var_name][yr]
                with xr.open_dataset(fpath) as ds:
                    # Variable name in file is the abbrev (tmmn, tmmx, pr, srad)
                    abbrev = _GRIDMET_VARS[var_name]
                    # Identify the actual variable name in the dataset
                    var_key = abbrev if abbrev in ds else list(ds.data_vars)[0]
                    data_3d = ds[var_key].values  # (time, lat, lon)

                    # Days for this year within our date range
                    yr_dates = pd.date_range(f"{yr}-01-01", f"{yr}-12-31", freq="D")
                    target_dates = yr_dates[(yr_dates >= start_date) &
                                            (yr_dates <= end_date)]
                    n_expected = len(target_dates)

                    n_time = data_3d.shape[0]
                    n_use = min(n_time, n_expected)

                    # Extract rows for chunk points
                    pt_vals = np.full((len(chunk_indices), n_use), np.nan)
                    for ci, (li, loi) in enumerate(zip(chunk_lat_idxs, chunk_lon_idxs)):
                        if li == -1:
                            continue
                        pt_vals[ci, :] = data_3d[:n_use, li, loi]

                    year_arrays.append(pt_vals)
                    del data_3d
                gc.collect()

            chunk_data[var_name] = np.concatenate(year_arrays, axis=1)
            del year_arrays
            gc.collect()

        # Sync column count across variables
        col_counts = [v.shape[1] for v in chunk_data.values()]
        min_cols = min(col_counts)
        chunk_dates = full_date_seq[:min_cols]
        for vn in chunk_data:
            if chunk_data[vn].shape[1] > min_cols:
                chunk_data[vn] = chunk_data[vn][:, :min_cols]

        # Write WTH files (serial)
        print("  -> Writing .WTH files...")
        for ci, global_idx in enumerate(chunk_indices):
            row = shapefile.iloc[global_idx]
            pid = str(row[id_col])
            lat = float(row[lat_col])
            lon = float(row[lon_col])

            out_f = os.path.join(output_dir, f"{pid}.WTH")
            if os.path.exists(out_f) or chunk_lat_idxs[ci] == -1:
                continue

            try:
                tmin_arr = chunk_data["TMIN"][ci] - 273.15  # K → °C
                tmax_arr = chunk_data["TMAX"][ci] - 273.15
                rain_arr = chunk_data["RAIN"][ci]
                srad_arr = chunk_data["SRAD"][ci] * 0.0864  # W/m² → MJ/m²/day

                if np.isnan(tmin_arr).any():
                    continue

                tdew_arr = tmin_arr - 2.5
                rh2m_arr = np.clip(100.0 - (tmax_arr - tmin_arr) * 2.0, 20.0, 100.0)
                wind_arr = np.full_like(tmin_arr, -99.0)

                tav = float(((tmax_arr + tmin_arr) / 2.0).mean())
                amp = _calc_amp(tmax_arr, tmin_arr, chunk_dates)

                header = (
                    f"$WEATHER DATA: GridMET Data (Point ID: {pid})\n"
                    f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
                    f" GMET  {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}   -99   -99\n"
                    f"@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND"
                )

                lines = []
                for di, dt in enumerate(chunk_dates):
                    date_str = f"{dt.year}{dt.day_of_year:03d}"
                    line = (
                        f"{date_str:>7s}"
                        f"{srad_arr[di]:6.1f}"
                        f"{tmax_arr[di]:6.1f}"
                        f"{tmin_arr[di]:6.1f}"
                        f"{rain_arr[di]:6.1f}"
                        f"{tdew_arr[di]:6.1f}"
                        f"{rh2m_arr[di]:6.1f}"
                        f"{wind_arr[di]:6.1f}"
                    )
                    line = line.replace(" -99.0", "   -99")
                    lines.append(line)

                with open(out_f, "w") as fh:
                    fh.write(header + "\n")
                    fh.write("\n".join(lines) + "\n")

            except Exception as exc:
                msg = f"Error on point {pid}: {exc}"
                print(f"  {msg}")
                with open(log_file, "a") as lf:
                    lf.write(msg + "\n")

        del chunk_data
        gc.collect()

    print(f"GridMET processing complete. Output: {output_dir}")
