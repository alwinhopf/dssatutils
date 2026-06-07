# File: weather_nasapower_chirps.py
# ---------------------------------------------------------------------------
# HYBRID weather source: NASA POWER (all variables) + CHIRPS (rainfall).
#
# WHY: NASA POWER is global and provides every variable DSSAT needs, but its
# precipitation is coarse (~0.5°, ~50 km). Rainfall is the single most
# spatially variable and yield-critical input for rainfed crops. CHIRPS
# (Climate Hazards Group InfraRed Precipitation with Station data) is daily,
# ~0.05° (~5.5 km), 1981–present, and station-blended — markedly better for
# precipitation, especially in the tropics / semi-arid regions (Africa, India).
#
# This module fetches NASA POWER per point (TMAX/TMIN/SRAD/RH/dewpoint/wind),
# then REPLACES the RAIN column with CHIRPS values sampled from the gridded
# CHIRPS daily netCDF. Where CHIRPS has no coverage (|lat| > 50°) or a no-data
# cell, it falls back to the NASA POWER rainfall, so the result is still global.
#
# CHIRPS data:  https://www.chc.ucsb.edu/data/chirps
# CHIRPS netCDF: https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/
# License: CHIRPS is public domain (U.S. Government work); cite Funk et al. 2015,
#          Sci. Data 2:150066. NASA POWER is freely available; cite per their docs.
#
# Mirrors process_weather_nasapower() but takes an extra `chirps_cache_dir`
# kwarg (like GridMET's cache dir), so the pipeline wires it the same way.
# ---------------------------------------------------------------------------

import os
import time
import logging
from datetime import date, timedelta
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

# Reuse NASA POWER fetch + climatology helpers (single source of truth).
from .weather_nasapower import _fetch_nasa_power, _calc_tav, _calc_amp

logger = logging.getLogger(__name__)

# CHIRPS spatial resolution: "p05" (~0.05°, ~5.5 km — recommended) or
# "p25" (~0.25°, ~28 km — lighter download, coarser). Overridable from the
# pipeline before calling the entry point.
CHIRPS_RESOLUTION: str = "p05"

# CHIRPS covers 50°S–50°N only; outside this band fall back to NASA POWER rain.
_CHIRPS_LAT_LIMIT = 50.0
_CHIRPS_NODATA = -9999.0
_CHIRPS_NC_URL = (
    "https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/"
    "{res}/chirps-v2.0.{year}.days_{res}.nc"
)


# ---------------------------------------------------------------------------
# CHIRPS download + extraction
# ---------------------------------------------------------------------------

def _download_chirps_year(year: int, res: str, cache_dir: str,
                          timeout: int = 7200) -> str | None:
    """Download one CHIRPS yearly netCDF to *cache_dir* (skip if present).

    Returns the local path, or None on failure.
    """
    import requests
    fname = f"chirps-v2.0.{year}.days_{res}.nc"
    dest = os.path.join(cache_dir, fname)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    url = _CHIRPS_NC_URL.format(res=res, year=year)
    tmp = dest + ".part"
    try:
        print(f"  Downloading CHIRPS {year} ({res})... (large file; cached after first run)")
        with requests.get(url, stream=True, timeout=timeout) as r:
            r.raise_for_status()
            with open(tmp, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
        os.replace(tmp, dest)
        return dest
    except Exception as exc:  # noqa: BLE001
        print(f"  CHIRPS {year} download failed: {exc}")
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass
        return None


def _extract_chirps_rain(nc_paths: list, ids: list, lats: np.ndarray,
                         lons: np.ndarray) -> dict:
    """Vectorised per-point CHIRPS rainfall extraction from yearly netCDFs.

    Returns {point_id: pandas.Series indexed by DATE-code 'YYYYDOY' → rain mm}.
    Points outside CHIRPS coverage / over no-data cells get NaN for those days
    (the caller falls back to NASA POWER there).
    """
    import xarray as xr

    # Build per-point empty series; fill from each yearly file.
    out = {pid: {} for pid in ids}
    pts_lat = xr.DataArray(lats, dims="points")
    pts_lon = xr.DataArray(lons, dims="points")

    for path in nc_paths:
        with xr.open_dataset(path) as ds:
            # CHIRPS uses 'precip' with coords latitude/longitude/time.
            var = "precip" if "precip" in ds else list(ds.data_vars)[0]
            latname = "latitude" if "latitude" in ds.coords else "lat"
            lonname = "longitude" if "longitude" in ds.coords else "lon"

            sel = ds[var].sel(
                {latname: pts_lat, lonname: pts_lon}, method="nearest"
            )  # dims: (time, points)
            vals = sel.values  # (ntime, npoints)
            times = pd.to_datetime(ds["time"].values)
            date_codes = [f"{t.year}{t.dayofyear:03d}" for t in times]

            vals = np.where(vals <= _CHIRPS_NODATA, np.nan, vals)
            date_codes = np.asarray(date_codes)
            # Vectorised: keep only non-NaN days per point (one pass per point,
            # no inner per-day Python loop).
            for j, pid in enumerate(ids):
                col = vals[:, j]
                good = ~np.isnan(col)
                if good.any():
                    out[pid].update(zip(date_codes[good].tolist(),
                                        col[good].astype(float).tolist()))

    return {pid: pd.Series(d, dtype="float64") for pid, d in out.items()}


# ---------------------------------------------------------------------------
# Per-point worker (NASA POWER fetch + CHIRPS rain merge + .WTH write)
# ---------------------------------------------------------------------------

def _process_single_point(args: dict) -> None:
    lat = args["latitude"]
    lon = args["longitude"]
    pid = args["point_id"]
    output_dir = args["output_dir"]
    start_date = args["start_date"]
    end_date = args["end_date"]
    log_file = args["log_file"]
    chirps_rain = args["chirps_rain"]   # dict DATE→mm (may be empty)
    res = args.get("chirps_resolution", CHIRPS_RESOLUTION)

    out_path = os.path.join(output_dir, f"{pid}.WTH")
    if os.path.exists(out_path):
        return

    try:
        df = _fetch_nasa_power(lat, lon, start_date, end_date)
        if df.empty:
            raise ValueError("No data returned from NASA POWER.")

        df = df.rename(columns={
            "ALLSKY_SFC_SW_DWN": "SRAD", "T2M_MAX": "TMAX", "T2M_MIN": "TMIN",
            "PRECTOTCORR": "RAIN", "T2MDEW": "TDEW", "RH2M": "RH2M",
            "WS2M": "WIND",
        })

        # --- Merge CHIRPS rainfall over NASA POWER rain ---
        n_chirps = 0
        if abs(lat) <= _CHIRPS_LAT_LIMIT and chirps_rain:
            cseries = pd.Series(chirps_rain, dtype="float64")
            mapped = df["DATE"].map(cseries)
            use = mapped.notna()
            df.loc[use, "RAIN"] = mapped[use].values
            n_chirps = int(use.sum())
        rain_source = (f"CHIRPS({res}) where available, "
                       f"{n_chirps} days; NASA-POWER otherwise"
                       if n_chirps else "NASA-POWER (CHIRPS unavailable here)")

        tav = _calc_tav(df)
        amp = _calc_amp(df)

        header = (
            f"$WEATHER DATA: NASA-POWER + CHIRPS rain (Point ID: {pid}) "
            f"[{rain_source}]\n"
            f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
            f"  NAPC {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}   2.0   2.0\n"
            f"@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND"
        )

        lines = []
        for _, row in df.iterrows():
            line = (
                f"{row['DATE']:>7s}"
                f"{row['SRAD']:6.1f}{row['TMAX']:6.1f}{row['TMIN']:6.1f}"
                f"{row['RAIN']:6.1f}{row['TDEW']:6.1f}{row['RH2M']:6.1f}"
                f"{row['WIND']:6.1f}"
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


# ---------------------------------------------------------------------------
# Public entry point (mirrors process_weather_nasapower + chirps_cache_dir)
# ---------------------------------------------------------------------------

def process_weather_nasapower_chirps(
    shapefile,           # GeoDataFrame
    start_year: int,
    end_year: int,
    output_dir: str,
    id_col: str,
    lat_col: str,
    lon_col: str,
    n_cores: int,
    log_file: str,
    chirps_cache_dir: str,
) -> None:
    """Download NASA POWER + CHIRPS and write hybrid DSSAT .WTH files.

    NASA POWER supplies all variables; CHIRPS replaces rainfall within its
    50°S–50°N coverage (NASA POWER rain is kept elsewhere). Requires xarray +
    netCDF4. The CHIRPS yearly netCDF is downloaded once into *chirps_cache_dir*
    and reused across points/runs.
    """
    today = date.today()
    current_year = today.year
    start_date = f"{start_year}0101"
    if end_year == current_year:
        safe_end = today - timedelta(days=2)
        end_date = safe_end.strftime("%Y%m%d")
        print(f"End year is current year. Fetching up to: {safe_end.isoformat()}")
    else:
        end_date = f"{end_year}1231"

    res = CHIRPS_RESOLUTION
    print(f"--- Starting NASA-POWER + CHIRPS({res}) Hybrid (Years: {start_year}–{end_year}) ---")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(chirps_cache_dir, exist_ok=True)

    ids = [str(r[id_col]) for _, r in shapefile.iterrows()]
    if hasattr(shapefile, "geometry") and shapefile.geometry.notna().any():
        try:
            g = shapefile.to_crs("EPSG:4326")
            lons = g.geometry.x.values.astype(float)
            lats = g.geometry.y.values.astype(float)
        except Exception:
            lats = np.array([float(r[lat_col]) for _, r in shapefile.iterrows()])
            lons = np.array([float(r[lon_col]) for _, r in shapefile.iterrows()])
    else:
        lats = np.array([float(r[lat_col]) for _, r in shapefile.iterrows()])
        lons = np.array([float(r[lon_col]) for _, r in shapefile.iterrows()])

    # --- 1. CHIRPS: download yearly netCDFs and extract per-point rain ---
    chirps_by_point: dict = {pid: {} for pid in ids}
    in_band = np.abs(lats) <= _CHIRPS_LAT_LIMIT
    if in_band.any():
        try:
            nc_paths = []
            for yr in range(start_year, end_year + 1):
                p = _download_chirps_year(yr, res, chirps_cache_dir)
                if p:
                    nc_paths.append(p)
            if nc_paths:
                print(f"  Extracting CHIRPS rainfall for {len(ids)} point(s) "
                      f"from {len(nc_paths)} year file(s)...")
                extracted = _extract_chirps_rain(nc_paths, ids, lats, lons)
                chirps_by_point = {pid: s.to_dict() for pid, s in extracted.items()}
            else:
                print("  No CHIRPS files available; falling back to NASA-POWER rain.")
        except Exception as exc:  # noqa: BLE001
            print(f"  CHIRPS extraction failed ({exc}); using NASA-POWER rain only.")
    else:
        print("  All points are outside CHIRPS coverage (|lat| > 50°); "
              "using NASA-POWER rain.")

    # --- 2. NASA POWER per point (parallel) + merge CHIRPS rain ---
    print(f"Registered {n_cores} cores for parallel NASA-POWER download.")
    tasks = []
    for _, row in shapefile.iterrows():
        pid = str(row[id_col])
        tasks.append(dict(
            latitude=float(row[lat_col]), longitude=float(row[lon_col]),
            point_id=pid, output_dir=output_dir,
            start_date=start_date, end_date=end_date, log_file=log_file,
            chirps_rain=chirps_by_point.get(pid, {}),
            chirps_resolution=res,
        ))

    with ProcessPoolExecutor(max_workers=n_cores) as pool:
        futures = {pool.submit(_process_single_point, t): t["point_id"] for t in tasks}
        for fut in as_completed(futures):
            pid = futures[fut]
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                print(f"ERROR (point {pid}): {exc}")

    print(f"\nNASA-POWER + CHIRPS processing complete. Check '{output_dir}'.\n")
