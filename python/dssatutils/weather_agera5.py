# File: weather_agera5.py
# ---------------------------------------------------------------------------
# Weather source: AgERA5 (ECMWF agrometeorological reanalysis) -> DSSAT .WTH.
#
# WHY: AgERA5 is ERA5 reprocessed *specifically for agriculture* — global,
# 0.1° (~10 km), daily, 1979–present, with the exact daily statistics crop
# models need (24 h max/min/mean temperature, daily solar radiation flux,
# precipitation flux, RH, wind, dewpoint). It covers the poles (unlike CHIRPS)
# and is higher-resolution than NASA POWER.
#
# ACCESS (requires a free key — NOT keyless like the other global sources):
#   1. Register at the Copernicus Climate Data Store: https://cds.climate.copernicus.eu/
#   2. Put your key in ~/.cdsapirc  (see https://cds.climate.copernicus.eu/how-to-api):
#        url: https://cds.climate.copernicus.eu/api
#        key: <your-personal-access-token>
#   3. pip install cdsapi xarray netcdf4
#   4. ONE-TIME: accept the dataset licence (otherwise requests 403 with
#      "required licences not accepted"):
#      https://cds.climate.copernicus.eu/datasets/sis-agrometeorological-indicators?tab=download#manage-licences
#   Dataset: "sis-agrometeorological-indicators"
#   Docs: https://cds.climate.copernicus.eu/datasets/sis-agrometeorological-indicators
#
# Requests are QUEUED by the CDS, so large/long runs can take a while; files are
# cached under `agera5_cache_dir` and reused. The network fetch is isolated in
# _download_agera5_var(); the .WTH formatting (_write_wth) is independently unit-
# testable with synthetic data (see tests/test_smoke.py pattern).
# ---------------------------------------------------------------------------

import os
import glob
import zipfile
import logging
from datetime import date, timedelta

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# AgERA5 CDS variable -> (cds 'variable', selector_kind, selector_value).
# selector_kind is "statistic" for most variables, but "time" for
# 2m_relative_humidity (provided at fixed hours, NOT as a 24-hour statistic),
# and None for fluxes that take no selector. Names follow the
# sis-agrometeorological-indicators catalogue.
_AGERA5_VARS = {
    "TMAX": ("2m_temperature", "statistic", "24_hour_maximum"),       # K
    "TMIN": ("2m_temperature", "statistic", "24_hour_minimum"),       # K
    "SRAD": ("solar_radiation_flux", None, None),                     # J/m²/day
    "RAIN": ("precipitation_flux", None, None),                       # mm/day
    "TDEW": ("2m_dewpoint_temperature", "statistic", "24_hour_mean"), # K
    "RH2M": ("2m_relative_humidity", "time", "15_00"),                # %  (mid-afternoon)
    "WIND": ("10m_wind_speed", "statistic", "24_hour_mean"),          # m/s
}
_CDS_DATASET = "sis-agrometeorological-indicators"


# ---------------------------------------------------------------------------
# Climatology helpers (shared convention with the other weather modules)
# ---------------------------------------------------------------------------

def _calc_tav(df: pd.DataFrame) -> float:
    return float(((df["TMAX"] + df["TMIN"]) / 2.0).mean())


def _calc_amp(df: pd.DataFrame) -> float:
    d = df.copy()
    d["TAVG"] = (d["TMAX"] + d["TMIN"]) / 2.0
    monthly = d.groupby(["YEAR", "MM"])["TAVG"].mean().reset_index()
    annual = monthly.groupby("YEAR")["TAVG"].agg(lambda x: x.max() - x.min())
    return float(annual.mean())


# ---------------------------------------------------------------------------
# Network fetch (ISOLATED so it can be mocked / debugged independently)
# ---------------------------------------------------------------------------

def _download_agera5_var(cds_var: str, sel_kind, sel_value, year: int,
                         area, cache_dir: str):
    """Download one AgERA5 variable-year over *area* via the CDS API.

    *sel_kind* is "statistic", "time", or None; *sel_value* the corresponding
    value. *area* = [north, west, south, east]. Returns the local zip path, or
    None on failure. Requires cdsapi + a configured ~/.cdsapirc key.
    """
    import cdsapi  # imported lazily so the module loads without the key/pkg

    tag = f"{cds_var}_{sel_value or 'na'}_{year}"
    dest = os.path.join(cache_dir, f"agera5_{tag}.zip")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest

    req = {
        "variable": cds_var,
        "year": str(year),
        "month": [f"{m:02d}" for m in range(1, 13)],
        "day": [f"{d:02d}" for d in range(1, 32)],
        "area": list(area),
        "version": "2_0",   # AgERA5 v2 (v1.1 deprecated from 2026-06-17)
    }
    if sel_kind is not None:
        req[sel_kind] = sel_value   # "statistic": "24_hour_mean"  or  "time": "15_00"

    try:
        cdsapi.Client().retrieve(_CDS_DATASET, req, dest)
        return dest if os.path.exists(dest) else None
    except Exception as exc:  # noqa: BLE001
        print(f"  AgERA5 download failed ({tag}): {exc}")
        return None


def _open_agera5(path: str):
    """Open an AgERA5 download (zip of daily ncs, or a single nc) as a dataset.

    AgERA5 delivers one netCDF per day; we eager-load each and concat along time.
    This avoids a hard dependency on dask (which xr.open_mfdataset requires) and
    is cheap because each file is a tiny per-grid-bbox subset.
    """
    import xarray as xr
    if path.endswith(".zip"):
        extract_dir = path[:-4] + "_nc"
        os.makedirs(extract_dir, exist_ok=True)
        with zipfile.ZipFile(path) as zf:
            zf.extractall(extract_dir)
        ncs = sorted(glob.glob(os.path.join(extract_dir, "*.nc")))
        if not ncs:
            return None
        dsets = []
        for f in ncs:
            with xr.open_dataset(f) as d:
                dsets.append(d.load())   # pull into memory, release the file handle
        if len(dsets) == 1:
            return dsets[0]
        combined = xr.concat(dsets, dim="time")
        return combined.sortby("time") if "time" in combined.coords else combined
    with xr.open_dataset(path) as d:
        return d.load()


# ---------------------------------------------------------------------------
# .WTH writer (TESTABLE with synthetic data; no network)
# ---------------------------------------------------------------------------

def _write_wth(df: pd.DataFrame, pid: str, lat: float, lon: float,
               output_dir: str) -> str:
    """Write one DSSAT .WTH from a daily DataFrame.

    *df* must contain DATE, SRAD, TMAX, TMIN, RAIN, TDEW, RH2M, WIND. Returns
    the output path. Shared formatting with the NASA POWER / Open-Meteo writers.
    """
    tav = _calc_tav(df)
    amp = _calc_amp(df)
    header = (
        f"$WEATHER DATA: AgERA5 (Point ID: {pid})\n"
        f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
        f"  AGE5 {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}   2.0  10.0\n"
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
    out_path = os.path.join(output_dir, f"{pid}.WTH")
    with open(out_path, "w") as fh:
        fh.write(header + "\n")
        fh.write("\n".join(lines) + "\n")
    return out_path


# ---------------------------------------------------------------------------
# Public entry point (mirrors the other weather sources + agera5_cache_dir)
# ---------------------------------------------------------------------------

def process_weather_agera5(
    shapefile,           # GeoDataFrame
    start_year: int,
    end_year: int,
    output_dir: str,
    id_col: str,
    lat_col: str,
    lon_col: str,
    n_cores: int,        # kept for API compatibility (extraction is serial I/O)
    log_file: str,
    agera5_cache_dir: str,
) -> None:
    """Download AgERA5 over the grid's bounding box and write DSSAT .WTH files.

    Requires a configured CDS API key (~/.cdsapirc) plus cdsapi + xarray. Data
    is downloaded once per variable-year (subset to the grid bbox) into
    *agera5_cache_dir* and reused. Unit conversions: temperature K→°C (−273.15),
    solar radiation J/m²/day → MJ/m²/day (×1e-6).
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(agera5_cache_dir, exist_ok=True)

    ids = [str(r[id_col]) for _, r in shapefile.iterrows()]
    lats = np.array([float(r[lat_col]) for _, r in shapefile.iterrows()])
    lons = np.array([float(r[lon_col]) for _, r in shapefile.iterrows()])

    # Bounding box (pad by ~0.2° so nearest-cell sampling always has coverage).
    pad = 0.2
    today = date.today()
    end_year = min(end_year, today.year)
    area = [float(lats.max() + pad), float(lons.min() - pad),
            float(lats.min() - pad), float(lons.max() + pad)]  # N, W, S, E

    print(f"--- Starting AgERA5 Download (Years: {start_year}–{end_year}) ---")
    print("  NOTE: AgERA5 requires a Copernicus CDS API key (~/.cdsapirc) and "
          "queues requests server-side; first run can be slow.")

    import xarray as xr  # noqa: F401  (fail early with a clear message if absent)
    from concurrent.futures import ThreadPoolExecutor

    # 1. Download every (variable, year) CONCURRENTLY. The CDS processes requests
    #    server-side, so submitting them in parallel overlaps the queue waits
    #    instead of paying them one after another (≈Nx faster). Cap concurrency
    #    at 4 to stay within the CDS per-user active-request limit.
    point_series = {pid: {v: {} for v in _AGERA5_VARS} for pid in ids}
    pts_lat = None

    jobs = [(dssat_var, cds_var, sel_kind, sel_value, year)
            for year in range(start_year, end_year + 1)
            for dssat_var, (cds_var, sel_kind, sel_value) in _AGERA5_VARS.items()]

    def _dl(job):
        dssat_var, cds_var, sel_kind, sel_value, year = job
        return (dssat_var, year,
                _download_agera5_var(cds_var, sel_kind, sel_value, year,
                                     area, agera5_cache_dir))

    paths = {}
    with ThreadPoolExecutor(max_workers=min(4, len(jobs))) as pool:
        for dssat_var, year, path in pool.map(_dl, jobs):
            if path:
                paths[(dssat_var, year)] = path

    # 2. Extract per-point daily series from each downloaded file (local + fast).
    for (dssat_var, year), path in paths.items():
        ds = _open_agera5(path)
        if ds is None:
            continue
        try:
            # Identify the data variable + coords.
            dv = list(ds.data_vars)[0]
            latname = "lat" if "lat" in ds.coords else "latitude"
            lonname = "lon" if "lon" in ds.coords else "longitude"
            if pts_lat is None:
                pts_lat = xr.DataArray(lats, dims="points")
                pts_lon = xr.DataArray(lons, dims="points")
            sel = ds[dv].sel({latname: pts_lat, lonname: pts_lon},
                             method="nearest")
            vals = np.asarray(sel.values)        # (time, points) or (points, time)
            times = pd.to_datetime(ds["time"].values)
            if vals.shape[0] != len(times):
                vals = vals.T
            # Unit conversions.
            if dssat_var in ("TMAX", "TMIN", "TDEW"):
                vals = vals - 273.15
            elif dssat_var == "SRAD":
                vals = vals * 1e-6
            date_codes = [f"{t.year}{t.dayofyear:03d}" for t in times]
            # Vectorised assignment: build a per-day dict in one pass per point.
            for j, pid in enumerate(ids):
                col = vals[:, j]
                good = ~np.isnan(col)
                point_series[pid][dssat_var].update(
                    {dc: float(v) for dc, v, g in zip(date_codes, col, good) if g})
        finally:
            ds.close()

    # 2. Assemble per-point frames and write .WTH.
    written = 0
    for pid, lat, lon in zip(ids, lats, lons):
        try:
            series = {v: pd.Series(point_series[pid][v], dtype="float64")
                      for v in _AGERA5_VARS}
            if series["TMAX"].empty:
                raise ValueError("No AgERA5 data extracted for this point.")
            frame = pd.DataFrame(series)
            frame.index.name = "DATE"
            frame = frame.reset_index()
            dts = pd.to_datetime(frame["DATE"], format="%Y%j")
            frame["YEAR"] = dts.dt.year
            frame["MM"] = dts.dt.month
            frame = frame.fillna(-99)
            _write_wth(frame, pid, lat, lon, output_dir)
            written += 1
        except Exception as exc:  # noqa: BLE001
            msg = f"\n--- ERROR ---\nAgERA5 point {pid} ({lat:.3f},{lon:.3f}): {exc}\n"
            print(msg)
            with open(log_file, "a") as lf:
                lf.write(msg)

    print(f"\nAgERA5 processing complete: {written}/{len(ids)} points written "
          f"to '{output_dir}'.\n")
