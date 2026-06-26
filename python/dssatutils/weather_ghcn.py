# Weather source: GHCN-Daily station observations (NOAA, live download).
#
# Global Historical Climatology Network - Daily is the standard ground-truth
# station archive. This backend snaps each grid point to the nearest station that
# actually has Tmax/Tmin (optionally precip) over the requested period and writes
# a DSSAT .WTH from the station record. SRAD/RH/wind are not GHCN-Daily core
# elements and are written DSSAT-missing (-99). Useful as an observational source
# and for bias-checking the gridded products.
#
# Access:
#   stations:  https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt
#   per-stn:   https://www.ncei.noaa.gov/access/services/data/v1  (server-side
#              filtered by station/period/element; returns metric units directly)

import io
import os

import numpy as np
import pandas as pd
import requests

from .weather_gridded_common import write_wth

_STATIONS_URL = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"
_DATA_SERVICE = "https://www.ncei.noaa.gov/access/services/data/v1"
# A descriptive User-Agent is good-citizen practice; NCEI also rate-limits by IP,
# so heavy gridded runs should reuse the cache and avoid re-querying in tight loops.
_HEADERS = {"User-Agent": "dssatutils/0.4 (DSSAT weather pipeline; research use)"}


def _load_stations(cache_dir: str) -> pd.DataFrame:
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, "ghcnd-stations.txt")
    if not os.path.exists(path):
        r = requests.get(_STATIONS_URL, headers=_HEADERS, timeout=300)
        r.raise_for_status()
        with open(path, "w") as fh:
            fh.write(r.text)
    ids, lats, lons = [], [], []
    with open(path) as fh:
        for line in fh:
            if len(line) < 31:
                continue
            try:
                ids.append(line[0:11].strip())
                lats.append(float(line[12:20]))
                lons.append(float(line[21:30]))
            except ValueError:
                continue
    return pd.DataFrame({"sid": ids, "lat": lats, "lon": lons})


def _haversine(lat0, lon0, lats, lons):
    R = 6371.0
    dlat = np.radians(lats - lat0)
    dlon = np.radians(lons - lon0)
    a = np.sin(dlat / 2) ** 2 + np.cos(np.radians(lat0)) * np.cos(np.radians(lats)) * np.sin(dlon / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def _fetch_station(sid: str, start_year: int, end_year: int) -> pd.DataFrame:
    """Return a daily frame (DATE, TMAX, TMIN, RAIN) for one station, or empty.

    Uses the NCEI data service with server-side station/period/element filtering
    and ``units=metric`` (Tmax/Tmin already in degC, precip in mm — no /10).
    """
    params = {
        "dataset": "daily-summaries", "stations": sid,
        "startDate": f"{start_year}-01-01", "endDate": f"{end_year}-12-31",
        "dataTypes": "TMAX,TMIN,PRCP", "units": "metric", "format": "csv",
    }
    try:
        r = requests.get(_DATA_SERVICE, params=params, headers=_HEADERS, timeout=120)
        if r.status_code != 200 or not r.text.strip():
            return pd.DataFrame()
        df = pd.read_csv(io.StringIO(r.text))
    except Exception:
        return pd.DataFrame()
    if "DATE" not in df or "TMAX" not in df or "TMIN" not in df:
        return pd.DataFrame()
    dts = pd.to_datetime(df["DATE"], errors="coerce")
    out = pd.DataFrame({
        "DATE": [f"{d.year}{d.dayofyear:03d}" for d in dts],
        "YEAR": dts.dt.year.values, "MM": dts.dt.month.values,
        "SRAD": -99.0,
        "TMAX": pd.to_numeric(df["TMAX"], errors="coerce").values,
        "TMIN": pd.to_numeric(df["TMIN"], errors="coerce").values,
        "RAIN": (pd.to_numeric(df["PRCP"], errors="coerce").values
                 if "PRCP" in df else -99.0),
        "TDEW": -99.0, "RH2M": -99.0, "WIND": -99.0,
    })
    return out[out["TMAX"].notna() & out["TMIN"].notna()]


def process_weather_ghcn(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    ghcn_cache_dir: str,
    max_candidates: int = 8,
) -> None:
    """Write DSSAT .WTH from the nearest GHCN-Daily station with valid Tmax/Tmin.

    For each grid point the nearest stations are tried in order until one returns
    Tmax/Tmin for the requested years. *ghcn_cache_dir* caches the station list.
    """
    os.makedirs(output_dir, exist_ok=True)
    stations = _load_stations(ghcn_cache_dir)
    pts = shapefile.copy()
    if hasattr(pts, "geometry"):
        pts = pts.to_crs("EPSG:4326")
        lats = pts.geometry.y.values.astype(float); lons = pts.geometry.x.values.astype(float)
    else:
        lats = pts[lat_col].astype(float).values; lons = pts[lon_col].astype(float).values
    ids = [str(r[id_col]) for _, r in pts.iterrows()]

    print(f"--- Starting GHCN-Daily Processing (Years: {start_year}-{end_year}) ---")
    written = 0
    for pid, lat, lon in zip(ids, lats, lons):
        try:
            order = np.argsort(_haversine(lat, lon, stations["lat"].values, stations["lon"].values))
            frame = pd.DataFrame()
            chosen = None
            for k in order[:max_candidates]:
                sid = stations["sid"].iloc[int(k)]
                frame = _fetch_station(sid, int(start_year), int(end_year))
                if not frame.empty:
                    chosen = sid
                    break
            if frame.empty:
                raise ValueError(f"no GHCN station with Tmax/Tmin near ({lat:.3f},{lon:.3f}) in {start_year}-{end_year}")
            write_wth(frame, pid, lat, lon, output_dir, f"GHCN-Daily {chosen}", "GHCN")
            written += 1
        except Exception as exc:  # noqa: BLE001
            msg = f"\n--- ERROR ---\nGHCN point {pid} ({lat:.3f},{lon:.3f}): {exc}\n"
            print(msg)
            with open(log_file, "a") as lf:
                lf.write(msg)
    print(f"\nGHCN-Daily processing complete: {written}/{len(ids)} point(s) written.\n")
