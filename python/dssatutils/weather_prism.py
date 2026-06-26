# Weather source: PRISM daily 4 km grids for the contiguous United States.

import os
import time
import zipfile
from datetime import date, timedelta

import numpy as np
import pandas as pd
import requests

from .weather_gridded_common import write_wth

_PRISM_URL = "https://services.nacse.org/prism/data/get/us/4km/{var}/{yyyymmdd}"
_VARS = {"ppt": "RAIN", "tmax": "TMAX", "tmin": "TMIN", "tdmean": "TDEW"}
# Polite spacing between NACSE requests (seconds) to avoid throttle responses.
_PRISM_REQUEST_DELAY = 1.0


def _download_grid(var: str, day: pd.Timestamp, cache_dir: str) -> str | None:
    ymd = day.strftime("%Y%m%d")
    out_dir = os.path.join(cache_dir, var, ymd)
    os.makedirs(out_dir, exist_ok=True)
    existing = [os.path.join(out_dir, f) for f in os.listdir(out_dir)
                if f.lower().endswith((".bil", ".tif", ".tiff"))]
    if existing:
        return existing[0]
    url = _PRISM_URL.format(var=var, yyyymmdd=ymd)
    zpath = os.path.join(out_dir, f"{var}_{ymd}.zip")
    try:
        # The NACSE PRISM service throttles rapid requests; a small polite delay
        # avoids being served a non-zip throttle page in place of the data.
        time.sleep(_PRISM_REQUEST_DELAY)
        with requests.get(url, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(zpath, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
        if not zipfile.is_zipfile(zpath):
            print(f"  PRISM {var} {ymd}: response was not a valid zip "
                  "(likely throttled); skipping this day.")
            return None
        with zipfile.ZipFile(zpath) as zf:
            zf.extractall(out_dir)
        existing = [os.path.join(out_dir, f) for f in os.listdir(out_dir)
                    if f.lower().endswith((".bil", ".tif", ".tiff"))]
        return existing[0] if existing else None
    except Exception as exc:  # noqa: BLE001
        print(f"  PRISM download failed for {var} {ymd}: {exc}")
        return None


def _sample_raster(path: str, lats, lons):
    import rasterio
    from pyproj import Transformer

    out = np.full(len(lats), np.nan, dtype=float)
    with rasterio.open(path) as src:
        dst = src.crs.to_string() if src.crs else "EPSG:4326"
        xs, ys = Transformer.from_crs("EPSG:4326", dst, always_xy=True).transform(lons, lats)
        nodata = src.nodata
        for i, cell in enumerate(src.sample(zip(xs, ys), masked=True)):
            v = cell[0]
            if v is np.ma.masked or np.ma.is_masked(v):
                continue
            if nodata is not None and float(v) == float(nodata):
                continue
            out[i] = float(v)
    return out


def process_weather_prism(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    prism_cache_dir: str,
) -> None:
    """Download/cache PRISM daily grids and write DSSAT .WTH files.

    PRISM provides precipitation, Tmax, Tmin, and mean dewpoint. Solar radiation,
    RH, and wind are not daily PRISM variables and are written as DSSAT missing
    values (-99).
    """
    latest_safe = pd.Timestamp(date.today() - timedelta(days=2))
    dates = pd.date_range(f"{start_year}-01-01", min(pd.Timestamp(f"{end_year}-12-31"), latest_safe), freq="D")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(prism_cache_dir, exist_ok=True)

    pts = shapefile.copy()
    if hasattr(pts, "geometry"):
        pts = pts.to_crs("EPSG:4326")
        pts[lat_col] = pts.geometry.y
        pts[lon_col] = pts.geometry.x
    ids = [str(r[id_col]) for _, r in pts.iterrows()]
    lats = np.array([float(r[lat_col]) for _, r in pts.iterrows()])
    lons = np.array([float(r[lon_col]) for _, r in pts.iterrows()])
    frames = {pid: [] for pid in ids}

    print(f"--- Starting PRISM Processing (Years: {start_year}-{end_year}) ---")
    for day in dates:
        day_vals = {}
        for var, dssat_name in _VARS.items():
            path = _download_grid(var, day, prism_cache_dir)
            if path:
                day_vals[dssat_name] = _sample_raster(path, lats, lons)
        for i, pid in enumerate(ids):
            frames[pid].append({
                "DATE": f"{day.year}{day.dayofyear:03d}",
                "YEAR": day.year, "MM": day.month,
                "SRAD": -99.0,
                "TMAX": day_vals.get("TMAX", np.full(len(ids), np.nan))[i],
                "TMIN": day_vals.get("TMIN", np.full(len(ids), np.nan))[i],
                "RAIN": day_vals.get("RAIN", np.full(len(ids), np.nan))[i],
                "TDEW": day_vals.get("TDEW", np.full(len(ids), -99.0))[i],
                "RH2M": -99.0,
                "WIND": -99.0,
            })

    written = 0
    for pid, lat, lon in zip(ids, lats, lons):
        df = pd.DataFrame(frames[pid])
        df = df[df["TMAX"].notna() & df["TMIN"].notna()].fillna(-99)
        if df.empty:
            with open(log_file, "a") as lf:
                lf.write(f"PRISM point {pid}: no valid TMAX/TMIN data extracted\n")
            continue
        write_wth(df, pid, lat, lon, output_dir, "PRISM 4km", "PRSM", refht=2.0, wndht=-99.0)
        written += 1
    print(f"\nPRISM processing complete: {written}/{len(ids)} point(s) written.\n")
