# Weather source: PRISM daily 4 km grids for the contiguous United States.

from __future__ import annotations

import json
import os
import time
import tempfile
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, timedelta
from typing import Optional

import numpy as np
import pandas as pd
import requests

from .weather_gridded_common import write_wth
from .weather_solar import estimate_srad_bristow_campbell
from .weather_validation import is_wth_valid

_PRISM_URL = "https://services.nacse.org/prism/data/get/us/4km/{var}/{yyyymmdd}"
_VARS = {"ppt": "RAIN", "tmax": "TMAX", "tmin": "TMIN", "tdmean": "TDEW"}
# Polite spacing between NACSE requests (seconds) to avoid throttle responses.
_PRISM_REQUEST_DELAY = 1.0

# ACIS Web Services GridData endpoint and PRISM grid code (grid 21)
_ACIS_GRIDDATA_URL = "https://data.rcc-acis.org/GridData"
_ACIS_PRISM_GRID = "21"


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


def _acis_request(lat, lon, sdate, edate):
    # Match the coordinates actually sent to ACIS in both languages.
    return {"version": 2, "loc": f"{lon:.4f},{lat:.4f}", "grid": _ACIS_PRISM_GRID,
            "sdate": sdate, "edate": edate, "elements": "maxt,mint,pcpn",
            "units": "degreeC,degreeC,mm"}


def _acis_cache_file(cache_dir, lat, lon, sdate, edate):
    return os.path.join(cache_dir, f"prism_acis_v2_21_{lon:.4f}_{lat:.4f}_{sdate}_{edate}.json")


def _complete_acis_records(raw, sdate, edate):
    """Require exact ordered coverage and observed core forcing before caching."""
    try:
        dates = pd.date_range(sdate, edate).strftime("%Y-%m-%d").tolist()
        if not dates or not isinstance(raw, list) or len(raw) != len(dates):
            return False
        for row, expected in zip(raw, dates):
            if not isinstance(row, (list, tuple)) or len(row) != 4 or row[0] != expected:
                return False
            tmax, tmin = float(row[1]), float(row[2])
            rain = 0.0 if str(row[3]).strip().lower() == "t" else float(row[3])
            if (not np.isfinite([tmax, tmin, rain]).all()
                    or not -90 <= tmin <= tmax <= 70 or not 0 <= rain <= 2000):
                return False
        return True
    except (ValueError, TypeError, OverflowError):
        return False


def _publish_acis_wth(df, pid, lat, lon, output_dir, source_label, sdate, edate, srad_method):
    required = ["TMAX", "TMIN", "RAIN"]
    if srad_method != "none":
        required.append("SRAD")
    # Never replace an existing usable file with a partial or invalid download.
    with tempfile.TemporaryDirectory(prefix=".prism-", dir=output_dir) as stage:
        path = write_wth(df, pid, lat, lon, stage, source_label, "PRSM", refht=2.0, wndht=-99.0)
        if not is_wth_valid(path, start_date=sdate, end_date=edate, required_columns=required):
            return False
        os.replace(path, os.path.join(output_dir, f"{pid}.WTH"))
    return True


def _fetch_acis_point(
    lat: float,
    lon: float,
    sdate: str,
    edate: str,
    cache_file: str | None = None,
    max_retries: int = 3,
) -> list | None:
    request = _acis_request(lat, lon, sdate, edate)
    if cache_file and os.path.exists(cache_file):
        try:
            with open(cache_file, "r", encoding="utf-8") as fh:
                data = json.load(fh)
                if (isinstance(data, dict) and data.get("request") == request
                        and _complete_acis_records(data.get("data"), sdate, edate)):
                    return data["data"]
        except Exception:
            pass

    payload = {
        "loc": f"{lon:.4f},{lat:.4f}",
        "sdate": sdate,
        "edate": edate,
        "grid": _ACIS_PRISM_GRID,
        "elems": [
            {"name": "maxt", "units": "degreeC"},
            {"name": "mint", "units": "degreeC"},
            {"name": "pcpn", "units": "mm"},
        ],
    }

    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.post(_ACIS_GRIDDATA_URL, json=payload, timeout=60)
            if resp.status_code == 200:
                txt = resp.text.strip()
                if txt:
                    data = resp.json()
                    if isinstance(data, dict) and _complete_acis_records(data.get("data"), sdate, edate):
                        if cache_file:
                            parent = os.path.dirname(os.path.abspath(cache_file))
                            os.makedirs(parent, exist_ok=True)
                            fd, stage = tempfile.mkstemp(prefix=".prism-", suffix=".json", dir=parent)
                            try:
                                with os.fdopen(fd, "w", encoding="utf-8") as fh:
                                    json.dump({"request": request, "data": data["data"]}, fh)
                                os.replace(stage, cache_file)
                            finally:
                                if os.path.exists(stage):
                                    os.unlink(stage)
                        return data["data"]
        except Exception:
            pass
        if attempt < max_retries:
            time.sleep(1.0 * attempt)
    return None


def _parse_acis_records(
    raw_data: list,
    lat: float,
    srad_method: str = "bristow_campbell",
) -> pd.DataFrame:
    if not raw_data:
        return pd.DataFrame()
    records = []
    for row in raw_data:
        if len(row) < 4:
            continue
        date_str, tmax_raw, tmin_raw, pcpn_raw = row[0], row[1], row[2], row[3]
        try:
            d = pd.to_datetime(date_str)
        except Exception:
            continue
        date_dssat = f"{d.year}{d.dayofyear:03d}"

        try:
            tmax = float(tmax_raw)
            if tmax == -999.0 or not np.isfinite(tmax):
                tmax = np.nan
        except (ValueError, TypeError):
            tmax = np.nan

        try:
            tmin = float(tmin_raw)
            if tmin == -999.0 or not np.isfinite(tmin):
                tmin = np.nan
        except (ValueError, TypeError):
            tmin = np.nan

        if np.isnan(tmax) or np.isnan(tmin):
            continue

        if isinstance(pcpn_raw, str) and pcpn_raw.strip().lower() == "t":
            pcpn = 0.0
        else:
            try:
                pcpn = float(pcpn_raw)
                if pcpn == -999.0 or not np.isfinite(pcpn):
                    pcpn = np.nan
            except (ValueError, TypeError):
                pcpn = np.nan

        records.append({
            "DATE": date_dssat,
            "YEAR": d.year,
            "MM": d.month,
            "SRAD": -99.0,
            "TMAX": tmax,
            "TMIN": tmin,
            "RAIN": -99.0 if np.isnan(pcpn) else pcpn,
            "TDEW": -99.0,
            "RH2M": -99.0,
            "WIND": -99.0,
        })
    df = pd.DataFrame(records)
    if df.empty:
        return df
    if srad_method == "bristow_campbell":
        df["SRAD"] = estimate_srad_bristow_campbell(df["DATE"], df["TMAX"], df["TMIN"], lat)
    else:
        df["SRAD"] = -99.0
    return df


def process_weather_prism(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    prism_cache_dir: Optional[str] = None,
    srad_method: str = "bristow_campbell",
    backend: str = "acis",
) -> None:
    """Download PRISM daily data and write DSSAT .WTH files.

    Supports two backends:
    - 'acis' (default): Rapid point-by-point time-series queries via NOAA RCC ACIS
      GridData API (grid 21). Eliminates national raster downloads and requires
      minimal disk cache.
    - 'nacse': Legacy automated download of whole-CONUS 4 km daily rasters from
      Oregon State University (NACSE) and local spatial sampling.
    """
    backend_norm = str(backend or "acis").strip().lower()
    if backend_norm not in ("acis", "point", "nacse", "grid"):
        raise ValueError(
            f"Unknown PRISM backend: {backend!r}. Expected 'acis' or 'nacse'."
        )
    if backend_norm in ("nacse", "grid") and not prism_cache_dir:
        raise ValueError("prism_cache_dir is required for the NACSE grid backend")
    if srad_method not in ("bristow_campbell", "none"):
        raise ValueError(
            f"Unknown srad_method: {srad_method!r}. Expected 'bristow_campbell' or 'none'."
        )

    os.makedirs(output_dir, exist_ok=True)
    if prism_cache_dir:
        os.makedirs(prism_cache_dir, exist_ok=True)

    pts = shapefile.copy()
    if hasattr(pts, "geometry"):
        pts = pts.to_crs("EPSG:4326")
        pts[lat_col] = pts.geometry.y
        pts[lon_col] = pts.geometry.x
    ids = [str(r[id_col]) for _, r in pts.iterrows()]
    lats = np.array([float(r[lat_col]) for _, r in pts.iterrows()])
    lons = np.array([float(r[lon_col]) for _, r in pts.iterrows()])

    source_label = (
        "PRISM 4km (SRAD estimated: Bristow-Campbell 1984)"
        if srad_method == "bristow_campbell"
        else "PRISM 4km"
    )

    if backend_norm in ("acis", "point"):
        print(f"--- Starting PRISM ACIS Point Processing (Years: {start_year}-{end_year}) ---")
        if srad_method == "bristow_campbell":
            print("  PRISM: Estimating solar radiation using Bristow-Campbell (1984)")
        else:
            print("  PRISM: Solar radiation estimation disabled (srad_method='none')")

        sdate = f"{int(start_year):04d}-01-01"
        latest_safe = pd.Timestamp(date.today() - timedelta(days=2))
        edate = min(pd.Timestamp(f"{int(end_year):04d}-12-31"), latest_safe).strftime("%Y-%m-%d")
        if sdate > edate:
            raise ValueError("PRISM start date is after the available end date")

        cache_dir = os.path.join(prism_cache_dir, "acis") if prism_cache_dir else None
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)

        def _worker(pid, lat, lon):
            cfile = _acis_cache_file(cache_dir, lat, lon, sdate, edate) if cache_dir else None
            raw = _fetch_acis_point(lat, lon, sdate, edate, cfile)
            if not _complete_acis_records(raw, sdate, edate):
                return pid, lat, lon, pd.DataFrame()
            return pid, lat, lon, _parse_acis_records(raw, lat, srad_method)

        written = 0
        workers = min(max(1, int(n_cores)), 4)
        if workers > 1 and len(ids) > 1:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(_worker, pid, lat, lon) for pid, lat, lon in zip(ids, lats, lons)]
                for fut in as_completed(futures):
                    pid, lat, lon, df = fut.result()
                    if df.empty:
                        with open(log_file, "a", encoding="utf-8") as lf:
                            lf.write(f"PRISM point {pid}: incomplete or invalid ACIS history for {sdate} through {edate}; not published\n")
                        continue
                    df = df.fillna(-99)
                    if _publish_acis_wth(df, pid, lat, lon, output_dir, source_label, sdate, edate, srad_method):
                        written += 1
                    else:
                        with open(log_file, "a", encoding="utf-8") as lf:
                            lf.write(f"PRISM point {pid}: invalid WTH forcing; not published\n")
        else:
            for pid, lat, lon in zip(ids, lats, lons):
                _, _, _, df = _worker(pid, lat, lon)
                if df.empty:
                    with open(log_file, "a", encoding="utf-8") as lf:
                        lf.write(f"PRISM point {pid}: incomplete or invalid ACIS history for {sdate} through {edate}; not published\n")
                    continue
                df = df.fillna(-99)
                if _publish_acis_wth(df, pid, lat, lon, output_dir, source_label, sdate, edate, srad_method):
                    written += 1
                else:
                    with open(log_file, "a", encoding="utf-8") as lf:
                        lf.write(f"PRISM point {pid}: invalid WTH forcing; not published\n")
        print(f"\nPRISM processing complete: {written}/{len(ids)} point(s) written.\n")
        return

    # Legacy NACSE CONUS grid download / raster extraction backend
    latest_safe = pd.Timestamp(date.today() - timedelta(days=2))
    dates = pd.date_range(f"{start_year}-01-01", min(pd.Timestamp(f"{end_year}-12-31"), latest_safe), freq="D")
    frames = {pid: [] for pid in ids}

    print(f"--- Starting PRISM Processing (Years: {start_year}-{end_year}) ---")
    if srad_method == "bristow_campbell":
        print("  PRISM: Estimating solar radiation using Bristow-Campbell (1984)")
    else:
        print("  PRISM: Solar radiation estimation disabled (srad_method='none')")

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
        df = df[df["TMAX"].notna() & df["TMIN"].notna()].copy()
        if df.empty:
            with open(log_file, "a", encoding="utf-8") as lf:
                lf.write(f"PRISM point {pid}: no valid TMAX/TMIN data extracted\n")
            continue
        if srad_method == "bristow_campbell":
            df["SRAD"] = estimate_srad_bristow_campbell(df["DATE"], df["TMAX"], df["TMIN"], lat)
        else:
            df["SRAD"] = -99.0
        df = df.fillna(-99)
        write_wth(df, pid, lat, lon, output_dir, source_label, "PRSM", refht=2.0, wndht=-99.0)
        written += 1
    print(f"\nPRISM processing complete: {written}/{len(ids)} point(s) written.\n")
