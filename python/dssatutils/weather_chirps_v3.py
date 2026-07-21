# File: weather_chirps_v3.py
# ---------------------------------------------------------------------------
# CHIRPS v3 daily rainfall extraction and hybrid weather helpers.
#
# CHIRPS v3 daily data are precipitation-only. This module keeps the rainfall
# layer reusable: extract_chirps_v3_rainfall() returns per-point daily rainfall
# keyed by DSSAT DATE (YYYYDOY), and process_weather_nasapower_chirps_v3()
# demonstrates the first full-weather hybrid by replacing NASA POWER RAIN.
#
# CHIRPS v3 daily products:
#   final/rnl: ERA5-ratio daily disaggregation, full 1981+ period.
#   final/sat: IMERG-ratio daily disaggregation, recent period only.
# ---------------------------------------------------------------------------

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from datetime import date, timedelta
from typing import Literal
from urllib.parse import urlparse

import numpy as np
import pandas as pd

from .weather_nasapower import _fetch_nasa_power, _calc_tav, _calc_amp
from .weather_rainfall_merge import merge_rainfall_into_weather
from .config import get_config_number, get_config_value

CHIRPS_V3_PRODUCT: str = get_config_value("weather.chirps_v3.product", "rnl")
CHIRPS_V3_STREAM: str = get_config_value("weather.chirps_v3.stream", "final")
CHIRPS_V3_FETCH_MODE: str = get_config_value(
    "weather.chirps_v3.fetch_mode", "monthly_netcdf"
)
CHIRPS_V3_RESOLUTION: str = get_config_value("weather.chirps_v3.resolution", "p05")

_CHIRPS_V3_BASE_URL = get_config_value(
    "weather.chirps_v3.base_url",
    "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/daily",
)
_CHIRPS_V3_LAT_LIMIT = get_config_number("weather.chirps_v3.latitude_limit", 60.0)
_CHIRPS_V3_NODATA = get_config_number("weather.chirps_v3.nodata", -9999.0)


def _normalize_chirps_v3_options(
    product: str | None = None,
    stream: str | None = None,
    fetch_mode: str | None = None,
    resolution: str | None = None,
) -> tuple[str, str, str, str]:
    product = (product or CHIRPS_V3_PRODUCT).lower()
    stream = (stream or CHIRPS_V3_STREAM).lower()
    fetch_mode = (fetch_mode or CHIRPS_V3_FETCH_MODE).lower()
    resolution = (resolution or CHIRPS_V3_RESOLUTION).lower()
    if product not in {"rnl", "sat"}:
        raise ValueError("chirps_product must be 'rnl' or 'sat'")
    if stream not in {"final", "prelim"}:
        raise ValueError("chirps_stream must be 'final' or 'prelim'")
    if fetch_mode not in {"monthly_netcdf", "yearly_netcdf", "gee", "remote_cog"}:
        raise ValueError("chirps_fetch_mode must be 'monthly_netcdf', 'yearly_netcdf', 'gee', or 'remote_cog'")
    if resolution != "p05":
        raise ValueError("CHIRPS v3 daily data is currently supported only at resolution 'p05'")
    if stream == "prelim" and product != "sat":
        raise ValueError("CHIRPS v3 preliminary daily data is available only for product 'sat'")
    if stream == "prelim" and fetch_mode not in {"yearly_netcdf"}:
        raise ValueError("CHIRPS v3 preliminary daily NetCDF is currently exposed only byYear")
    return product, stream, fetch_mode, resolution


def _chirps_v3_file_info(
    year: int,
    month: int | None = None,
    product: str | None = None,
    stream: str | None = None,
    fetch_mode: str | None = None,
    resolution: str | None = None,
) -> tuple[str, str]:
    """Return (filename, url) for a CHIRPS v3 daily NetCDF."""
    product, stream, fetch_mode, resolution = _normalize_chirps_v3_options(
        product, stream, fetch_mode, resolution
    )
    if fetch_mode == "monthly_netcdf":
        if month is None:
            raise ValueError("month is required when chirps_fetch_mode='monthly_netcdf'")
        fname = f"chirps-v3.0.{int(year)}.{int(month):02d}.days_{resolution}.nc"
        url = f"{_CHIRPS_V3_BASE_URL}/{stream}/{product}/netcdf/byMonth/{fname}"
    else:
        fname = f"chirps-v3.0.{product}.{int(year)}.days_{resolution}.nc"
        url = f"{_CHIRPS_V3_BASE_URL}/{stream}/{product}/netcdf/byYear/{fname}"
    return fname, url


def _chirps_v3_cache_path(cache_dir: str, fname: str, product: str,
                          stream: str, fetch_mode: str) -> str:
    subdir = os.path.join(cache_dir, f"v3_{stream}_{product}_{fetch_mode}")
    os.makedirs(subdir, exist_ok=True)
    return os.path.join(subdir, fname)


def _validate_chirps_v3_nc(path: str) -> bool:
    """Return True only if a CHIRPS v3 NetCDF opens and has daily precip data."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    try:
        import xarray as xr
        with xr.open_dataset(path) as ds:
            var = "precip" if "precip" in ds else (list(ds.data_vars)[0] if ds.data_vars else None)
            if var is None or "time" not in ds[var].dims:
                return False
            if ds.sizes.get("time", 0) < 1:
                return False
            sample = ds[var].isel(time=0).load()
            if sample.size == 0:
                return False
        return True
    except Exception as exc:  # noqa: BLE001
        print(f"  Invalid/corrupt CHIRPS v3 NetCDF {os.path.basename(path)}: {exc}")
        return False


def _download_chirps_v3_file(
    year: int,
    month: int | None,
    cache_dir: str,
    product: str,
    stream: str,
    fetch_mode: str,
    resolution: str,
    timeout: int | None = None,
) -> str | None:
    """Download one CHIRPS v3 daily NetCDF into a product-specific cache dir."""
    timeout = int(timeout if timeout is not None else get_config_number(
        "weather.chirps_v3.download_timeout_seconds", 14400
    ))
    import requests

    fname, url = _chirps_v3_file_info(year, month, product, stream, fetch_mode, resolution)
    dest = _chirps_v3_cache_path(cache_dir, fname, product, stream, fetch_mode)
    if os.path.exists(dest) and _validate_chirps_v3_nc(dest):
        return dest
    if os.path.exists(dest):
        print(f"  Removing corrupt cached CHIRPS v3 file: {fname}")
        try:
            os.remove(dest)
        except OSError:
            pass

    tmp = dest + ".part"
    if os.path.exists(tmp):
        try:
            os.remove(tmp)
        except OSError:
            pass
    label = f"{year}-{month:02d}" if month else str(year)
    try:
        print(f"  Downloading CHIRPS v3 {stream}/{product} {label} ({resolution})...")
        with requests.get(url, stream=True, timeout=timeout) as r:
            if r.status_code == 404:
                print(f"  CHIRPS v3 file unavailable: {os.path.basename(urlparse(url).path)}")
                return None
            r.raise_for_status()
            with open(tmp, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    if chunk:
                        fh.write(chunk)
        if not _validate_chirps_v3_nc(tmp):
            raise RuntimeError("downloaded file failed NetCDF validation")
        os.replace(tmp, dest)
        return dest
    except Exception as exc:  # noqa: BLE001
        print(f"  CHIRPS v3 {label} download failed: {exc}")
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass
        return None


def _point_arrays(shapefile, id_col: str, lat_col: str, lon_col: str):
    ids = [str(r[id_col]) for _, r in shapefile.iterrows()]
    if hasattr(shapefile, "geometry") and shapefile.geometry.notna().any():
        try:
            g = shapefile.to_crs("EPSG:4326")
            lons = g.geometry.x.values.astype(float)
            lats = g.geometry.y.values.astype(float)
            return ids, lats, lons
        except Exception:
            pass
    lats = np.array([float(r[lat_col]) for _, r in shapefile.iterrows()])
    lons = np.array([float(r[lon_col]) for _, r in shapefile.iterrows()])
    return ids, lats, lons


def _months_for_range(
    start_year: int,
    end_year: int,
    months: list[int] | tuple[int, ...] | set[int] | None = None,
) -> list[tuple[int, int]]:
    today = date.today()
    month_filter = None
    if months is not None:
        month_filter = {int(m) for m in months}
        bad = [m for m in month_filter if m < 1 or m > 12]
        if bad:
            raise ValueError("chirps_months values must be calendar months 1..12")
    out = []
    for yr in range(int(start_year), int(end_year) + 1):
        last_month = 12
        if yr == today.year:
            last_month = today.month
        elif yr > today.year:
            continue
        for mo in range(1, last_month + 1):
            if month_filter is not None and mo not in month_filter:
                continue
            out.append((yr, mo))
    return out


def _extract_chirps_v3_rain(nc_paths: list[str], ids: list[str],
                            lats: np.ndarray, lons: np.ndarray) -> dict[str, pd.Series]:
    """Extract daily CHIRPS v3 rainfall from already-downloaded NetCDF files."""
    import xarray as xr

    out = {pid: {} for pid in ids}
    pts_lat = xr.DataArray(lats, dims="points")

    for path in nc_paths:
        with xr.open_dataset(path) as ds:
            var = "precip" if "precip" in ds else list(ds.data_vars)[0]
            latname = "latitude" if "latitude" in ds.coords else "lat"
            lonname = "longitude" if "longitude" in ds.coords else "lon"

            ds_lons = np.asarray(ds[lonname].values)
            sample_lons = np.asarray(lons, dtype=float)
            if np.nanmin(ds_lons) >= 0 and np.nanmin(sample_lons) < 0:
                sample_lons = np.where(sample_lons < 0, sample_lons + 360.0, sample_lons)
            pts_lon = xr.DataArray(sample_lons, dims="points")

            sel = ds[var].sel({latname: pts_lat, lonname: pts_lon}, method="nearest")
            vals = np.asarray(sel.values)
            times = pd.to_datetime(ds["time"].values)
            if vals.shape[0] != len(times):
                vals = vals.T
            date_codes = np.asarray([f"{t.year}{t.dayofyear:03d}" for t in times])
            vals = np.where(vals <= _CHIRPS_V3_NODATA, np.nan, vals)
            for j, pid in enumerate(ids):
                col = vals[:, j]
                good = ~np.isnan(col)
                if good.any():
                    out[pid].update(zip(date_codes[good].tolist(),
                                        col[good].astype(float).tolist()))

    return {pid: pd.Series(values, dtype="float64") for pid, values in out.items()}


def _fetch_single_cog_val(url: str, lon: float, lat: float) -> float:
    """Read one COG cell, allowing callers to distinguish fetch errors from nodata."""
    import rasterio

    with rasterio.open(url) as src:
        val = next(src.sample([(lon, lat)]))[0]
        if val <= _CHIRPS_V3_NODATA:
            return np.nan
        return float(val)


def _extract_point_cog(pid: str, lat: float, lon: float, dates: list[date], product: str, stream: str) -> pd.Series:
    results: dict[str, float] = {}
    failures = 0
    urls = []
    for d in dates:
        url = f"{_CHIRPS_V3_BASE_URL}/{stream}/{product}/cogs/{d.year}/chirps-v3.0.{product}.{d.year}.{d.month:02d}.{d.day:02d}.cog"
        urls.append((d, url))

    with ThreadPoolExecutor(max_workers=20) as executor:
        future_to_date = {executor.submit(_fetch_single_cog_val, url, lon, lat): d for d, url in urls}
        for future in future_to_date:
            d = future_to_date[future]
            date_code = f"{d.year}{d.timetuple().tm_yday:03d}"
            try:
                val = future.result()
                results[date_code] = val
            except Exception:
                failures += 1
                results[date_code] = np.nan
    series = pd.Series(results, dtype=float)
    series.attrs["fetch_failures"] = failures
    return series


def _cog_dates_for_year(year: int, months: list[int] | None) -> list[date]:
    months_set = None if months is None else {int(month) for month in months}
    current = date(year, 1, 1)
    end = date(year, 12, 31)
    dates = []
    while current <= end:
        if months_set is None or current.month in months_set:
            dates.append(current)
        current += timedelta(days=1)
    return dates


def _cog_cache_tag(months: list[int] | None) -> str:
    if months is None:
        return "all"
    return "m" + "-".join(f"{month:02d}" for month in sorted({int(m) for m in months}))


def _extract_chirps_v3_rain_remote_cog(
    ids: list[str],
    lats: np.ndarray,
    lons: np.ndarray,
    start_year: int,
    end_year: int,
    product: str,
    stream: str,
    chirps_cache_dir: str,
    months: list[int] | None = None,
) -> dict[str, pd.Series]:
    out = {}
    for pid, lat, lon in zip(ids, lats, lons):
        all_series = []
        for yr in range(start_year, end_year + 1):
            dates_yr = _cog_dates_for_year(yr, months)
            expected_codes = [f"{d.year}{d.timetuple().tm_yday:03d}" for d in dates_yr]
            cache_fn = (
                f"cog_cache_{lat:.5f}_{lon:.5f}_{product}_{stream}_{yr}_"
                f"{_cog_cache_tag(months)}.csv"
            )
            cache_path = os.path.join(chirps_cache_dir, cache_fn)
            if os.path.exists(cache_path):
                try:
                    df_cache = pd.read_csv(cache_path, dtype={"DATE": str})
                    df_cache['DATE'] = df_cache['DATE'].astype(str)
                    if (df_cache['DATE'].tolist() == expected_codes
                            and not df_cache['DATE'].duplicated().any()):
                        s_yr = pd.Series(
                            df_cache['RAIN'].values,
                            index=df_cache['DATE'].values,
                            dtype=float,
                        ).dropna()
                        all_series.append(s_yr)
                        continue
                    print(
                        f"  Warning: Incomplete CHIRPS COG cache for point {pid} "
                        f"year {yr}; re-fetching."
                    )
                except Exception as e:
                    print(f"  Warning: Failed to read cache for point {pid} year {yr}: {e}. Re-fetching.")

            s_yr = _extract_point_cog(pid, lat, lon, dates_yr, product, stream)
            all_series.append(s_yr.dropna())

            if s_yr.attrs.get("fetch_failures", 0) == 0:
                df_cache = pd.DataFrame({'DATE': s_yr.index, 'RAIN': s_yr.values})
                df_cache.to_csv(cache_path, index=False)
            else:
                print(
                    f"  Warning: {s_yr.attrs['fetch_failures']} CHIRPS COG request(s) "
                    f"failed for point {pid}, year {yr}; incomplete results were not cached."
                )

        if all_series:
            s_combined = pd.concat(all_series)
            out[pid] = s_combined
        else:
            out[pid] = pd.Series(dtype=float)
    return out


def _init_gee(project: str | None = None) -> None:
    try:
        import ee
        project = project or get_config_value("weather.chirps_v3.gee_project", "")
        if project:
            ee.Initialize(project=project)
        else:
            ee.Initialize()
    except Exception as exc:
        try:
            import ee
            ee.Initialize(project='earthengine-public')
        except Exception:
            raise RuntimeError(
                "Earth Engine initialization failed. Please run 'earthengine authenticate' "
                "in your terminal/shell to set up GEE credentials, or check your GEE account."
            ) from exc


def _extract_point_gee(pid: str, lat: float, lon: float, year: int, product: str, stream: str) -> pd.Series:
    import ee
    collection_id = f"UCSB-CHC/CHIRPS/V3/DAILY_{product.upper()}"
    col = ee.ImageCollection(collection_id).select('precipitation')
    point = ee.Geometry.Point([float(lon), float(lat)])
    start_date = f"{year}-01-01"
    end_date = f"{year}-12-31"

    info = col.filterDate(start_date, end_date).getRegion(point, scale=5500).getInfo()
    if not info or len(info) <= 1:
        return pd.Series(dtype=float)

    df_pt = pd.DataFrame(info[1:], columns=info[0])
    times = pd.to_datetime(df_pt['time'], unit='ms')
    date_codes = [f"{t.year}{t.dayofyear:03d}" for t in times]
    precip = df_pt['precipitation'].astype(float).values
    precip = np.where(precip <= _CHIRPS_V3_NODATA, np.nan, precip)

    return pd.Series(precip, index=date_codes, dtype=float).dropna()


def _extract_chirps_v3_rain_gee(
    ids: list[str],
    lats: np.ndarray,
    lons: np.ndarray,
    start_year: int,
    end_year: int,
    product: str,
    stream: str,
    chirps_cache_dir: str,
    months: list[int] | None = None,
    gee_project: str | None = None,
) -> dict[str, pd.Series]:
    _init_gee(gee_project)
    out = {}
    for pid, lat, lon in zip(ids, lats, lons):
        all_series = []
        for yr in range(start_year, end_year + 1):
            cache_fn = f"gee_cache_{lat:.5f}_{lon:.5f}_{product}_{stream}_{yr}.csv"
            cache_path = os.path.join(chirps_cache_dir, cache_fn)
            if os.path.exists(cache_path):
                try:
                    df_cache = pd.read_csv(cache_path)
                    df_cache['DATE'] = df_cache['DATE'].astype(str)
                    s_yr = pd.Series(df_cache['RAIN'].values, index=df_cache['DATE'].values, dtype=float)
                    all_series.append(s_yr)
                    continue
                except Exception as e:
                    print(f"  Warning: Failed to read GEE cache for point {pid} year {yr}: {e}. Re-fetching.")

            try:
                s_yr = _extract_point_gee(pid, lat, lon, yr, product, stream)
                all_series.append(s_yr)
                if not s_yr.empty:
                    df_cache = pd.DataFrame({'DATE': s_yr.index, 'RAIN': s_yr.values})
                    df_cache.to_csv(cache_path, index=False)
            except Exception as e:
                print(f"  Error: Failed to fetch GEE data for point {pid} year {yr}: {e}")

        if all_series:
            s_combined = pd.concat(all_series)
            if months is not None:
                months_set = set(months)
                keep = []
                for date_code in s_combined.index:
                    yr_val = int(date_code[:4])
                    doy_val = int(date_code[4:])
                    d_obj = date(yr_val, 1, 1) + timedelta(days=doy_val - 1)
                    if d_obj.month in months_set:
                        keep.append(True)
                    else:
                        keep.append(False)
                s_combined = s_combined[keep]
            out[pid] = s_combined
        else:
            out[pid] = pd.Series(dtype=float)
    return out


def _extract_chirps_v3_rain_cli(
    temp_csv: str,
    start_year: int,
    end_year: int,
    product: str,
    stream: str,
    fetch_mode: str,
    chirps_cache_dir: str,
    temp_out: str,
    gee_project: str | None = None,
) -> None:
    df_pts = pd.read_csv(temp_csv)
    ids = df_pts['id'].astype(str).tolist()
    lats = df_pts['lat'].values
    lons = df_pts['lon'].values

    if fetch_mode == "gee":
        extracted = _extract_chirps_v3_rain_gee(
            ids, lats, lons, start_year, end_year, product, stream, chirps_cache_dir, gee_project=gee_project
        )
    elif fetch_mode == "remote_cog":
        extracted = _extract_chirps_v3_rain_remote_cog(
            ids, lats, lons, start_year, end_year, product, stream, chirps_cache_dir
        )
    else:
        raise ValueError(f"Unsupported fetch_mode for CLI: {fetch_mode}")

    out_dict = {pid: s.to_dict() for pid, s in extracted.items()}
    import json
    with open(temp_out, 'w') as f:
        json.dump(out_dict, f)


def extract_chirps_v3_rainfall(
    shapefile,
    start_year: int,
    end_year: int,
    id_col: str,
    lat_col: str,
    lon_col: str,
    chirps_cache_dir: str,
    product: Literal["rnl", "sat"] | None = None,
    stream: Literal["final", "prelim"] | None = None,
    fetch_mode: Literal["monthly_netcdf", "yearly_netcdf", "gee", "remote_cog"] | None = None,
    resolution: str | None = None,
    months: list[int] | tuple[int, ...] | set[int] | None = None,
    gee_project: str | None = None,
) -> dict[str, dict[str, float]]:
    """Download/extract CHIRPS v3 daily rainfall for points.

    Returns ``{point_id: {"YYYYDOY": rain_mm_day}}``. Points or dates outside
    CHIRPS coverage are omitted so callers can fall back to their base source.
    ``months`` is an optional 1..12 filter for monthly NetCDF extraction; it is
    mainly useful for live smoke tests and partial backfills.
    """
    product, stream, fetch_mode, resolution = _normalize_chirps_v3_options(
        product, stream, fetch_mode, resolution
    )
    if months is not None:
        months = sorted({int(month) for month in months})
        if any(month < 1 or month > 12 for month in months):
            raise ValueError("chirps_months/months values must be calendar months 1..12")
    os.makedirs(chirps_cache_dir, exist_ok=True)
    ids, lats, lons = _point_arrays(shapefile, id_col, lat_col, lon_col)
    in_band = np.abs(lats) <= _CHIRPS_V3_LAT_LIMIT
    if not in_band.any():
        print("  All points outside CHIRPS v3 coverage (|lat| > 60); using fallback rainfall.")
        return {pid: {} for pid in ids}

    extract_ids = [pid for pid, keep in zip(ids, in_band) if keep]
    extract_lats = lats[in_band]
    extract_lons = lons[in_band]
    empty_result = {pid: {} for pid in ids}

    if fetch_mode == "gee":
        extracted = _extract_chirps_v3_rain_gee(
            extract_ids, extract_lats, extract_lons, start_year, end_year,
            product, stream, chirps_cache_dir, months, gee_project=gee_project
        )
        empty_result.update({pid: series.to_dict() for pid, series in extracted.items()})
        return empty_result

    elif fetch_mode == "remote_cog":
        extracted = _extract_chirps_v3_rain_remote_cog(
            extract_ids, extract_lats, extract_lons, start_year, end_year,
            product, stream, chirps_cache_dir, months
        )
        empty_result.update({pid: series.to_dict() for pid, series in extracted.items()})
        return empty_result

    nc_paths: list[str] = []
    if fetch_mode == "monthly_netcdf":
        pieces = _months_for_range(start_year, end_year, months=months)
    else:
        if months is not None:
            raise ValueError("chirps_months/months is only supported with monthly_netcdf")
        current_year = date.today().year
        pieces = [(yr, None) for yr in range(int(start_year), min(int(end_year), current_year) + 1)]

    for yr, mo in pieces:
        path = _download_chirps_v3_file(
            yr, mo, chirps_cache_dir, product, stream, fetch_mode, resolution
        )
        if path:
            nc_paths.append(path)
    if not nc_paths:
        print("  No CHIRPS v3 files available; using fallback rainfall.")
        return {pid: {} for pid in ids}

    print(f"  Extracting CHIRPS v3 {stream}/{product} rainfall for {len(ids)} point(s) "
          f"from {len(nc_paths)} file(s)...")
    extracted = _extract_chirps_v3_rain(
        nc_paths, extract_ids, extract_lats, extract_lons
    )
    empty_result.update({pid: series.to_dict() for pid, series in extracted.items()})
    return empty_result


def _process_single_nasapower_chirps_v3(args: dict) -> None:
    lat = args["latitude"]
    lon = args["longitude"]
    pid = args["point_id"]
    output_dir = args["output_dir"]
    start_date = args["start_date"]
    end_date = args["end_date"]
    log_file = args["log_file"]
    chirps_rain = args["chirps_rain"]
    product = args["chirps_product"]
    stream = args["chirps_stream"]

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
        n_chirps = 0
        if abs(lat) <= _CHIRPS_V3_LAT_LIMIT:
            n_chirps = merge_rainfall_into_weather(df, chirps_rain)
        rain_source = (
            f"CHIRPS-v3 {stream}/{product} where available, {n_chirps} days; "
            "NASA-POWER otherwise"
            if n_chirps else "NASA-POWER (CHIRPS-v3 unavailable here)"
        )

        tav = _calc_tav(df)
        amp = _calc_amp(df)
        header = (
            f"$WEATHER DATA: NASA-POWER + CHIRPS-v3 rain (Point ID: {pid}) "
            f"[{rain_source}]\n"
            f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
            f"  NCV3 {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}   2.0   2.0\n"
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
            lines.append(line.replace(" -99.0", "   -99"))
        with open(out_path, "w") as fh:
            fh.write(header + "\n")
            fh.write("\n".join(lines) + "\n")
    except Exception as exc:  # noqa: BLE001
        msg = (f"\n--- ERROR ---\nFailed: Point ID {pid} | "
               f"Lat {lat:.3f}, Lon {lon:.3f}\nError: {exc}\n")
        print(msg)
        with open(log_file, "a") as lf:
            lf.write(msg)


def process_weather_nasapower_chirps_v3(
    shapefile,
    start_year: int,
    end_year: int,
    output_dir: str,
    id_col: str,
    lat_col: str,
    lon_col: str,
    n_cores: int,
    log_file: str,
    chirps_cache_dir: str,
    chirps_product: Literal["rnl", "sat"] | None = None,
    chirps_stream: Literal["final", "prelim"] | None = None,
    chirps_fetch_mode: Literal["monthly_netcdf", "yearly_netcdf", "gee", "remote_cog"] | None = None,
    chirps_resolution: str | None = None,
    chirps_months: list[int] | tuple[int, ...] | set[int] | None = None,
    chirps_gee_project: str | None = None,
) -> None:
    """Write DSSAT .WTH files using NASA POWER with CHIRPS v3 rainfall.

    NASA POWER supplies all non-rain variables. CHIRPS v3 replaces ``RAIN``
    where daily values are available; dates outside CHIRPS v3 coverage fall back
    to NASA POWER rainfall.
    """
    product, stream, fetch_mode, resolution = _normalize_chirps_v3_options(
        chirps_product, chirps_stream, chirps_fetch_mode, chirps_resolution
    )
    today = date.today()
    current_year = today.year
    start_date = f"{int(start_year)}0101"
    if int(end_year) == current_year:
        safe_end = today - timedelta(days=2)
        end_date = safe_end.strftime("%Y%m%d")
        print(f"End year is current year. Fetching NASA POWER up to: {safe_end.isoformat()}")
    else:
        end_date = f"{int(end_year)}1231"

    print(f"--- Starting NASA-POWER + CHIRPS-v3 {stream}/{product} "
          f"({fetch_mode}, Years: {start_year}-{end_year}) ---")
    os.makedirs(output_dir, exist_ok=True)

    chirps_by_point = extract_chirps_v3_rainfall(
        shapefile=shapefile,
        start_year=start_year,
        end_year=end_year,
        id_col=id_col,
        lat_col=lat_col,
        lon_col=lon_col,
        chirps_cache_dir=chirps_cache_dir,
        product=product,
        stream=stream,
        fetch_mode=fetch_mode,
        resolution=resolution,
        months=chirps_months,
        gee_project=chirps_gee_project,
    )
    ids, lats, lons = _point_arrays(shapefile, id_col, lat_col, lon_col)
    tasks = []
    for pid, lat, lon in zip(ids, lats, lons):
        tasks.append(dict(
            latitude=float(lat), longitude=float(lon), point_id=pid,
            output_dir=output_dir, start_date=start_date, end_date=end_date,
            log_file=log_file, chirps_rain=chirps_by_point.get(pid, {}),
            chirps_product=product, chirps_stream=stream,
        ))

    print(f"Registered {n_cores} cores for parallel NASA-POWER download.")
    with ProcessPoolExecutor(max_workers=n_cores) as pool:
        futures = {pool.submit(_process_single_nasapower_chirps_v3, t): t["point_id"] for t in tasks}
        for fut in as_completed(futures):
            pid = futures[fut]
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                print(f"ERROR (point {pid}): {exc}")

    print(f"\nNASA-POWER + CHIRPS-v3 processing complete. Check '{output_dir}'.\n")
