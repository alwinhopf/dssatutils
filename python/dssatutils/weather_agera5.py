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

from .credentials import (
    cdsapirc_candidates as _shared_cdsapirc_candidates,
    make_cds_client as _shared_make_cds_client,
    read_cdsapirc as _shared_read_cdsapirc,
)

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
_CDS_TIMESERIES_DATASET = "sis-agrometeorological-indicators-timeseries"
_AGERA5_CDS_REQUEST_CAP = 4
_AGERA5_TIMESERIES_MAX_EXTENT_DEG = 5.0
_AGERA5_TIMESERIES_DEFAULT_CHUNK_DEG = 4.5
_AGERA5_TIMESERIES_PAD_DEG = 0.2
_AGERA5_TIMESERIES_VARS = {
    "TMAX": ("2m_temperature_24_hour_maximum", "Temperature_Air_2m_Max_24h"),
    "TMIN": ("2m_temperature_24_hour_minimum", "Temperature_Air_2m_Min_24h"),
    "SRAD": ("solar_radiation_flux", "Solar_Radiation_Flux"),
    "RAIN": ("precipitation_flux", "Precipitation_Flux"),
    "TDEW": ("2m_dewpoint_temperature_24_hour_mean", "Dew_Point_Temperature_2m_Mean_24h"),
    "RH2M": ("2m_relative_humidity_at_15_00", "Relative_Humidity_2m_15h"),
    "WIND": ("10m_wind_speed_24_hour_mean", "Wind_Speed_10m_Mean_24h"),
}


# ---------------------------------------------------------------------------
# CDS credential helpers
# ---------------------------------------------------------------------------

def _cdsapirc_candidates():
    yield from (str(p) for p in _shared_cdsapirc_candidates())


def _read_cdsapirc():
    return _shared_read_cdsapirc()


def _make_cds_client(cdsapi):
    """Create a cdsapi client from env vars or a discovered .cdsapirc file."""
    return _shared_make_cds_client(cdsapi)


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
    None on failure. Requires cdsapi + a configured CDS token.
    """
    import cdsapi  # imported lazily so the module loads without the key/pkg

    area_tag = "_".join(_slug_float(v) for v in area)
    tag = f"{cds_var}_{sel_value or 'na'}_{year}_{area_tag}"
    dest = os.path.join(cache_dir, f"agera5_{tag}.zip")
    if _valid_agera5_download(dest):
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
        partial = dest + ".partial"
        if os.path.exists(partial):
            os.remove(partial)
        _make_cds_client(cdsapi).retrieve(_CDS_DATASET, req, partial)
        if not _valid_agera5_download(partial):
            raise ValueError("CDS response is not a valid, non-empty AgERA5 archive")
        os.replace(partial, dest)
        return dest
    except Exception as exc:  # noqa: BLE001
        print(f"  AgERA5 download failed ({tag}): {exc}")
        return None


def _valid_agera5_download(path: str) -> bool:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    if path.endswith((".zip", ".partial")):
        try:
            with zipfile.ZipFile(path) as archive:
                members = [m for m in archive.infolist()
                           if not m.is_dir() and m.filename.lower().endswith(".nc")]
                return bool(members) and all(m.file_size > 0 for m in members)
        except (OSError, zipfile.BadZipFile):
            return False
    return True


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


def _slug_float(value: float) -> str:
    return f"{float(value):.4f}".replace("-", "m").replace(".", "p")


def _date_bounds_for_year(year: int):
    start = date(int(year), 1, 1)
    end = date(int(year), 12, 31)
    # AgERA5 is near-real-time but not truly same-day. The CDS form exposes a
    # moving max date; keep current-year requests away from unavailable future
    # dates while preserving complete past years.
    latest_safe = date.today() - timedelta(days=10)
    if end > latest_safe:
        end = latest_safe
    if end < start:
        return None
    return start.isoformat(), end.isoformat()


def _split_agera5_timeseries_chunks(lats, lons, chunk_degrees=None, pad=None):
    """Return point-assignment chunks plus padded CDS request areas.

    The time-series CDS form allows a maximum area of 5 x 5 degrees. Points are
    assigned to non-overlapping unpadded chunks, while each CDS request is padded
    slightly so nearest-cell sampling remains covered.
    """
    pad = _AGERA5_TIMESERIES_PAD_DEG if pad is None else float(pad)
    raw_chunk = (_AGERA5_TIMESERIES_DEFAULT_CHUNK_DEG if chunk_degrees is None
                 else float(chunk_degrees))
    max_raw = max(0.1, _AGERA5_TIMESERIES_MAX_EXTENT_DEG - 2 * pad)
    step = min(max(0.1, raw_chunk), max_raw)
    lat_min, lat_max = float(np.min(lats)), float(np.max(lats))
    lon_min, lon_max = float(np.min(lons)), float(np.max(lons))
    chunks = []
    south = lat_min
    while south <= lat_max + 1e-12:
        north = min(lat_max, south + step)
        west = lon_min
        while west <= lon_max + 1e-12:
            east = min(lon_max, west + step)
            mask = ((lats >= south - 1e-12) & (lats <= north + 1e-12) &
                    (lons >= west - 1e-12) & (lons <= east + 1e-12))
            idx = np.where(mask)[0]
            if len(idx):
                area = [
                    min(90.0, north + pad),
                    max(-180.0, west - pad),
                    max(-90.0, south - pad),
                    min(179.9, east + pad),
                ]
                chunks.append({
                    "idx": idx,
                    "bounds": (south, west, north, east),
                    "area": area,
                })
            if east >= lon_max:
                break
            west = east
        if north >= lat_max:
            break
        south = north
    return chunks


def _agera5_timeseries_cache_path(cache_dir: str, year: int, area, data_format: str) -> str:
    ext = "csv" if data_format == "csv" else "nc"
    n, w, s, e = area
    tag = "_".join(_slug_float(v) for v in (n, w, s, e))
    return os.path.join(cache_dir, f"agera5_timeseries_{year}_{tag}.{ext}")


def _download_agera5_timeseries(year: int, area, cache_dir: str,
                                data_format: str = "csv"):
    """Download one all-variable AgERA5 time-series chunk for a year."""
    data_format = str(data_format or "csv").lower()
    if data_format != "csv":
        raise ValueError("AgERA5 time-series backend currently supports data_format='csv'.")

    bounds = _date_bounds_for_year(year)
    if bounds is None:
        return None
    dest = _agera5_timeseries_cache_path(cache_dir, year, area, data_format)
    if _valid_agera5_timeseries_csv(dest):
        return dest

    import cdsapi  # imported lazily so the module loads without the key/pkg

    req = {
        "variable": [v[0] for v in _AGERA5_TIMESERIES_VARS.values()],
        "date": [bounds[0], bounds[1]],
        "data_format": data_format,
        "area": list(area),
    }
    try:
        partial = dest + ".partial"
        if os.path.exists(partial):
            os.remove(partial)
        _make_cds_client(cdsapi).retrieve(_CDS_TIMESERIES_DATASET, req, partial)
        if not _valid_agera5_timeseries_csv(partial):
            raise ValueError("CDS response is not a valid AgERA5 time-series CSV")
        os.replace(partial, dest)
        return dest
    except Exception as exc:  # noqa: BLE001
        print(f"  AgERA5 time-series download failed ({year}, area={area}): {exc}")
        return None


def _valid_agera5_timeseries_csv(path: str) -> bool:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    try:
        sample = pd.read_csv(path, nrows=5)
    except Exception:  # noqa: BLE001
        return False
    required = {"valid_time", "latitude", "longitude"}
    return bool(len(sample)) and required.issubset({str(c).lower() for c in sample.columns})


def _find_timeseries_column(df: pd.DataFrame, expected: str) -> str:
    if expected in df.columns:
        return expected
    lowered = {str(c).lower(): c for c in df.columns}
    if expected.lower() in lowered:
        return lowered[expected.lower()]
    raise KeyError(f"AgERA5 time-series CSV missing expected column {expected!r}")


def _read_agera5_timeseries_csv(path: str) -> pd.DataFrame:
    raw = pd.read_csv(path)
    date_col = _find_timeseries_column(raw, "valid_time")
    lat_col = _find_timeseries_column(raw, "latitude")
    lon_col = _find_timeseries_column(raw, "longitude")
    out = pd.DataFrame({
        "valid_time": pd.to_datetime(raw[date_col]),
        "latitude": raw[lat_col].astype(float),
        "longitude": raw[lon_col].astype(float),
    })
    for dssat_var, (_, col) in _AGERA5_TIMESERIES_VARS.items():
        values = raw[_find_timeseries_column(raw, col)].astype(float).to_numpy()
        if dssat_var in ("TMAX", "TMIN", "TDEW"):
            values = values - 273.15
        elif dssat_var == "SRAD":
            values = values * 1e-6
        out[dssat_var] = values
    out["DATE"] = [f"{d.year}{d.dayofyear:03d}" for d in out["valid_time"]]
    return out


def _add_timeseries_chunk_to_points(path: str, point_indices, ids, lats, lons,
                                    point_series):
    df = _read_agera5_timeseries_csv(path)
    grids = df[["latitude", "longitude"]].drop_duplicates().reset_index(drop=True)
    for j in point_indices:
        dist = ((grids["latitude"].to_numpy() - lats[j]) ** 2 +
                (grids["longitude"].to_numpy() - lons[j]) ** 2)
        nearest = grids.iloc[int(np.argmin(dist))]
        sub = df[
            np.isclose(df["latitude"], nearest["latitude"]) &
            np.isclose(df["longitude"], nearest["longitude"])
        ].sort_values("valid_time")
        pid = ids[j]
        for dssat_var in _AGERA5_TIMESERIES_VARS:
            point_series[pid][dssat_var].update(
                {dc: float(v) for dc, v in zip(sub["DATE"], sub[dssat_var])
                 if pd.notna(v)}
            )


def _process_weather_agera5_timeseries(
    shapefile,
    start_year: int,
    end_year: int,
    output_dir: str,
    id_col: str,
    lat_col: str,
    lon_col: str,
    n_cores: int,
    log_file: str,
    agera5_cache_dir: str,
    agera5_data_format: str = "csv",
    agera5_timeseries_chunk_degrees: float = _AGERA5_TIMESERIES_DEFAULT_CHUNK_DEG,
) -> None:
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(agera5_cache_dir, exist_ok=True)

    ids = [str(r[id_col]) for _, r in shapefile.iterrows()]
    lats = np.array([float(r[lat_col]) for _, r in shapefile.iterrows()])
    lons = np.array([float(r[lon_col]) for _, r in shapefile.iterrows()])
    end_year = min(end_year, date.today().year)
    chunks = _split_agera5_timeseries_chunks(
        lats, lons, chunk_degrees=agera5_timeseries_chunk_degrees
    )

    print(f"--- Starting AgERA5 Time-Series Download (Years: {start_year}-{end_year}) ---")
    print("  Backend: sis-agrometeorological-indicators-timeseries "
          f"({len(chunks)} area chunk(s), format={agera5_data_format}).")

    point_series = {pid: {v: {} for v in _AGERA5_TIMESERIES_VARS} for pid in ids}
    jobs = [(year, chunk) for year in range(start_year, end_year + 1)
            for chunk in chunks if _date_bounds_for_year(year) is not None]

    try:
        requested_workers = int(n_cores)
    except Exception:  # noqa: BLE001
        requested_workers = 1
    workers = max(1, min(requested_workers, _AGERA5_CDS_REQUEST_CAP, len(jobs) or 1))
    print(
        f"  AgERA5 time-series cache/download phase: {len(jobs)} year-area job(s); "
        f"using {workers} concurrent CDS request(s) (cap={_AGERA5_CDS_REQUEST_CAP})."
    )

    def _dl(job):
        year, chunk = job
        return (
            year,
            chunk,
            _download_agera5_timeseries(
                year, chunk["area"], agera5_cache_dir, agera5_data_format
            ),
        )

    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for year, chunk, path in pool.map(_dl, jobs):
            if not path:
                msg = f"  AgERA5 time-series missing ({year}, area={chunk['area']})"
                print(msg)
                if log_file:
                    with open(log_file, "a") as lf:
                        lf.write(msg + "\n")
                continue
            try:
                _add_timeseries_chunk_to_points(path, chunk["idx"], ids, lats, lons, point_series)
            except Exception as exc:  # noqa: BLE001
                msg = f"  AgERA5 time-series parse failed ({path}): {exc}"
                print(msg)
                if log_file:
                    with open(log_file, "a") as lf:
                        lf.write(msg + "\n")

    written = 0
    for pid, lat, lon in zip(ids, lats, lons):
        try:
            series = {v: pd.Series(point_series[pid][v], dtype="float64")
                      for v in _AGERA5_TIMESERIES_VARS}
            if series["TMAX"].empty:
                raise ValueError("No AgERA5 time-series data extracted for this point.")
            frame = pd.DataFrame(series)
            frame.index.name = "DATE"
            frame = frame.reset_index().sort_values("DATE")
            dts = pd.to_datetime(frame["DATE"], format="%Y%j")
            frame["YEAR"] = dts.dt.year
            frame["MM"] = dts.dt.month
            frame = frame.fillna(-99)
            _write_wth(frame, pid, lat, lon, output_dir)
            written += 1
        except Exception as exc:  # noqa: BLE001
            msg = f"\n--- ERROR ---\nAgERA5 time-series point {pid} ({lat:.3f},{lon:.3f}): {exc}\n"
            print(msg)
            if log_file:
                with open(log_file, "a") as lf:
                    lf.write(msg)

    print(f"\nAgERA5 time-series processing complete: {written}/{len(ids)} points "
          f"written to '{output_dir}'.\n")


# ---------------------------------------------------------------------------
# .WTH writer (TESTABLE with synthetic data; no network)
# ---------------------------------------------------------------------------

def _write_wth(df: pd.DataFrame, pid: str, lat: float, lon: float,
               output_dir: str) -> str:
    """Write one DSSAT .WTH from a daily DataFrame.

    *df* must contain DATE, SRAD, TMAX, TMIN, RAIN, TDEW, RH2M, WIND. Returns
    the output path. Shared formatting with the NASA POWER / Open-Meteo writers.
    """
    climatology = df.copy()
    climatology[["TMAX", "TMIN"]] = climatology[["TMAX", "TMIN"]].where(
        climatology[["TMAX", "TMIN"]] > -90
    )
    tav = _calc_tav(climatology)
    amp = _calc_amp(climatology)
    if not np.isfinite(tav) or not np.isfinite(amp):
        raise ValueError("No valid AgERA5 temperature climatology for point.")
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
    agera5_backend: str = "gridded",
    agera5_data_format: str = "csv",
    agera5_timeseries_chunk_degrees: float = _AGERA5_TIMESERIES_DEFAULT_CHUNK_DEG,
) -> None:
    """Download AgERA5 over the grid's bounding box and write DSSAT .WTH files.

    Requires a configured CDS API key (~/.cdsapirc) plus cdsapi + xarray. Data
    is downloaded once per variable-year (subset to the grid bbox) into
    *agera5_cache_dir* and reused. Unit conversions: temperature K→°C (−273.15),
    solar radiation J/m²/day → MJ/m²/day (×1e-6).
    """
    backend = str(agera5_backend or "gridded").lower().replace("-", "_")
    if backend in ("timeseries", "time_series", "ts"):
        return _process_weather_agera5_timeseries(
            shapefile=shapefile,
            start_year=start_year,
            end_year=end_year,
            output_dir=output_dir,
            id_col=id_col,
            lat_col=lat_col,
            lon_col=lon_col,
            n_cores=n_cores,
            log_file=log_file,
            agera5_cache_dir=agera5_cache_dir,
            agera5_data_format=agera5_data_format,
            agera5_timeseries_chunk_degrees=agera5_timeseries_chunk_degrees,
        )
    if backend not in ("gridded", "grid", "classic"):
        raise ValueError("agera5_backend must be 'gridded' or 'timeseries'.")

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
    #    instead of paying them one after another. Cap concurrency at 4 to avoid
    #    hammering the CDS per-user active-request queue.
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

    try:
        requested_workers = int(n_cores)
    except Exception:  # noqa: BLE001
        requested_workers = 1
    workers = max(1, min(requested_workers, _AGERA5_CDS_REQUEST_CAP, len(jobs)))
    print(
        f"  AgERA5 cache/download phase: {len(jobs)} variable-year job(s); "
        f"using {workers} concurrent CDS request(s) (cap={_AGERA5_CDS_REQUEST_CAP})."
    )

    paths = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
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
