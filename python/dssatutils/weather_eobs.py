# File: weather_eobs.py
# Python port of weather_eobs.R
#
# Weather source: E-OBS (ECA&D European gridded daily observations) -> DSSAT .WTH.
#
# WHY: E-OBS is the standard high-resolution (0.1°, ~11 km) daily gauge-based
# gridded dataset for Europe, 1950-present. Unlike most gridded products it
# includes daily GLOBAL RADIATION (qq), so DSSAT's SRAD comes straight from the
# data rather than being estimated — a good European peer to Daymet/PRISM and a
# gridded complement to the DWD station product.
#
# E-OBS variables (one NetCDF per variable, single data variable inside):
#   tx = max temp (°C)   tn = min temp (°C)   tg = mean temp (°C)
#   rr = precipitation (mm)   qq = global radiation (W/m²)
#   fg = mean wind speed (m/s)   hu = relative humidity (%)
#
# ACCESS — two modes:
#   (A) LOCAL (default, no key): point *eobs_nc_dir* at a folder of pre-downloaded
#       E-OBS NetCDF files (e.g. tx_ens_mean_0.1deg_reg_v*.nc). Download from
#       https://www.ecad.eu/download/ensembles/download.php (free registration).
#   (B) CDS (optional): set eobs_use_cds=True to fetch an area/time SUBSET via the
#       Copernicus CDS dataset "insitu-gridded-observations-europe" (needs a
#       ~/.cdsapirc key + cdsapi, exactly like the AgERA5 module).
#
# The NetCDF extraction (_extract_points) and the .WTH formatting (_write_wth) are
# isolated from any network/credential dependency, so they are unit-testable with
# a synthetic E-OBS-structured dataset (see tests/).

import glob
import math
import os
from datetime import date

import numpy as np
import pandas as pd

# DSSAT var -> (E-OBS variable name, filename token used to locate the file).
_EOBS_VARS = {
    "TMAX": "tx",
    "TMIN": "tn",
    "TMEAN": "tg",
    "RAIN": "rr",
    "SRAD": "qq",     # W/m² -> MJ/m²/day (×0.0864)
    "WIND": "fg",
    "RH2M": "hu",
}
_CDS_DATASET = "insitu-gridded-observations-europe"


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


def _tdew_from_rh(tmean_c, rh_pct):
    """Dew-point (°C) from mean temperature + relative humidity (Magnus)."""
    t = np.asarray(tmean_c, dtype=float)
    rh = np.clip(np.asarray(rh_pct, dtype=float), 1.0, 100.0)
    a, b = 17.625, 243.04
    with np.errstate(invalid="ignore", divide="ignore"):
        gamma = np.log(rh / 100.0) + (a * t) / (b + t)
        td = (b * gamma) / (a - gamma)
    return td


# ---------------------------------------------------------------------------
# File resolution + optional CDS fetch
# ---------------------------------------------------------------------------

def _find_var_file(eobs_nc_dir: str, token: str):
    """Locate the E-OBS NetCDF for a variable token (e.g. 'tx') in a folder."""
    if not eobs_nc_dir or not os.path.isdir(eobs_nc_dir):
        return None
    # Match files whose name starts with '<token>_' or contains '_<token>_'.
    for f in sorted(glob.glob(os.path.join(eobs_nc_dir, "*.nc"))):
        base = os.path.basename(f).lower()
        if base.startswith(f"{token}_") or f"_{token}_" in base or f"{token}_ens" in base:
            return f
    return None


def _download_eobs_cds(token: str, year_start: int, year_end: int, area,
                       cache_dir: str):
    """Optional: fetch an E-OBS variable subset via the Copernicus CDS."""
    import cdsapi
    dest = os.path.join(cache_dir, f"eobs_{token}_{year_start}_{year_end}.nc")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    name_map = {"tx": "maximum_temperature", "tn": "minimum_temperature",
                "tg": "mean_temperature", "rr": "precipitation_amount",
                "qq": "surface_shortwave_downwelling_radiation",
                "fg": "wind_speed", "hu": "relative_humidity"}
    req = {
        "product_type": "ensemble_mean",
        "variable": name_map.get(token, token),
        "grid_resolution": "0.1deg",
        "period": "full_period",
        "version": "30.0e",
        "format": "netcdf",
        "area": list(area),  # N, W, S, E
    }
    try:
        cdsapi.Client().retrieve(_CDS_DATASET, req, dest)
        return dest if os.path.exists(dest) else None
    except Exception as exc:  # noqa: BLE001
        print(f"  E-OBS CDS download failed ({token}): {exc}")
        return None


# ---------------------------------------------------------------------------
# NetCDF point extraction (TESTABLE with a synthetic dataset; no network)
# ---------------------------------------------------------------------------

def _extract_points(path: str, dssat_var: str, ids, lats, lons,
                    year_start: int, year_end: int) -> dict:
    """Extract a per-point daily series ({pid: {YYYYDOY: value}}) for one var."""
    import xarray as xr
    out = {pid: {} for pid in ids}
    with xr.open_dataset(path) as ds:
        ds = ds.load()
        dv = next((v for v in ds.data_vars), None)
        if dv is None:
            return out
        latname = "latitude" if "latitude" in ds.coords else ("lat" if "lat" in ds.coords else None)
        lonname = "longitude" if "longitude" in ds.coords else ("lon" if "lon" in ds.coords else None)
        if latname is None or lonname is None or "time" not in ds.coords:
            return out
        da = ds[dv]
        # Restrict to requested years to keep memory small.
        times = pd.to_datetime(ds["time"].values)
        ysel = (times.year >= year_start) & (times.year <= year_end)
        da = da.isel(time=np.where(ysel)[0])
        times = times[ysel]
        pts_lat = xr.DataArray(np.asarray(lats), dims="points")
        pts_lon = xr.DataArray(np.asarray(lons), dims="points")
        sel = da.sel({latname: pts_lat, lonname: pts_lon}, method="nearest")
        vals = np.asarray(sel.values)
        if vals.ndim == 1:                       # single point -> (time,)
            vals = vals.reshape(len(times), 1) if vals.shape[0] == len(times) else vals.reshape(1, -1).T
        elif vals.shape[0] != len(times):
            vals = vals.T
        # Unit conversion: global radiation W/m² -> MJ/m²/day.
        if dssat_var == "SRAD":
            vals = vals * 0.0864
        date_codes = [f"{t.year}{t.dayofyear:03d}" for t in times]
        for j, pid in enumerate(ids):
            col = vals[:, j]
            good = np.isfinite(col)
            out[pid].update({dc: float(v) for dc, v, g in zip(date_codes, col, good) if g})
    return out


# ---------------------------------------------------------------------------
# .WTH writer (TESTABLE with synthetic data; no network)
# ---------------------------------------------------------------------------

def _write_wth(df: pd.DataFrame, pid: str, lat: float, lon: float,
               output_dir: str) -> str:
    """Write one DSSAT .WTH from a daily DataFrame (E-OBS header/INSI)."""
    tav = _calc_tav(df)
    amp = _calc_amp(df)
    header = (
        f"$WEATHER DATA: E-OBS (Point ID: {pid})\n"
        f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
        f"  EOBS {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}   2.0  10.0\n"
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


def _assemble_and_write(per_var, ids, lats, lons, output_dir, log_file) -> int:
    """Join per-variable series into per-point frames and write .WTH files."""
    written = 0
    for pid, lat, lon in zip(ids, lats, lons):
        try:
            cols = {}
            for dssat_var in _EOBS_VARS:
                series = per_var.get(dssat_var, {}).get(pid, {})
                cols[dssat_var] = pd.Series(series, dtype="float64")
            tmax = cols["TMAX"]
            if tmax.empty:
                raise ValueError("No E-OBS data extracted for this point (outside grid / land mask).")
            frame = pd.DataFrame(cols)
            frame.index.name = "DATE"
            frame = frame.reset_index()
            dts = pd.to_datetime(frame["DATE"], format="%Y%j")
            frame["YEAR"] = dts.dt.year
            frame["MM"] = dts.dt.month
            # Derive TDEW from mean temp + RH where both exist, else -99.
            if "TMEAN" in frame and "RH2M" in frame:
                frame["TDEW"] = _tdew_from_rh(frame["TMEAN"].values, frame["RH2M"].values)
            else:
                frame["TDEW"] = -99.0
            for need in ("SRAD", "RAIN", "RH2M", "WIND", "TDEW"):
                if need not in frame:
                    frame[need] = -99.0
            frame = frame[frame["TMAX"].notna() & frame["TMIN"].notna()].reset_index(drop=True)
            frame = frame.fillna(-99)
            _write_wth(frame, pid, lat, lon, output_dir)
            written += 1
        except Exception as exc:  # noqa: BLE001
            msg = f"\n--- ERROR ---\nE-OBS point {pid} ({lat:.3f},{lon:.3f}): {exc}\n"
            print(msg)
            with open(log_file, "a") as lf:
                lf.write(msg)
    return written


# ---------------------------------------------------------------------------
# Public entry point (mirrors the other weather sources)
# ---------------------------------------------------------------------------

def process_weather_eobs(
    shapefile,           # GeoDataFrame
    start_year: int,
    end_year: int,
    output_dir: str,
    id_col: str,
    lat_col: str,
    lon_col: str,
    n_cores: int,        # kept for API compatibility (extraction is serial I/O)
    log_file: str,
    eobs_nc_dir: str = "",
    eobs_cache_dir: str = "",
    eobs_use_cds: bool = False,
) -> None:
    """Build DSSAT .WTH from E-OBS for every grid point.

    Provide pre-downloaded E-OBS NetCDFs via *eobs_nc_dir* (default), or set
    *eobs_use_cds=True* to fetch an area subset through the Copernicus CDS
    (requires ~/.cdsapirc + cdsapi). TMAX/TMIN/RAIN are required; SRAD comes from
    E-OBS global radiation (qq); WIND/RH/TDEW are filled where available, else -99.
    """
    os.makedirs(output_dir, exist_ok=True)
    end_year = min(end_year, date.today().year)

    pts = shapefile.copy()
    if hasattr(pts, "geometry"):
        pts = pts.to_crs("EPSG:4326")
        pts[lat_col] = pts.geometry.y
        pts[lon_col] = pts.geometry.x
    ids = [str(r[id_col]) for _, r in pts.iterrows()]
    lats = np.array([float(r[lat_col]) for _, r in pts.iterrows()])
    lons = np.array([float(r[lon_col]) for _, r in pts.iterrows()])

    print(f"--- Starting E-OBS Processing (Years: {start_year}–{end_year}) ---")

    # Resolve a NetCDF file path per variable (local folder or CDS subset).
    var_paths = {}
    if eobs_use_cds:
        if not eobs_cache_dir:
            eobs_cache_dir = os.path.join(output_dir, "eobs_cache")
        os.makedirs(eobs_cache_dir, exist_ok=True)
        pad = 0.2
        area = [float(lats.max() + pad), float(lons.min() - pad),
                float(lats.min() - pad), float(lons.max() + pad)]
        print("  E-OBS via Copernicus CDS (requires ~/.cdsapirc + cdsapi).")
        for dssat_var, token in _EOBS_VARS.items():
            var_paths[dssat_var] = _download_eobs_cds(token, start_year, end_year, area, eobs_cache_dir)
    else:
        if not eobs_nc_dir or not os.path.isdir(eobs_nc_dir):
            raise FileNotFoundError(
                "E-OBS local mode needs eobs_nc_dir pointing at a folder of E-OBS "
                "NetCDF files (tx/tn/rr/qq...). Download from www.ecad.eu, or set "
                "eobs_use_cds=True for the Copernicus CDS subset path.")
        for dssat_var, token in _EOBS_VARS.items():
            var_paths[dssat_var] = _find_var_file(eobs_nc_dir, token)

    if not var_paths.get("TMAX") or not var_paths.get("TMIN"):
        raise FileNotFoundError("E-OBS requires at least tx (TMAX) and tn (TMIN) NetCDF files.")

    per_var = {}
    for dssat_var, path in var_paths.items():
        if not path:
            if dssat_var in ("SRAD", "WIND", "RH2M", "TMEAN", "RAIN"):
                print(f"  E-OBS: no file for {dssat_var} ({_EOBS_VARS[dssat_var]}); it will be written as -99.")
            continue
        per_var[dssat_var] = _extract_points(path, dssat_var, ids, lats, lons, start_year, end_year)

    written = _assemble_and_write(per_var, ids, lats, lons, output_dir, log_file)
    print(f"\nE-OBS processing complete: {written}/{len(ids)} points written to '{output_dir}'.\n")
