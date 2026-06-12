# File: weather_xavier.py
# Python port of weather_xavier.R
#
# Weather source: BR-DWGD (Xavier et al.) Brazilian Daily Weather Gridded Data
# -> DSSAT .WTH.
#
# WHY: BR-DWGD is the standard high-resolution (0.1°) daily gauge-interpolated
# weather product for BRAZIL, 1961-present. It already provides daily solar
# radiation in MJ/m^2/day (no estimation needed) plus Tmax/Tmin, precipitation,
# relative humidity and 2 m wind — exactly DSSAT's inputs.
#
# Variables (one NetCDF per variable; single data variable inside):
#   Tmax (°C)  Tmin (°C)  pr (mm)  Rs (MJ/m^2/day)  RH (%)  u2 (m/s @2m)  ETo (mm)
#
# ACCESS (open, NO key): download the NetCDF files (per variable, partitioned in
# ~20-year blocks) from the BR-DWGD site and point *xavier_nc_dir* at the folder:
#   https://sites.google.com/site/alexandrecandidoxavierufes/brazilian-daily-weather-gridded-data
# (also mirrored as a Google Earth Engine ImageCollection). Coverage: Brazil.
#
# The NetCDF extraction (_extract_points) and .WTH writer (_write_wth) are
# isolated from any download so they are unit-testable with synthetic data.

import glob
import os
from datetime import date

import numpy as np
import pandas as pd

# DSSAT var -> (E-OBS-style filename token / variable name in the Xavier files).
_XAVIER_VARS = {
    "TMAX": "Tmax",
    "TMIN": "Tmin",
    "RAIN": "pr",
    "SRAD": "Rs",      # already MJ/m^2/day — no conversion
    "RH2M": "RH",
    "WIND": "u2",      # 2 m wind
}


def _calc_tav(df):
    return float(((df["TMAX"] + df["TMIN"]) / 2.0).mean())


def _calc_amp(df):
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
        return (b * gamma) / (a - gamma)


def _find_var_files(nc_dir, token):
    """All Xavier NetCDFs for a variable token (filenames may be year-partitioned)."""
    if not nc_dir or not os.path.isdir(nc_dir):
        return []
    hits = []
    for f in sorted(glob.glob(os.path.join(nc_dir, "*.nc"))):
        base = os.path.basename(f)
        # Xavier names start with the variable, e.g. 'Tmax_...nc' / 'pr_...nc'.
        if base.startswith(f"{token}_") or f"_{token}_" in base or base == f"{token}.nc":
            hits.append(f)
    return hits


def _extract_points(paths, dssat_var, ids, lats, lons, year_start, year_end):
    """Per-point daily series ({pid: {YYYYDOY: value}}) across one variable's files."""
    import xarray as xr
    out = {pid: {} for pid in ids}
    pts_lat = xr.DataArray(np.asarray(lats), dims="points")
    pts_lon = xr.DataArray(np.asarray(lons), dims="points")
    for path in paths:
        with xr.open_dataset(path) as ds:
            ds = ds.load()
            dv = next((v for v in ds.data_vars), None)
            if dv is None or "time" not in ds.coords:
                continue
            latname = "latitude" if "latitude" in ds.coords else ("lat" if "lat" in ds.coords else None)
            lonname = "longitude" if "longitude" in ds.coords else ("lon" if "lon" in ds.coords else None)
            if latname is None or lonname is None:
                continue
            times = pd.to_datetime(ds["time"].values)
            ysel = (times.year >= year_start) & (times.year <= year_end)
            if not ysel.any():
                continue
            da = ds[dv].isel(time=np.where(ysel)[0])
            times = times[ysel]
            sel = da.sel({latname: pts_lat, lonname: pts_lon}, method="nearest")
            vals = np.asarray(sel.values)
            if vals.ndim == 1:
                vals = vals.reshape(len(times), 1) if vals.shape[0] == len(times) else vals.reshape(1, -1).T
            elif vals.shape[0] != len(times):
                vals = vals.T
            date_codes = [f"{t.year}{t.dayofyear:03d}" for t in times]
            for j, pid in enumerate(ids):
                col = vals[:, j]
                good = np.isfinite(col)
                out[pid].update({dc: float(v) for dc, v, g in zip(date_codes, col, good) if g})
    return out


def _write_wth(df, pid, lat, lon, output_dir):
    tav = _calc_tav(df)
    amp = _calc_amp(df)
    header = (
        f"$WEATHER DATA: BR-DWGD/Xavier (Point ID: {pid})\n"
        f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
        f"  XAVR {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}   2.0   2.0\n"
        f"@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND"
    )
    lines = []
    for _, row in df.iterrows():
        line = (
            f"{row['DATE']:>7s}"
            f"{row['SRAD']:6.1f}{row['TMAX']:6.1f}{row['TMIN']:6.1f}"
            f"{row['RAIN']:6.1f}{row['TDEW']:6.1f}{row['RH2M']:6.1f}{row['WIND']:6.1f}"
        )
        line = line.replace(" -99.0", "   -99")
        lines.append(line)
    out_path = os.path.join(output_dir, f"{pid}.WTH")
    with open(out_path, "w") as fh:
        fh.write(header + "\n")
        fh.write("\n".join(lines) + "\n")
    return out_path


def process_weather_xavier(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    xavier_nc_dir="",
):
    """Build DSSAT .WTH from the BR-DWGD/Xavier Brazilian grids for each point.

    Point *xavier_nc_dir* at the folder of pre-downloaded Xavier NetCDFs
    (Tmax/Tmin/pr/Rs/RH/u2). TMAX/TMIN/RAIN required; SRAD is taken directly
    (already MJ/m^2/day); TDEW is derived from mean temperature + RH.
    """
    os.makedirs(output_dir, exist_ok=True)
    end_year = min(end_year, date.today().year)
    if not xavier_nc_dir or not os.path.isdir(xavier_nc_dir):
        raise FileNotFoundError(
            "Xavier needs xavier_nc_dir pointing at a folder of BR-DWGD NetCDF files "
            "(Tmax/Tmin/pr/Rs/RH/u2). Download from the BR-DWGD site (no key).")

    pts = shapefile.copy()
    if hasattr(pts, "geometry"):
        pts = pts.to_crs("EPSG:4326")
        pts[lat_col] = pts.geometry.y
        pts[lon_col] = pts.geometry.x
    ids = [str(r[id_col]) for _, r in pts.iterrows()]
    lats = np.array([float(r[lat_col]) for _, r in pts.iterrows()])
    lons = np.array([float(r[lon_col]) for _, r in pts.iterrows()])

    print(f"--- Starting BR-DWGD/Xavier Processing (Years: {start_year}-{end_year}) ---")
    var_files = {v: _find_var_files(xavier_nc_dir, tok) for v, tok in _XAVIER_VARS.items()}
    if not var_files["TMAX"] or not var_files["TMIN"]:
        raise FileNotFoundError("Xavier requires at least Tmax and Tmin NetCDF files.")

    per_var = {}
    for v, paths in var_files.items():
        if not paths:
            print(f"  Xavier: no file for {v} ({_XAVIER_VARS[v]}); it will be written as -99.")
            continue
        per_var[v] = _extract_points(paths, v, ids, lats, lons, start_year, end_year)

    written = 0
    for pid, lat, lon in zip(ids, lats, lons):
        try:
            cols = {v: pd.Series(per_var.get(v, {}).get(pid, {}), dtype="float64") for v in _XAVIER_VARS}
            if cols["TMAX"].empty:
                raise ValueError("No Xavier data extracted (point outside the Brazil grid).")
            frame = pd.DataFrame(cols)
            frame.index.name = "DATE"
            frame = frame.reset_index()
            dts = pd.to_datetime(frame["DATE"], format="%Y%j")
            frame["YEAR"] = dts.dt.year
            frame["MM"] = dts.dt.month
            tmean = (frame["TMAX"] + frame["TMIN"]) / 2.0
            frame["TDEW"] = _tdew_from_rh(tmean.values, frame["RH2M"].values) if "RH2M" in frame else -99.0
            for need in ("SRAD", "RAIN", "RH2M", "WIND", "TDEW"):
                if need not in frame:
                    frame[need] = -99.0
            frame = frame[frame["TMAX"].notna() & frame["TMIN"].notna()].reset_index(drop=True)
            frame = frame.fillna(-99)
            _write_wth(frame, pid, lat, lon, output_dir)
            written += 1
        except Exception as exc:  # noqa: BLE001
            msg = f"\n--- ERROR ---\nXavier point {pid} ({lat:.3f},{lon:.3f}): {exc}\n"
            print(msg)
            with open(log_file, "a") as lf:
                lf.write(msg)
    print(f"\nXavier processing complete: {written}/{len(ids)} points written to '{output_dir}'.\n")
