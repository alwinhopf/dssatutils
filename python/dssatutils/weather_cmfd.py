# File: weather_cmfd.py
# Python port of weather_cmfd.R
#
# Weather source: CMFD (China Meteorological Forcing Dataset) -> DSSAT .WTH.
#
# WHY: CMFD is the standard high-resolution (0.1°) near-surface forcing dataset
# for CHINA, 1979-2018, fusing station, reanalysis and satellite data. It is
# 3-HOURLY, so this module aggregates each variable to the daily statistics DSSAT
# needs (Tmax/Tmin from sub-daily temperature, daily precip total, mean solar
# radiation, etc.).
#
# CMFD variables (one NetCDF per variable per month; var name == token):
#   temp (K, 2 m air T)   prec (mm/hr, precip rate)   srad (W/m^2, dwn shortwave)
#   shum (kg/kg, specific humidity)   pres (Pa)   wind (m/s, 10 m)
#
# Daily reduction:  TMAX/TMIN = max/min of temp-273.15;  RAIN = mean(prec)*24;
#   SRAD = mean(srad)*0.0864 MJ/m^2/day;  WIND = mean(wind);
#   RH/TDEW from daily-mean shum + temp + pres.
#
# ACCESS: free, but requires a (free) account at the National Tibetan Plateau
# Data Center (data.tpdc.ac.cn) to DOWNLOAD the NetCDF files. Point *cmfd_nc_dir*
# at the downloaded folder. Coverage: China.

import glob
import os
from datetime import date

import numpy as np
import pandas as pd

# DSSAT-relevant CMFD variables + filename/variable token.
_CMFD_TOKENS = {"temp": "temp", "prec": "prec", "srad": "srad",
                "shum": "shum", "pres": "pres", "wind": "wind"}


def _calc_tav(df):
    return float(((df["TMAX"] + df["TMIN"]) / 2.0).mean())


def _calc_amp(df):
    d = df.copy()
    d["TAVG"] = (d["TMAX"] + d["TMIN"]) / 2.0
    monthly = d.groupby(["YEAR", "MM"])["TAVG"].mean().reset_index()
    annual = monthly.groupby("YEAR")["TAVG"].agg(lambda x: x.max() - x.min())
    return float(annual.mean())


def _rh_tdew_from_shum(shum, temp_c, pres_pa):
    """Relative humidity (%) and dew point (°C) from specific humidity."""
    q = np.asarray(shum, dtype=float)
    p_hpa = np.asarray(pres_pa, dtype=float) / 100.0
    tc = np.asarray(temp_c, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        e = q * p_hpa / (0.622 + 0.378 * q)            # actual vapour pressure (hPa)
        es = 6.112 * np.exp(17.67 * tc / (tc + 243.5))  # saturation VP (hPa)
        rh = np.clip(100.0 * e / es, 0.0, 100.0)
        ln = np.log(np.where(e > 0, e / 6.112, np.nan))
        tdew = (243.5 * ln) / (17.67 - ln)
    return rh, tdew


def _find_var_files(nc_dir, token):
    if not nc_dir or not os.path.isdir(nc_dir):
        return []
    hits = []
    for f in sorted(glob.glob(os.path.join(nc_dir, "*.nc"))):
        base = os.path.basename(f).lower()
        if base.startswith(f"{token}_") or f"_{token}_" in base or f"{token}-" in base:
            hits.append(f)
    return hits


def _extract_daily(paths, token, ids, lats, lons):
    """Return {pid: DataFrame(date-indexed daily reductions)} for one variable."""
    import xarray as xr
    parts = {pid: [] for pid in ids}
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
            sel = ds[dv].sel({latname: pts_lat, lonname: pts_lon}, method="nearest")
            vals = np.asarray(sel.values)
            if vals.ndim == 1:
                vals = vals.reshape(len(times), 1) if vals.shape[0] == len(times) else vals.reshape(1, -1).T
            elif vals.shape[0] != len(times):
                vals = vals.T
            for j, pid in enumerate(ids):
                s = pd.Series(vals[:, j], index=times)
                parts[pid].append(s)
    out = {}
    for pid in ids:
        if not parts[pid]:
            out[pid] = pd.DataFrame()
            continue
        s = pd.concat(parts[pid]).sort_index()
        s = s[~s.index.duplicated()]
        g = s.groupby(s.index.normalize())
        if token == "temp":                       # K -> daily max/min/mean (°C)
            daily = pd.DataFrame({"TMAX": g.max() - 273.15, "TMIN": g.min() - 273.15,
                                  "TMEAN": g.mean() - 273.15})
        elif token == "prec":                     # mm/hr -> daily total mm
            daily = pd.DataFrame({"RAIN": g.mean() * 24.0})
        elif token == "srad":                     # W/m^2 -> MJ/m^2/day
            daily = pd.DataFrame({"SRAD": g.mean() * 0.0864})
        elif token == "wind":
            daily = pd.DataFrame({"WIND": g.mean()})
        elif token == "shum":
            daily = pd.DataFrame({"SHUM": g.mean()})
        elif token == "pres":
            daily = pd.DataFrame({"PRES": g.mean()})
        else:
            daily = pd.DataFrame({token.upper(): g.mean()})
        out[pid] = daily
    return out


def _write_wth(df, pid, lat, lon, output_dir):
    tav = _calc_tav(df)
    amp = _calc_amp(df)
    header = (
        f"$WEATHER DATA: CMFD (Point ID: {pid})\n"
        f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
        f"  CMFD {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}   2.0  10.0\n"
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


def process_weather_cmfd(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    cmfd_nc_dir="",
):
    """Build DSSAT .WTH from CMFD 3-hourly grids (aggregated to daily) per point.

    Point *cmfd_nc_dir* at the folder of pre-downloaded CMFD NetCDFs
    (temp/prec/srad/shum/pres/wind). temp is required; RH/TDEW are derived from
    specific humidity + pressure + temperature when those files are present.
    """
    os.makedirs(output_dir, exist_ok=True)
    end_year = min(end_year, date.today().year)
    if not cmfd_nc_dir or not os.path.isdir(cmfd_nc_dir):
        raise FileNotFoundError(
            "CMFD needs cmfd_nc_dir pointing at a folder of CMFD NetCDF files "
            "(temp/prec/srad/shum/pres/wind). Download from the TPDC (free account).")

    pts = shapefile.copy()
    if hasattr(pts, "geometry"):
        pts = pts.to_crs("EPSG:4326")
        pts[lat_col] = pts.geometry.y
        pts[lon_col] = pts.geometry.x
    ids = [str(r[id_col]) for _, r in pts.iterrows()]
    lats = np.array([float(r[lat_col]) for _, r in pts.iterrows()])
    lons = np.array([float(r[lon_col]) for _, r in pts.iterrows()])

    print(f"--- Starting CMFD Processing (Years: {start_year}-{end_year}) ---")
    var_files = {tok: _find_var_files(cmfd_nc_dir, tok) for tok in _CMFD_TOKENS}
    if not var_files["temp"]:
        raise FileNotFoundError("CMFD requires at least the temp NetCDF files.")

    daily = {tok: _extract_daily(paths, tok, ids, lats, lons)
             for tok, paths in var_files.items() if paths}

    written = 0
    for pid, lat, lon in zip(ids, lats, lons):
        try:
            tdf = daily.get("temp", {}).get(pid, pd.DataFrame())
            if tdf.empty:
                raise ValueError("No CMFD temperature data extracted (point outside the China grid).")
            frame = tdf.copy()
            for tok, col in (("prec", "RAIN"), ("srad", "SRAD"), ("wind", "WIND")):
                d = daily.get(tok, {}).get(pid)
                frame = frame.join(d[[col]], how="left") if (d is not None and not d.empty) else frame
            # RH / TDEW from specific humidity + pressure + daily-mean temp.
            shum = daily.get("shum", {}).get(pid)
            pres = daily.get("pres", {}).get(pid)
            if shum is not None and not shum.empty and pres is not None and not pres.empty:
                j = frame.join(shum, how="left").join(pres, how="left")
                rh, tdew = _rh_tdew_from_shum(j["SHUM"].values, j["TMEAN"].values, j["PRES"].values)
                frame["RH2M"] = rh
                frame["TDEW"] = tdew
            frame = frame.reset_index().rename(columns={"index": "DT"})
            frame["DT"] = pd.to_datetime(frame["DT"])
            frame = frame[(frame["DT"].dt.year >= start_year) & (frame["DT"].dt.year <= end_year)]
            if frame.empty:
                raise ValueError("No CMFD data in the requested year range.")
            frame["DATE"] = frame["DT"].apply(lambda t: f"{t.year}{t.dayofyear:03d}")
            frame["YEAR"] = frame["DT"].dt.year
            frame["MM"] = frame["DT"].dt.month
            for need in ("SRAD", "RAIN", "RH2M", "WIND", "TDEW"):
                if need not in frame:
                    frame[need] = -99.0
            frame = frame[frame["TMAX"].notna() & frame["TMIN"].notna()].fillna(-99)
            _write_wth(frame, pid, lat, lon, output_dir)
            written += 1
        except Exception as exc:  # noqa: BLE001
            msg = f"\n--- ERROR ---\nCMFD point {pid} ({lat:.3f},{lon:.3f}): {exc}\n"
            print(msg)
            with open(log_file, "a") as lf:
                lf.write(msg)
    print(f"\nCMFD processing complete: {written}/{len(ids)} points written to '{output_dir}'.\n")
