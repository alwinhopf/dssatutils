# Internal helpers for local/cache-backed gridded weather sources.

import glob
import os
import re
from datetime import date

import numpy as np
import pandas as pd
from typing import Optional


def calc_tav(df: pd.DataFrame) -> float:
    return float(((df["TMAX"] + df["TMIN"]) / 2.0).mean())


def calc_amp(df: pd.DataFrame) -> float:
    d = df.copy()
    d["TAVG"] = (d["TMAX"] + d["TMIN"]) / 2.0
    monthly = d.groupby(["YEAR", "MM"])["TAVG"].mean().reset_index()
    annual = monthly.groupby("YEAR")["TAVG"].agg(lambda x: x.max() - x.min())
    return float(annual.mean())


def tdew_from_rh(tmean_c, rh_pct):
    t = np.asarray(tmean_c, dtype=float)
    rh = np.clip(np.asarray(rh_pct, dtype=float), 1.0, 100.0)
    a, b = 17.625, 243.04
    with np.errstate(invalid="ignore", divide="ignore"):
        gamma = np.log(rh / 100.0) + (a * t) / (b + t)
        return (b * gamma) / (a - gamma)


def write_wth(df: pd.DataFrame, pid: str, lat: float, lon: float,
              output_dir: str, source_label: str, insi: str,
              refht: float = 2.0, wndht: float = 2.0) -> str:
    tav = calc_tav(df)
    amp = calc_amp(df)
    header = (
        f"$WEATHER DATA: {source_label} (Point ID: {pid})\n"
        f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
        f"  {insi:<4s} {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}"
        f" {refht:5.1f} {wndht:5.1f}\n"
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
    os.makedirs(output_dir, exist_ok=True)
    out = os.path.join(output_dir, f"{pid}.WTH")
    with open(out, "w") as fh:
        fh.write(header + "\n")
        fh.write("\n".join(lines) + "\n")
    return out


def find_nc_files(nc_dir: str, tokens) -> list[str]:
    if not nc_dir or not os.path.isdir(nc_dir):
        return []
    tokens = [str(t).lower() for t in tokens]
    matches = []
    for f in sorted(glob.glob(os.path.join(nc_dir, "*.nc"))):
        stem = os.path.splitext(os.path.basename(f).lower())[0]
        components = set(filter(None, re.split(r"[^a-z0-9]+", stem)))
        if any(t in components or stem == t for t in tokens):
            matches.append(f)
    return matches


def find_nc_file(nc_dir: str, tokens) -> Optional[str]:
    """Compatibility wrapper returning the first exact-token match."""
    matches = find_nc_files(nc_dir, tokens)
    return matches[0] if matches else None


def _pick_var(ds, aliases):
    aliases = [a.lower() for a in aliases]
    lower = {v.lower(): v for v in ds.data_vars}
    for a in aliases:
        if a in lower:
            return lower[a]
    for v in ds.data_vars:
        lv = v.lower()
        if any(a in lv for a in aliases):
            return v
    return next(iter(ds.data_vars), None)


def _coord_names(ds):
    lat = next((c for c in ("lat", "latitude", "y") if c in ds.coords), None)
    lon = next((c for c in ("lon", "longitude", "x") if c in ds.coords), None)
    return lat, lon


def convert_units(values, units: str, kind: str):
    units_l = (units or "").lower().replace("**", "^")
    arr = np.asarray(values, dtype=float)
    if kind == "temp":
        if "k" in units_l or np.nanmedian(arr) > 100:
            arr = arr - 273.15
    elif kind == "rain":
        if "s-1" in units_l or "/s" in units_l:
            arr = arr * 86400.0
    elif kind == "srad":
        if "w" in units_l and "m" in units_l:
            arr = arr * 0.0864
        elif "j" in units_l and "m" in units_l:
            arr = arr / 1_000_000.0
    elif kind == "wind":
        # Reanalysis/gridded wind is reported at 10 m; DSSAT wants 2 m.
        # FAO-56 log-profile factor u2 = u10 * 4.87 / ln(67.8*10 - 5.42) ~ 0.748.
        arr = arr * 0.748
    elif kind == "vp":
        # Vapour pressure (hPa) -> dewpoint (degC) via the inverse Magnus formula.
        e = np.clip(arr, 1e-3, None)
        ln = np.log(e / 6.1094)
        arr = (243.04 * ln) / (17.625 - ln)
    return arr


def extract_netcdf_series(path, aliases, ids, lats, lons,
                          start_year: int, end_year: int, kind: str) -> dict:
    import xarray as xr

    out = {pid: {} for pid in ids}
    paths = [path] if isinstance(path, (str, os.PathLike)) else list(path)
    datasets = [xr.open_dataset(p) for p in paths]
    try:
        ds = xr.concat(datasets, dim="time").sortby("time") if len(datasets) > 1 else datasets[0]
        var = _pick_var(ds, aliases)
        latname, lonname = _coord_names(ds)
        if var is None or latname is None or lonname is None or "time" not in ds.coords:
            return out
        da = ds[var]
        units = str(da.attrs.get("units", ""))
        times = pd.to_datetime(ds["time"].values)
        keep = (times.year >= start_year) & (times.year <= end_year)
        if not keep.any():
            return out
        da = da.isel(time=np.where(keep)[0])
        times = times[keep]
        qlons = np.asarray(lons, dtype=float)
        grid_lons = np.asarray(ds[lonname].values, dtype=float)
        if np.nanmin(grid_lons) >= 0 and np.nanmax(qlons) <= 180:
            qlons = np.where(qlons < 0, qlons + 360, qlons)
        pts_lat = xr.DataArray(np.asarray(lats, dtype=float), dims="points")
        pts_lon = xr.DataArray(qlons, dims="points")
        sel = da.sel({latname: pts_lat, lonname: pts_lon}, method="nearest")
        vals = np.asarray(sel.values)
        if vals.ndim == 1:
            vals = vals.reshape(len(times), 1)
        elif vals.shape[0] != len(times):
            vals = vals.T
        vals = convert_units(vals, units, kind)
        date_codes = [f"{t.year}{t.dayofyear:03d}" for t in times]
        for j, pid in enumerate(ids):
            col = vals[:, j]
            good = np.isfinite(col)
            out[pid].update({dc: float(v) for dc, v, g in zip(date_codes, col, good) if g})
    finally:
        for dataset in datasets:
            dataset.close()
    return out


def process_local_netcdf_weather(shapefile, start_year, end_year, output_dir,
                                 id_col, lat_col, lon_col, log_file,
                                 nc_dir, var_specs, source_label, insi,
                                 refht=2.0, wndht=2.0) -> int:
    if not nc_dir or not os.path.isdir(nc_dir):
        raise FileNotFoundError(f"{source_label} needs a local NetCDF directory: {nc_dir}")
    os.makedirs(output_dir, exist_ok=True)
    end_year = min(int(end_year), date.today().year)

    pts = shapefile.copy()
    if hasattr(pts, "geometry"):
        pts = pts.to_crs("EPSG:4326")
        pts[lat_col] = pts.geometry.y
        pts[lon_col] = pts.geometry.x
    ids = [str(r[id_col]) for _, r in pts.iterrows()]
    lats = np.array([float(r[lat_col]) for _, r in pts.iterrows()])
    lons = np.array([float(r[lon_col]) for _, r in pts.iterrows()])

    per_var = {}
    for dssat_var, spec in var_specs.items():
        paths = find_nc_files(nc_dir, spec["tokens"])
        if not paths:
            if spec.get("required", False):
                raise FileNotFoundError(f"{source_label} required variable {dssat_var} not found in {nc_dir}")
            print(f"  {source_label}: no NetCDF for {dssat_var}; writing -99 where needed.")
            continue
        per_var[dssat_var] = extract_netcdf_series(
            paths, spec.get("aliases", spec["tokens"]), ids, lats, lons,
            start_year, end_year, spec["kind"])

    for required in ("TMAX", "TMIN", "RAIN", "SRAD"):
        if required not in per_var:
            raise FileNotFoundError(
                f"{source_label} requires {required}; refusing to write a WTH with missing forcing"
            )

    written = 0
    for pid, lat, lon in zip(ids, lats, lons):
        try:
            cols = {}
            for dssat_var in var_specs:
                cols[dssat_var] = pd.Series(per_var.get(dssat_var, {}).get(pid, {}), dtype="float64")
            expected_end = pd.Timestamp(end_year, 12, 31)
            if end_year == date.today().year:
                expected_end = pd.Timestamp(date.today())
            expected = {f"{d.year}{d.dayofyear:03d}" for d in
                        pd.date_range(pd.Timestamp(start_year, 1, 1), expected_end, freq="D")}
            for required in ("TMAX", "TMIN", "RAIN", "SRAD"):
                actual = set(cols[required].dropna().index.astype(str))
                missing = expected - actual
                if missing:
                    raise ValueError(
                        f"{required} is incomplete for requested period ({len(missing)} missing day(s))"
                    )
            frame = pd.DataFrame(cols)
            frame.index.name = "DATE"
            frame = frame.reset_index()
            dts = pd.to_datetime(frame["DATE"], format="%Y%j")
            frame["YEAR"] = dts.dt.year
            frame["MM"] = dts.dt.month
            if "TDEW" not in frame:
                if "TMEAN" in frame and "RH2M" in frame:
                    frame["TDEW"] = tdew_from_rh(frame["TMEAN"].values, frame["RH2M"].values)
                else:
                    frame["TDEW"] = -99.0
            for need in ("SRAD", "RAIN", "RH2M", "WIND"):
                if need not in frame:
                    frame[need] = -99.0
            frame = frame[frame["TMAX"].notna() & frame["TMIN"].notna()].fillna(-99)
            write_wth(frame, pid, lat, lon, output_dir, source_label, insi, refht, wndht)
            written += 1
        except Exception as exc:  # noqa: BLE001
            msg = f"\n--- ERROR ---\n{source_label} point {pid} ({lat:.3f},{lon:.3f}): {exc}\n"
            print(msg)
            with open(log_file, "a") as lf:
                lf.write(msg)
    return written
