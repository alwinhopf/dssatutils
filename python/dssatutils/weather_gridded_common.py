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
    matches = [v for v in ds.data_vars if any(a in v.lower() for a in aliases)]
    if len(matches) > 1:
        raise ValueError(f"Ambiguous weather variable aliases {aliases}: {matches}")
    return matches[0] if matches else None


def _coord_names(ds):
    lat = next((c for c in ("lat", "latitude") if c in ds.coords), None)
    lon = next((c for c in ("lon", "longitude") if c in ds.coords), None)
    return lat, lon


def convert_units(values, units: str, kind: str, wind_height_m=None):
    """Convert declared daily-forcing units; reject ambiguous units/height."""
    u = re.sub(r"[\s_^*]", "", (units or "").lower()).replace("per", "/")
    arr = np.asarray(values, dtype=float)
    def unsupported():
        raise ValueError(f"Unsupported or missing {kind} units: {units!r}")
    if kind == "temp":
        if u in {"k", "kelvin", "degk", "degreek", "degreeskelvin"}:
            return arr - 273.15
        if u in {"c", "degc", "celsius", "degreecelsius", "degreescelsius", "°c"}:
            return arr
        unsupported()
    if kind == "rain":
        if u in {"mm", "mm/day", "mmday-1", "mmd-1", "kgm-2", "kg/m2"}: return arr
        if u in {"kgm-2s-1", "kg/m2/s", "mms-1", "mm/s"}: return arr * 86400
        if u in {"m", "m/day", "mday-1"}: return arr * 1000
        unsupported()
    if kind == "srad":
        if u in {"wm-2", "w/m2"}: return arr * 0.0864
        for prefix, scale in (("mj", 1.0), ("kj", 0.001), ("j", 1e-6)):
            if u in {prefix + tail for tail in ("m-2", "/m2", "m-2day-1", "m-2d-1", "/m2/day", "m-2/day")}:
                return arr * scale
        unsupported()
    if kind == "wind":
        if u in {"km/h", "kmh-1"}: arr = arr / 3.6
        elif u not in {"m/s", "ms-1"}: unsupported()
        if wind_height_m is None or not np.isfinite(float(wind_height_m)) or float(wind_height_m) < 0.1:
            raise ValueError("Wind requires an explicit height_m in metres")
        h = float(wind_height_m)
        return arr if h == 2 else arr * np.log(67.8 * 2 - 5.42) / np.log(67.8 * h - 5.42)
    if kind == "vp":
        if u == "pa": arr = arr / 100
        elif u == "kpa": arr = arr * 10
        elif u not in {"hpa", "mbar"}: unsupported()
        e = np.where(arr > 0, arr, np.nan)
        ln = np.log(e / 6.1094)
        return (243.04 * ln) / (17.625 - ln)
    if kind in {"rh", "rh2m"}:
        if u in {"1", "fraction"}: return arr * 100
        if u in {"%", "percent", "percentage"}: return arr
        unsupported()
    return arr


def extract_netcdf_series(path, aliases, ids, lats, lons,
                          start_year: int, end_year: int, kind: str, wind_height_m=None) -> dict:
    import xarray as xr

    out = {pid: {} for pid in ids}
    paths = [path] if isinstance(path, (str, os.PathLike)) else list(path)
    datasets = [xr.open_dataset(p) for p in paths]
    try:
        ds = xr.concat(datasets, dim="time").sortby("time") if len(datasets) > 1 else datasets[0]
        var = _pick_var(ds, aliases)
        if var is None:
            raise KeyError(
                f"No matching variable for aliases {aliases} found in dataset. "
                f"Available variables: {list(ds.data_vars)}"
            )
        latname, lonname = _coord_names(ds)
        if latname is None or lonname is None or "time" not in ds.coords:
            raise ValueError(
                f"Dataset missing required coordinates (lat, lon, time). Found: {list(ds.coords)}"
            )
        da = ds[var]
        for coordinate, accepted in ((latname, {"degrees_north", "degree_north", "degrees_n", "degree_n"}),
                                     (lonname, {"degrees_east", "degree_east", "degrees_e", "degree_e"})):
            axis = ds[coordinate]
            axis_units = str(axis.attrs.get("units", "")).lower()
            if axis.ndim != 1 or (axis_units and axis_units not in accepted):
                raise ValueError("Weather coordinates must be one-dimensional geographic degrees")
        if da.attrs.get("grid_mapping"):
            mapping = ds.get(da.attrs["grid_mapping"])
            if mapping is not None and mapping.attrs.get("grid_mapping_name", "latitude_longitude") != "latitude_longitude":
                raise ValueError("Projected weather grids must be reprojected to geographic latitude/longitude")
        units = str(da.attrs.get("units", ""))
        times = pd.to_datetime(ds["time"].values)
        keep = (times.year >= start_year) & (times.year <= end_year)
        if not keep.any():
            return out
        da = da.isel(time=np.where(keep)[0])
        times = times[keep]
        # These adapters consume daily forcing. Without interval bounds and
        # accumulation metadata, subdaily integration is ambiguous.
        if times.has_duplicates:
            raise ValueError("Duplicate weather timestamps; remove overlapping records before extraction")
        if times.normalize().has_duplicates:
            raise ValueError("Subdaily weather data is unsupported; provide one daily record per date")
        qlons = np.asarray(lons, dtype=float).copy()
        grid_lons = np.asarray(ds[lonname].values, dtype=float)
        grid_lats = np.asarray(ds[latname].values, dtype=float)
        if np.nanmin(grid_lons) >= 0:
            qlons = np.where(qlons < 0, qlons + 360, qlons)
        elif np.nanmax(grid_lons) <= 180:
            qlons = np.where(qlons > 180, qlons - 360, qlons)

        min_lat, max_lat = float(np.nanmin(grid_lats)), float(np.nanmax(grid_lats))
        min_lon, max_lon = float(np.nanmin(grid_lons)), float(np.nanmax(grid_lons))
        lat_res = 0.0
        lon_res = 0.0
        for lat_val, lon_val, pid in zip(lats, qlons, ids):
            if not (min_lat - lat_res <= lat_val <= max_lat + lat_res) or \
               not (min_lon - lon_res <= lon_val <= max_lon + lon_res):
                raise ValueError(
                    f"Point {pid} ({lat_val:.4f}, {lon_val:.4f}) is outside grid domain: "
                    f"lat [{min_lat:.4f}, {max_lat:.4f}], lon [{min_lon:.4f}, {max_lon:.4f}]"
                )

        pts_lat = xr.DataArray(np.asarray(lats, dtype=float), dims="points")
        pts_lon = xr.DataArray(qlons, dims="points")
        sel = da.sel({latname: pts_lat, lonname: pts_lon}, method="nearest")
        vals = np.asarray(sel.values)
        if vals.ndim == 1:
            vals = vals.reshape(len(times), 1)
        elif vals.shape[0] != len(times):
            vals = vals.T
        height = wind_height_m if wind_height_m is not None else da.attrs.get("height_m", da.attrs.get("height"))
        vals = convert_units(vals, units, kind, wind_height_m=height)
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
            start_year, end_year, spec["kind"], wind_height_m=spec.get("height_m"))

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
