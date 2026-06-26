# Weather source: TerraClimate monthly climate, disaggregated to daily DSSAT input.
#
# TerraClimate (Abatzoglou et al. 2018) is a MONTHLY product (~4 km), so it cannot
# be written straight to a daily DSSAT .WTH. This backend expands each month into
# continuous daily records: Tmax/Tmin/SRAD/wind are held constant at the monthly
# value and the monthly precipitation TOTAL is spread uniformly across the days of
# the month. The result is a runnable but deliberately smooth weather file —
# intended for screening / climatology only, NOT for analyses that depend on
# day-to-day variability (rainfall intensity, dry spells, heat-stress days).

import calendar

import numpy as np
import pandas as pd

from .weather_gridded_common import (
    extract_netcdf_series, find_nc_file, write_wth,
)

# DSSAT var -> (filename tokens, NetCDF var aliases, unit-conversion kind)
_VARS = {
    "TMAX": (["tmmx", "tmax"], ["tmmx", "tmax"], "temp"),
    "TMIN": (["tmmn", "tmin"], ["tmmn", "tmin"], "temp"),
    "RAIN": (["ppt", "precip", "pr"], ["ppt", "precip", "pr"], "rain"),
    "SRAD": (["srad", "rsds"], ["srad", "rsds"], "srad"),
    "WIND": (["ws", "wind"], ["ws", "wind"], "wind"),
}


def _coords(shapefile, id_col, lat_col, lon_col):
    pts = shapefile.copy()
    if hasattr(pts, "geometry"):
        pts = pts.to_crs("EPSG:4326")
        lats = pts.geometry.y.values.astype(float)
        lons = pts.geometry.x.values.astype(float)
    else:
        lats = pts[lat_col].astype(float).values
        lons = pts[lon_col].astype(float).values
    ids = [str(r[id_col]) for _, r in pts.iterrows()]
    return ids, lats, lons


def _expand_month_to_daily(monthly_by_var):
    """monthly_by_var: {DSSAT_VAR: {YYYYDOY(month-start): value}} -> daily DataFrame.

    Returns continuous daily rows for every month for which TMAX and TMIN exist.
    """
    months = sorted(set(monthly_by_var.get("TMAX", {})) & set(monthly_by_var.get("TMIN", {})))
    rows = []
    for code in months:
        year = int(code[:4])
        month = (pd.Timestamp(f"{year}-01-01") + pd.Timedelta(days=int(code[4:]) - 1)).month
        ndays = calendar.monthrange(year, month)[1]
        tmax = monthly_by_var["TMAX"][code]
        tmin = monthly_by_var["TMIN"][code]
        srad = monthly_by_var.get("SRAD", {}).get(code, -99.0)
        wind = monthly_by_var.get("WIND", {}).get(code, -99.0)
        ppt_total = monthly_by_var.get("RAIN", {}).get(code, np.nan)
        rain_daily = (ppt_total / ndays) if np.isfinite(ppt_total) else -99.0
        for dom in range(1, ndays + 1):
            d = pd.Timestamp(year=year, month=month, day=dom)
            rows.append({
                "DATE": f"{d.year}{d.dayofyear:03d}", "YEAR": d.year, "MM": d.month,
                "SRAD": srad, "TMAX": tmax, "TMIN": tmin, "RAIN": rain_daily,
                "TDEW": -99.0, "RH2M": -99.0, "WIND": wind,
            })
    return pd.DataFrame(rows)


def process_weather_terraclimate(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    terraclimate_nc_dir: str,
) -> None:
    """Write screening-grade daily DSSAT .WTH from TerraClimate monthly NetCDF.

    TerraClimate is monthly; this backend disaggregates each month to continuous
    daily records (constant Tmax/Tmin/SRAD/wind, monthly precip spread evenly).
    Use for screening/climatology only — it carries no day-to-day variability.
    """
    import os
    if not terraclimate_nc_dir or not os.path.isdir(terraclimate_nc_dir):
        raise FileNotFoundError(
            f"TerraClimate needs a local NetCDF directory: {terraclimate_nc_dir}")
    os.makedirs(output_dir, exist_ok=True)
    print("WARNING: TerraClimate is monthly; disaggregated to daily for screening/climatology only.")
    print(f"--- Starting TerraClimate Processing (Years: {start_year}-{end_year}) ---")

    ids, lats, lons = _coords(shapefile, id_col, lat_col, lon_col)

    per_var = {}
    for dssat_var, (tokens, aliases, kind) in _VARS.items():
        path = find_nc_file(terraclimate_nc_dir, tokens)
        if not path:
            if dssat_var in ("TMAX", "TMIN", "RAIN"):
                raise FileNotFoundError(
                    f"TerraClimate required variable {dssat_var} not found in {terraclimate_nc_dir}")
            continue
        per_var[dssat_var] = extract_netcdf_series(
            path, aliases, ids, lats, lons, int(start_year), int(end_year), kind)

    written = 0
    for pid, lat, lon in zip(ids, lats, lons):
        try:
            monthly_by_var = {v: per_var.get(v, {}).get(pid, {}) for v in _VARS}
            frame = _expand_month_to_daily(monthly_by_var)
            if frame.empty:
                raise ValueError("No overlapping monthly TMAX/TMIN extracted for this point.")
            write_wth(frame, pid, lat, lon, output_dir, "TerraClimate monthly->daily", "TCLM")
            written += 1
        except Exception as exc:  # noqa: BLE001
            msg = f"\n--- ERROR ---\nTerraClimate point {pid} ({lat:.3f},{lon:.3f}): {exc}\n"
            print(msg)
            with open(log_file, "a") as lf:
                lf.write(msg)
    print(f"\nTerraClimate processing complete: {written} point(s) written.\n")
