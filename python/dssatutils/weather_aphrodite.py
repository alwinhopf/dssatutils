# Hybrid weather source: NASA POWER variables + local APHRODITE rainfall NetCDF.
#
# APHRODITE (Asian Precipitation - Highly-Resolved Observational Data) is a dense
# rain-gauge daily product for monsoon Asia / South & SE Asia / Himalaya / Middle
# East (0.25 or 0.5 deg, 1951-). Its core product is precipitation (the separate
# temperature product is daily MEAN only), so this backend mirrors the CHIRPS/
# MSWEP hybrids: T/SRAD/RH/wind come from NASA POWER, rainfall from APHRODITE.

from .weather_format import format_wth_value, wind_run

import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

from .weather_gridded_common import extract_netcdf_series, find_nc_file
from .weather_nasapower import _fetch_nasa_power, _calc_tav, _calc_amp


def _worker(args):
    pid = args["pid"]; lat = args["lat"]; lon = args["lon"]
    out = os.path.join(args["output_dir"], f"{pid}.WTH")
    if os.path.exists(out):
        return
    df = _fetch_nasa_power(lat, lon, args["start"], args["end"]).rename(columns={
        "ALLSKY_SFC_SW_DWN": "SRAD", "T2M_MAX": "TMAX", "T2M_MIN": "TMIN",
        "PRECTOTCORR": "RAIN", "T2MDEW": "TDEW", "RH2M": "RH2M", "WS2M": "WIND",
    })
    rain = pd.Series(args["aphro_rain"], dtype="float64")
    mapped = df["DATE"].map(rain)
    use = mapped.notna()
    df.loc[use, "RAIN"] = mapped[use].values
    tav, amp = _calc_tav(df), _calc_amp(df)
    header = (
        f"$WEATHER DATA: NASA-POWER + APHRODITE rain (Point ID: {pid})\n"
        f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
        f"  APHR {lat:8.4f} {lon:8.4f}   -99 {tav:5.1f} {amp:5.1f}   2.0   2.0\n"
        f"@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND"
    )
    lines = []
    for _, r in df.iterrows():
        line = (f"{r['DATE']:>7s}{format_wth_value(r['SRAD'], 1)}{format_wth_value(r['TMAX'], 1)}{format_wth_value(r['TMIN'], 1)}"
                f"{format_wth_value(r['RAIN'], 1)}{format_wth_value(r['TDEW'], 1)}{format_wth_value(r['RH2M'], 1)}{format_wth_value(wind_run(r['WIND']), 1)}")
        lines.append(line.replace(" -99.0", "   -99"))
    with open(out, "w") as fh:
        fh.write(header + "\n" + "\n".join(lines) + "\n")


def process_weather_aphrodite(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    aphrodite_nc_dir: str,
) -> None:
    """Write DSSAT .WTH using NASA POWER plus local APHRODITE (Asia) rainfall."""
    os.makedirs(output_dir, exist_ok=True)
    pts = shapefile.to_crs("EPSG:4326") if hasattr(shapefile, "geometry") else shapefile
    ids = [str(r[id_col]) for _, r in pts.iterrows()]
    if hasattr(pts, "geometry"):
        lats = pts.geometry.y.values.astype(float); lons = pts.geometry.x.values.astype(float)
    else:
        lats = pts[lat_col].astype(float).values; lons = pts[lon_col].astype(float).values
    path = find_nc_file(aphrodite_nc_dir, ["aphro", "precip", "rain", "pr"])
    if not path:
        raise FileNotFoundError(f"APHRODITE rainfall NetCDF not found in {aphrodite_nc_dir}")
    rain = extract_netcdf_series(path, ["precip", "rain", "pr"], ids, lats, lons,
                                 int(start_year), int(end_year), "rain")
    tasks = []
    for pid, lat, lon in zip(ids, lats, lons):
        tasks.append({"pid": pid, "lat": float(lat), "lon": float(lon),
                      "output_dir": output_dir, "start": f"{int(start_year)}0101",
                      "end": f"{int(end_year)}1231", "aphro_rain": rain.get(pid, {})})
    print(f"--- Starting NASA-POWER + APHRODITE Processing (Years: {start_year}-{end_year}) ---")
    with ProcessPoolExecutor(max_workers=n_cores) as pool:
        for fut in as_completed([pool.submit(_worker, t) for t in tasks]):
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                with open(log_file, "a") as lf:
                    lf.write(f"APHRODITE hybrid worker failed: {exc}\n")
    print(f"\nNASA-POWER + APHRODITE processing complete. Check '{output_dir}'.\n")
