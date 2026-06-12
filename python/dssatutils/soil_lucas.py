# File: soil_lucas.py
# Python port of soil_lucas.R
#
# Soil source: LUCAS Topsoil (via ESDAC, European Soil Data Centre) -> DSSAT .SOL.
#
# WHY: LUCAS is the harmonised, ground-measured topsoil survey for the EU
# (tens of thousands of georeferenced samples: texture, organic carbon, pH, N,
# CaCO3, ...). For European sites it is real lab-measured soil rather than a
# global model prediction.
#
# IMPORTANT LIMITATION: LUCAS samples the TOPSOIL only (0-20 cm). DSSAT needs a
# profile, so this module writes the measured 0-20 cm layer and carries those
# properties down to LUCAS_ROOTING_MAX_CM as an explicit extrapolation (flagged
# in the .SOL header). Treat the subsoil as approximate.
#
# Bulk density is NOT in the standard LUCAS table; it is estimated from the
# Saxton & Rawls saturation (porosity) as BD = (1 - SSAT) * 2.65 g/cm^3.
#
# ACCESS: free, but ESDAC distributes LUCAS behind a one-off request form
# (https://esdac.jrc.ec.europa.eu/ -> LUCAS Topsoil). Download the point table
# (CSV/XLSX), then point *lucas_csv* at it. Coverage: EU.

import math
import os
from typing import Optional

import numpy as np
import pandas as pd

from .soil_ssurgo import _failure, _saxton_rawls

_LUCAS_ROOTING_MAX_CM = 150

# Column aliases across LUCAS releases (2009 / 2015 / 2018). First match wins.
_ALIASES = {
    "id":   ["POINTID", "POINT_ID", "Point_ID", "id", "ID"],
    "lat":  ["TH_LAT", "GPS_LAT", "lat", "Latitude", "latitude", "Y", "POINT_Y"],
    "lon":  ["TH_LONG", "GPS_LONG", "lon", "Longitude", "longitude", "X", "POINT_X"],
    "clay": ["clay", "Clay", "Clay_content", "clay_content"],
    "sand": ["sand", "Sand", "Sand_content", "sand_content"],
    "silt": ["silt", "Silt", "Silt_content", "silt_content"],
    "oc":   ["OC", "oc", "OC_gkg", "organic_carbon"],     # g/kg
    "ph":   ["pH_H2O", "pH_in_H2O", "pH_CaCl2", "pH", "ph"],
}


def _resolve_cols(df: pd.DataFrame, col_map: Optional[dict]) -> dict:
    cols = {}
    user = col_map or {}
    for key, aliases in _ALIASES.items():
        if key in user and user[key] in df.columns:
            cols[key] = user[key]
            continue
        cols[key] = next((a for a in aliases if a in df.columns), None)
    return cols


def _load_lucas(lucas_csv: str, col_map: Optional[dict]) -> pd.DataFrame:
    if lucas_csv.lower().endswith((".xlsx", ".xls")):
        raw = pd.read_excel(lucas_csv)
    else:
        raw = pd.read_csv(lucas_csv, sep=None, engine="python")
    cols = _resolve_cols(raw, col_map)
    for need in ("lat", "lon", "clay", "sand"):
        if cols[need] is None:
            raise ValueError(
                f"LUCAS table is missing a '{need}' column (looked for {_ALIASES[need]}). "
                f"Pass col_map={{'{need}': '<your column>'}} to override.")
    out = pd.DataFrame({
        "src_id": raw[cols["id"]].astype(str) if cols["id"] else [f"L{i}" for i in range(len(raw))],
        "lat": pd.to_numeric(raw[cols["lat"]], errors="coerce"),
        "lon": pd.to_numeric(raw[cols["lon"]], errors="coerce"),
        "clay": pd.to_numeric(raw[cols["clay"]], errors="coerce"),
        "sand": pd.to_numeric(raw[cols["sand"]], errors="coerce"),
    })
    out["silt"] = pd.to_numeric(raw[cols["silt"]], errors="coerce") if cols["silt"] else (100 - out["clay"] - out["sand"])
    out["oc"] = pd.to_numeric(raw[cols["oc"]], errors="coerce") if cols["oc"] else np.nan
    out["ph"] = pd.to_numeric(raw[cols["ph"]], errors="coerce") if cols["ph"] else np.nan
    out = out.dropna(subset=["lat", "lon", "clay", "sand"]).reset_index(drop=True)
    if out.empty:
        raise ValueError("No usable LUCAS rows after parsing (check column mapping / units).")
    return out


def _write_sol(profile: pd.DataFrame, output_dir: str) -> None:
    soil_id = str(profile["ID"].iloc[0])
    lat = profile["latitude"].iloc[0]
    lon = profile["longitude"].iloc[0]
    path = os.path.join(output_dir, f"{soil_id}.SOL")
    if os.path.exists(path):
        return
    lines = [
        "*SOILS: Europe LUCAS Topsoil Profiles",
        "! Generated from LUCAS topsoil (0-20 cm MEASURED; subsoil EXTRAPOLATED), Saxton & Rawls",
        "",
        f"*{soil_id:<6s}  LUCAS         {lat:9.3f} {lon:9.3f}",
        "@SITE        COUNTRY          LAT     LONG SCS FAMILY",
        f" {soil_id:<11s} EU          {lat:9.3f} {lon:9.3f} ",
        "@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE",
        "    BN   .13     6    .6    73     1     1 IB001 IB001 IB001",
        "@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC",
    ]
    for _, layer in profile.sort_values("depth_bottom").iterrows():
        def _f3(v):
            s = f"{v:5.3f}"
            return (" " + s[1:]) if s.startswith("0.") else s
        slll, sdul, ssat = _f3(layer["SLLL"]), _f3(layer["SDUL"]), _f3(layer["SSAT"])
        depth = int(layer["depth_bottom"])
        om_sloc = layer["om_pct"] / 1.724
        lines.append(
            f"{depth:6d}   -99 {slll} {sdul} {ssat}  1.00   -99"
            f" {layer['bulk_density']:5.2f} {om_sloc:5.2f}"
            f" {layer['clay_pct']:5.1f} {layer['silt_pct']:5.1f}"
            f"   -99   -99   -99   -99   -99   -99"
        )
    lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def process_soils_lucas(
    grid_points,              # GeoDataFrame
    output_dir_csv: str,
    output_dir_individual: str,
    n_cores: int,             # kept for API compat (nearest-point join is vectorised)
    id_col: str,
    lat_col: str,
    long_col: str,
    format_sql_func=None,     # kept for API compat; unused
    lucas_csv: str = "",
    max_dist_km: float = 50.0,
    col_map: Optional[dict] = None,
) -> bool:
    """Assign each grid point the nearest LUCAS topsoil sample and write a .SOL.

    Mirrors ``process_soils_ssurgo`` (signature + smart-resume + failure log),
    plus *lucas_csv* (the downloaded ESDAC LUCAS table) and *max_dist_km*.
    """
    print("Starting LUCAS Topsoil Processing (Smart Resume Mode)...")
    if not lucas_csv or not os.path.exists(lucas_csv):
        raise FileNotFoundError(
            "LUCAS needs lucas_csv pointing at the downloaded ESDAC LUCAS topsoil "
            "table (CSV/XLSX). Request it free at esdac.jrc.ec.europa.eu.")
    os.makedirs(output_dir_individual, exist_ok=True)
    lucas = _load_lucas(lucas_csv, col_map)
    print(f"LUCAS: loaded {len(lucas)} topsoil samples from {os.path.basename(lucas_csv)}")

    gdf = grid_points.copy()
    if hasattr(gdf, "geometry"):
        gdf = gdf.to_crs("EPSG:4326")
        gdf[lat_col] = gdf.geometry.y
        gdf[long_col] = gdf.geometry.x

    existing = {os.path.splitext(f)[0] for f in os.listdir(output_dir_individual) if f.endswith(".SOL")}
    todo = gdf[[str(p) not in existing for p in gdf[id_col].astype(str)]].reset_index(drop=True)
    print(f"Resume Check: Found {len(gdf) - len(todo)} existing profiles. Processing {len(todo)} remaining.")
    if todo.empty:
        print("All soil profiles already exist. Skipping LUCAS processing.")
        return True

    s_lat = lucas["lat"].to_numpy()
    s_lon = lucas["lon"].to_numpy()
    results, failures = [], []
    csv_header_written = os.path.exists(output_dir_csv)

    for _, prow in todo.iterrows():
        ID = str(prow[id_col]); lat = float(prow[lat_col]); lon = float(prow[long_col])
        if os.path.exists(os.path.join(output_dir_individual, f"{ID}.SOL")):
            continue
        dlat = np.radians(s_lat - lat); dlon = np.radians(s_lon - lon)
        a = np.sin(dlat / 2) ** 2 + math.cos(math.radians(lat)) * np.cos(np.radians(s_lat)) * np.sin(dlon / 2) ** 2
        dist = 6371.0 * 2 * np.arcsin(np.sqrt(a))
        k = int(np.argmin(dist))
        if dist[k] > max_dist_km:
            failures.append(_failure(ID, lat, lon,
                f"no-coverage: nearest LUCAS sample is {dist[k]:.0f} km away (> {max_dist_km:.0f} km; outside EU survey)"))
            continue
        rec = lucas.iloc[k]
        clay = float(rec["clay"]); sand = float(rec["sand"])
        silt = float(rec["silt"]) if np.isfinite(rec["silt"]) else max(0.0, 100 - clay - sand)
        om = float(rec["oc"]) / 10.0 * 1.724 if np.isfinite(rec["oc"]) else 1.0  # OC g/kg -> OC% -> OM%
        SLLL, SDUL, SSAT = _saxton_rawls(sand, clay, om)
        bd = max(0.9, min(1.8, (1.0 - SSAT) * 2.65))   # estimate BD from porosity (no measured BD in LUCAS)
        rows = []
        for top, bottom in [(0, 20), (20, _LUCAS_ROOTING_MAX_CM)]:
            rows.append({
                "ID": ID, "latitude": lat, "longitude": lon,
                "depth_top": top, "depth_bottom": bottom,
                "clay_pct": clay, "sand_pct": sand, "silt_pct": silt,
                "om_pct": om, "bulk_density": bd,
                "SLLL": SLLL, "SDUL": SDUL, "SSAT": SSAT,
            })
        profile_df = pd.DataFrame(rows)
        _write_sol(profile_df, output_dir_individual)
        results.append(profile_df)

    if results:
        pd.concat(results, ignore_index=True).to_csv(
            output_dir_csv, mode="a", index=False, header=not csv_header_written)
    if failures:
        fail_df = pd.DataFrame(failures)[["ID", "latitude", "longitude", "reason"]]
        base = os.path.splitext(os.path.basename(output_dir_csv))[0]
        failure_log = os.path.join(os.path.dirname(output_dir_csv), f"{base}_download_failures.csv")
        fail_df.to_csv(failure_log, index=False)
        print(f"[LUCAS] {len(fail_df)} of {len(todo)} point(s) had no LUCAS sample within "
              f"{max_dist_km:.0f} km. Details -> {failure_log}")
    print("LUCAS Topsoil Processing Complete.")
    return True
