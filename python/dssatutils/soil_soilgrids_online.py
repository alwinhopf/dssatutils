# File: soil_soilgrids_online.py
# Python port of soil_soilgrids_online.R
#
# Fetches SoilGrids 2.0 data via the ISRIC REST API or VRT/GDAL virtual
# rasters, computes DSSAT soil physics (Saxton & Rawls 2006), and writes
# individual per-point .SOL files plus a mapping CSV.
#
# REST API docs: https://rest.isric.org/soilgrids/v2.0/
# VRT data root: https://files.isric.org/soilgrids/latest/data/

import os
import re
import time
import warnings
from datetime import datetime
from typing import Optional

import numpy as np
import pandas as pd
import requests

# USE_REST_API = True   → JSON REST API (interactive / local; rate-limited)
# USE_REST_API = False  → VRT/GDAL virtual rasters (HPC / batch)
# Can be overridden from dssat_main_pipeline.py before calling the function.
USE_REST_API: bool = False

_ISRIC_REST_URL = "https://rest.isric.org/soilgrids/v2.0/properties/query"
_ISRIC_VRT_ROOT = "https://files.isric.org/soilgrids/latest/data/"

_PROPS  = ["clay", "sand", "silt", "soc", "bdod", "cfvo"]
_DEPTHS = ["0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm"]
_DEPTH_CENTERS = [2.5, 10.0, 22.5, 45.0, 80.0, 150.0]
_DEPTH_BOTTOMS = [5,   15,   30,   60,   100,   200]


# ---------------------------------------------------------------------------
# 1. Physics (Saxton & Rawls 2006)
# ---------------------------------------------------------------------------

def _calculate_soil_physics(sand_pct: float, clay_pct: float,
                             om_pct: float) -> dict[str, float]:
    """
    Return SLLL (wilting point), SDUL (field capacity), SSAT (saturation)
    as volumetric fractions.  Mirrors ``calculate_soil_physics`` in R.
    """
    S  = sand_pct / 100.0
    C  = clay_pct / 100.0
    OM = om_pct   / 100.0

    theta_1500t = (-0.024 * S + 0.487 * C + 0.006 * OM
                   + 0.005 * S * OM - 0.013 * C * OM
                   + 0.068 * S * C + 0.031)
    SLLL = theta_1500t + (0.14 * theta_1500t - 0.02)

    theta_33t = (-0.251 * S + 0.195 * C + 0.011 * OM
                 + 0.006 * S * OM - 0.027 * C * OM
                 + 0.452 * S * C + 0.299)
    SDUL = theta_33t + (1.283 * theta_33t**2 - 0.374 * theta_33t - 0.015)

    theta_s33t = (0.278 * S + 0.034 * C + 0.022 * OM
                  - 0.018 * S * OM - 0.027 * C * OM
                  - 0.584 * S * C + 0.078)
    theta_s33 = theta_s33t + (0.636 * theta_s33t - 0.107)
    SSAT = SDUL + theta_s33 - 0.097 * S + 0.043

    return {"SLLL": float(SLLL), "SDUL": float(SDUL), "SSAT": float(SSAT)}


# ---------------------------------------------------------------------------
# 2. DSSAT .SOL formatter
# ---------------------------------------------------------------------------

def _format_dssat_sol_file(site_data: pd.DataFrame, output_dir: str,
                           source_tag: str = "REST API",
                           source_name: str = "ISRIC SoilGrids 2.0") -> None:
    """Write one DSSAT .SOL file.  Mirrors ``format_dssat_sol_file`` in R.

    *source_name* is the product family (e.g. "ISRIC SoilGrids 2.0",
    "FAO HWSD v2.0") and *source_tag* the access detail (e.g. "REST API",
    "VRT", "HWSD2"); both are recorded in the header so provenance is visible.
    """
    if site_data.empty:
        raise ValueError("No soil layers found for this ID.")

    critical = ["sand", "clay", "bdod", "soc_pct"]
    if site_data[critical].isna().any().any():
        raise ValueError(
            "Critical soil physical data (Sand, Clay, Bulk Density, or OC) contains NAs."
        )

    soil_id = str(site_data["ID"].iloc[0])
    lat     = float(site_data["latitude"].iloc[0])
    lon     = float(site_data["longitude"].iloc[0])
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, f"{soil_id}.SOL")

    lines = [
        f"*SOILS: {source_name}",
        f"! Source: {source_name} ({source_tag})",
        "",
        f"*{soil_id[:10]:<10s}  SOILGRIDS     {lat:9.3f} {lon:9.3f}",
        "@SITE        COUNTRY          LAT     LONG SCS FAMILY",
        f" {soil_id[:11]:<11s} World         {lat:9.3f} {lon:9.3f} ",
        "@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE",
        "    BN   .13     6    .6    73     1     1 IB001 IB001 IB001",
        "@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC",
    ]

    site_data = site_data.sort_values("depth_bottom").reset_index(drop=True)

    for _, layer in site_data.iterrows():
        dc = float(layer["depth_center"])
        srgf = max(0.0, 1.0 * np.exp(-0.02 * dc))
        if srgf < 0.02:
            srgf = 0.0

        sand = float(layer["sand"])
        clay = float(layer["clay"])
        ssks = min(999.0, 60.96 * (10 ** (0.0126 * sand - 0.0064 * clay - 0.6)))

        lines.append(
            f"{int(layer['depth_bottom']):6d}   -99"
            f" {layer['SLLL']:5.3f} {layer['SDUL']:5.3f} {layer['SSAT']:5.3f}"
            f" {srgf:5.2f} {ssks:5.1f} {layer['bdod']:5.2f} {layer['soc_pct']:5.2f}"
            f" {layer['clay']:5.1f} {layer['silt']:5.1f} {layer['cfvo']:5.1f}"
            f"   -99   -99   -99   -99   -99"
        )
    lines.append("")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# 3. REST API fetcher (with exponential back-off)
# ---------------------------------------------------------------------------

def _fetch_soilgrids_rest(lat: float, lon: float,
                           max_retries: int = 5, base_wait: float = 2.0
                           ) -> Optional[pd.DataFrame]:
    """
    Fetch SoilGrids mean values for one point via the REST API.
    Returns a tidy DataFrame: prop, depth_label, depth_bottom, depth_center, value.
    Mirrors ``fetch_soilgrids_rest`` in R.
    """
    params: list[tuple] = [
        ("lat", lat), ("lon", lon), ("value", "mean"),
    ]
    for p in _PROPS:
        params.append(("property", p))
    for d in _DEPTHS:
        params.append(("depth", d))

    for attempt in range(max_retries):
        try:
            time.sleep(1)  # polite delay
            r = requests.get(_ISRIC_REST_URL, params=params, timeout=90)
            status = r.status_code

            if status == 200:
                data = r.json()
                layers = data.get("properties", {}).get("layers", [])
                if not layers:
                    return None
                records = []
                for layer in layers:
                    prop_name = layer["name"]
                    for depth_info in layer["depths"]:
                        label = depth_info["label"]
                        val   = depth_info["values"].get("mean")
                        nums  = [int(x) for x in re.findall(r"\d+", label)]
                        records.append({
                            "prop": prop_name,
                            "depth_label": label,
                            "depth_bottom": nums[1] if len(nums) >= 2 else None,
                            "depth_center": (nums[0] + nums[1]) / 2 if len(nums) >= 2 else None,
                            "value": val,
                        })
                return pd.DataFrame(records)

            elif status in (429, 503, 504):
                wait = base_wait * (2 ** attempt)
                print(f"Rate limit (HTTP {status}). Retrying in {wait:.0f}s…")
                time.sleep(wait)
            else:
                warnings.warn(f"SoilGrids REST fatal error: HTTP {status}")
                return None

        except Exception as exc:
            warnings.warn(f"REST fetch failed: {exc}")
            time.sleep(2)

    return None


# ---------------------------------------------------------------------------
# 4. VRT fetcher
# ---------------------------------------------------------------------------

def _fetch_soilgrids_vrt(gridfile, id_col: str) -> pd.DataFrame:
    """
    Extract SoilGrids values for all points using GDAL VRT virtual rasters.
    Requires ``rasterio`` (with GDAL >= 3.x) or ``rioxarray``.
    Mirrors ``fetch_soilgrids_vrt`` in R.
    """
    try:
        import rasterio
        from pyproj import Transformer
    except ImportError as exc:
        raise ImportError(
            "rasterio and pyproj are required for VRT mode. "
            "Install with: pip install rasterio pyproj"
        ) from exc

    # Project points to IGH (Interrupted Goode's Homolosine) — SoilGrids native CRS
    IGH_PROJ = ("+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs")

    import geopandas as gpd
    gdf = gridfile.to_crs("EPSG:4326")
    lons = gdf.geometry.x.values
    lats = gdf.geometry.y.values

    transformer = Transformer.from_crs("EPSG:4326", IGH_PROJ, always_xy=True)
    xs, ys = transformer.transform(lons, lats)

    records = []
    for prop in _PROPS:
        print(f"VRT: Extracting {prop}...")
        for i, depth in enumerate(_DEPTHS):
            vrt_url = (
                f"/vsicurl/{_ISRIC_VRT_ROOT}"
                f"{prop}/{prop}_{depth}_mean.vrt"
            )
            try:
                with rasterio.open(vrt_url) as src:
                    # Sample at each point. masked=True returns a numpy masked
                    # array so SoilGrids nodata (-32768, e.g. over ocean) becomes
                    # masked rather than a raw sentinel. We convert masked/nodata
                    # entries to NaN so such points are skipped downstream instead
                    # of producing garbage soil values. (terra::extract in the R
                    # version masks nodata automatically; rasterio.sample does not
                    # unless masked=True is passed.)
                    vals = list(src.sample(zip(xs, ys), masked=True))
                    vals_flat = []
                    for v in vals:
                        cell = v[0] if v is not None else None
                        if cell is None or cell is np.ma.masked or \
                                (np.ma.is_masked(cell)):
                            vals_flat.append(np.nan)
                        else:
                            vals_flat.append(float(cell))

                for j, pid in enumerate(gridfile[id_col]):
                    records.append({
                        "ID": pid,
                        "prop": prop,
                        "depth_label": depth,
                        "depth_bottom": _DEPTH_BOTTOMS[i],
                        "depth_center": _DEPTH_CENTERS[i],
                        "value": vals_flat[j],
                        "VRT": True,
                    })
            except Exception as exc:
                print(f"  Skip {prop} {depth}: {exc}")

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# 5. Main pipeline function
# ---------------------------------------------------------------------------

def process_soils_soilgrids_online(
    gridfile,               # GeoDataFrame
    soilfile_csv_path: str,
    output_sol_dir: str,
    id_col: str,
) -> None:
    """
    Fetch SoilGrids data, compute DSSAT soil physics, write per-point .SOL
    files and a mapping CSV.  Mirrors ``process_soils_soilgrids_online`` in R.

    The global variable ``USE_REST_API`` (default False) switches between:
      True  → JSON REST API (one request per point; rate-limited)
      False → VRT/GDAL virtual rasters (batch-friendly; requires GDAL)
    """
    use_rest = USE_REST_API
    mode = "REST API" if use_rest else "VRT"
    print(f"--- Starting SoilGrids Extraction (Mode: {mode}) ---")

    if not hasattr(gridfile, "geometry"):
        raise TypeError("gridfile must be a GeoDataFrame (with .geometry column).")

    grid_wgs84 = gridfile.to_crs("EPSG:4326")
    grid_wgs84 = grid_wgs84.copy()
    grid_wgs84["lon_wgs84"] = grid_wgs84.geometry.x
    grid_wgs84["lat_wgs84"] = grid_wgs84.geometry.y

    full_df: Optional[pd.DataFrame] = None

    # --- Branch: REST API ---
    if use_rest:
        results = []
        n = len(grid_wgs84)
        for i, (_, row) in enumerate(grid_wgs84.iterrows()):
            if i % 10 == 0:
                print(f"  Processed {i} / {n} points...")
            pid = row[id_col]
            res = _fetch_soilgrids_rest(float(row["lat_wgs84"]), float(row["lon_wgs84"]))
            if res is not None:
                res["ID"] = pid
                results.append(res)
        if results:
            full_df = pd.concat(results, ignore_index=True)

    # --- Branch: VRT ---
    else:
        full_df = _fetch_soilgrids_vrt(gridfile, id_col)

    if full_df is None or full_df.empty:
        raise RuntimeError(
            "No soil data extracted! Check internet connection or coordinates."
        )

    # -----------------------------------------------------------------------
    # Restructure to wide format and compute physics
    # -----------------------------------------------------------------------
    print("Restructuring and calculating soil physics...")

    wide_df = full_df.pivot_table(
        index=["ID", "depth_label", "depth_bottom", "depth_center"],
        columns="prop",
        values="value",
        aggfunc="first",
    ).reset_index()

    # pivot_table silently DROPS a property column when every value for it is
    # null, and produces ZERO rows when a point returns all-null values (e.g. a
    # point over water or a SoilGrids coverage gap). Guarantee all expected
    # property columns exist (filled with NaN) so the unit conversions below
    # never raise KeyError, and bail out cleanly if nothing usable came back.
    for _prop in _PROPS:
        if _prop not in wide_df.columns:
            wide_df[_prop] = np.nan

    _input_ids = set(full_df["ID"].unique())

    if not wide_df.empty:
        # Drop points whose sand/clay are entirely NaN — no usable physics.
        _usable = wide_df.groupby("ID")[["sand", "clay"]].transform(
            lambda s: s.notna().any()
        ).any(axis=1)
        wide_df = wide_df[_usable].reset_index(drop=True)

    _kept_ids    = set(wide_df["ID"].unique()) if not wide_df.empty else set()
    _skipped_ids = sorted(_input_ids - _kept_ids)
    if _skipped_ids:
        print(f"  [warn] {len(_skipped_ids)} point(s) had no SoilGrids data "
              f"and were skipped: {', '.join(map(str, _skipped_ids[:10]))}"
              + (" ..." if len(_skipped_ids) > 10 else ""))

    if wide_df.empty:
        raise RuntimeError(
            "No usable soil data for any point (all coordinates returned null). "
            "Check coordinates or try SOILGRIDS_10K / SSURGO."
        )

    # Unit conversions (matching R script)
    wide_df["clay"]    = wide_df["clay"]  / 10.0    # g/kg → %
    wide_df["sand"]    = wide_df["sand"]  / 10.0
    wide_df["silt"]    = wide_df["silt"]  / 10.0
    wide_df["soc_pct"] = wide_df["soc"]   / 100.0   # dg/kg → %
    wide_df["om_pct"]  = wide_df["soc_pct"] * 1.724
    wide_df["bdod"]    = wide_df["bdod"]  / 100.0   # cg/cm³ → g/cm³
    wide_df["cfvo"]    = wide_df["cfvo"]  / 10.0    # cm³/dm³ → %

    # Calculate physics row-by-row
    physics_rows = wide_df.apply(
        lambda r: pd.Series(
            _calculate_soil_physics(
                float(r["sand"]) if not pd.isna(r["sand"]) else 40.0,
                float(r["clay"]) if not pd.isna(r["clay"]) else 20.0,
                float(r["om_pct"]) if not pd.isna(r["om_pct"]) else 1.0,
            )
        ),
        axis=1,
    )
    processed_df = pd.concat([wide_df, physics_rows], axis=1)

    # Join coordinates back
    coords_df = grid_wgs84[[id_col, "lon_wgs84", "lat_wgs84"]].rename(
        columns={"lon_wgs84": "longitude", "lat_wgs84": "latitude"}
    )
    final_df = processed_df.merge(coords_df, on=id_col, how="inner")
    final_df = final_df.rename(columns={id_col: "ID"})

    # -----------------------------------------------------------------------
    # Write outputs
    # -----------------------------------------------------------------------
    # 1. Mapping CSV
    mapping_df = pd.DataFrame({
        id_col: grid_wgs84[id_col].values,
        "SOIL_ID": grid_wgs84[id_col].values,
    })
    mapping_df.to_csv(soilfile_csv_path, index=False)

    # 2. Individual SOL files with error logging
    os.makedirs(output_sol_dir, exist_ok=True)
    log_path = os.path.join(output_sol_dir, "soil_processing_errors.log")
    with open(log_path, "w") as lf:
        lf.write(f"Log started: {datetime.now()}\n")

    unique_ids = final_df["ID"].unique()
    print(f"Writing {len(unique_ids)} potential .SOL files to: {output_sol_dir}")

    success = errors = 0
    for uid in unique_ids:
        subset = final_df[final_df["ID"] == uid].copy()
        try:
            _format_dssat_sol_file(subset, output_sol_dir, source_tag=mode)
            success += 1
        except Exception as exc:
            msg = f"ID: {uid} | Error: {exc}"
            print(f"SKIPPED: {msg}")
            with open(log_path, "a") as lf:
                lf.write(msg + "\n")
            errors += 1

    print(f"SoilGrids processing complete. Success: {success}, Errors: {errors}")
    print(f"Check {log_path} for details on skipped points.")
