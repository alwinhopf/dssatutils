# File: soil_isdasoil.py
# Python port of soil_isdasoil.R
#
# Soil source: iSDAsoil (Innovative Solutions for Decision Agriculture) ->
# DSSAT .SOL. 30 m predicted soil properties for AFRICA at two depths
# (0-20 cm, 20-50 cm).
#
# WHY: iSDAsoil is the highest-resolution open soil product for Africa (30 m vs
# SoilGrids' 250 m), so for African sites it resolves field-scale soil variation
# that the global products miss.
#
# ACCESS (fully open, NO key): cloud-optimised GeoTIFFs in a public S3 bucket,
# streamed per-point via GDAL /vsicurl. Each property file has 4 bands:
#   band 1 = mean 0-20 cm, band 2 = mean 20-50 cm, bands 3-4 = std-dev.
# Native CRS EPSG:3857 (Web Mercator), nodata = 255.
#   https://isdasoil.s3.amazonaws.com/soil_data/<property>/<property>.tif
#
# Stored uint8 values are back-transformed to real units (verified against the
# iSDA / Google Earth Engine catalogue):
#   clay_content / sand_content : value as-is (%)
#   carbon_organic              : exp(x/10) - 1   (g/kg)  -> OC% -> OM%
#   bulk_density                : x / 100          (g/cm^3)
#   ph                          : x / 10           (-)
#
# DSSAT physics (SLLL/SDUL/SSAT via Saxton & Rawls) and the .SOL layout are the
# same as soil_ssurgo.py, so an iSDAsoil profile is comparable to a SSURGO one.
# Coverage: Africa.

import math
import os
import warnings
from typing import Optional

import numpy as np
import pandas as pd

from .soil_ssurgo import _failure, _saxton_rawls

_ISDA_BASE = "https://isdasoil.s3.amazonaws.com/soil_data"
_ISDA_NODATA = 255

# DSSAT depth layers built from the two iSDAsoil depths. The 20-50 cm prediction
# is also carried down to ROOTING_MAX_CM so DSSAT has a usable rooting profile
# (iSDAsoil itself only predicts to 50 cm).
_ROOTING_MAX_CM = 150
_LAYERS = [
    # (depth_bottom_cm, iSDA depth-band index: 0 = 0-20cm, 1 = 20-50cm)
    (20, 0),
    (50, 1),
    (_ROOTING_MAX_CM, 1),
]


def _back_transform(prop: str, raw: np.ndarray) -> np.ndarray:
    """Convert stored uint8 iSDAsoil values to real units."""
    v = raw.astype(float)
    v[raw == _ISDA_NODATA] = np.nan
    if prop == "carbon_organic":
        return np.expm1(v / 10.0)        # g/kg
    if prop == "bulk_density":
        return v / 100.0                  # g/cm^3
    if prop == "ph":
        return v / 10.0
    return v                              # clay_content / sand_content: %


def _sample_property(prop: str, xs, ys) -> Optional[np.ndarray]:
    """Sample the mean bands (0-20, 20-50 cm) of one iSDAsoil property.

    *xs*, *ys* are coordinates already in the raster CRS (EPSG:3857). Returns an
    (n_points, 2) array of real-unit values, or None if the COG can't be read.
    """
    import rasterio
    url = f"/vsicurl/{_ISDA_BASE}/{prop}/{prop}.tif"
    try:
        with rasterio.open(url) as ds:
            samples = np.array([list(v)[:2] for v in ds.sample(list(zip(xs, ys)))])
    except Exception as exc:  # noqa: BLE001
        warnings.warn(f"iSDAsoil: could not read {prop} COG: {exc}")
        return None
    out = np.column_stack([
        _back_transform(prop, samples[:, 0]),
        _back_transform(prop, samples[:, 1]),
    ])
    return out


# ---------------------------------------------------------------------------
# DSSAT .SOL writer (iSDAsoil-labelled; same column layout as SSURGO)
# ---------------------------------------------------------------------------

def _write_sol(profile: pd.DataFrame, output_dir: str) -> None:
    soil_id = str(profile["ID"].iloc[0])
    lat = profile["latitude"].iloc[0]
    lon = profile["longitude"].iloc[0]
    path = os.path.join(output_dir, f"{soil_id}.SOL")
    if os.path.exists(path):
        return
    lines = [
        "*SOILS: Africa iSDAsoil Soil Profiles",
        "! Generated from iSDAsoil 30 m (0-20 & 20-50 cm), Saxton & Rawls physics",
        "",
        f"*{soil_id:<6s}  ISDA          {lat:9.3f} {lon:9.3f}",
        "@SITE        COUNTRY          LAT     LONG SCS FAMILY",
        f" {soil_id:<11s} -99         {lat:9.3f} {lon:9.3f} ",
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


# ---------------------------------------------------------------------------
# Public entry point — signature mirrors process_soils_ssurgo
# ---------------------------------------------------------------------------

def process_soils_isdasoil(
    grid_points,              # GeoDataFrame
    output_dir_csv: str,
    output_dir_individual: str,
    n_cores: int,             # kept for API compat (COG sampling is one batched pass)
    id_col: str,
    lat_col: str,
    long_col: str,
    format_sql_func=None,     # kept for API compat; unused
) -> bool:
    """Sample iSDAsoil for every point, write per-point .SOL + a mapping CSV.

    Mirrors ``process_soils_ssurgo`` (same signature, smart-resume, failure log).
    All points are sampled in a few batched COG passes (one per property), which
    is far faster than per-point HTTP for SoilGrids-style sources.
    """
    import geopandas as gpd
    from pyproj import Transformer

    print("Starting iSDAsoil Processing (Smart Resume Mode)...")
    os.makedirs(output_dir_individual, exist_ok=True)

    gdf = grid_points.copy()
    if hasattr(gdf, "geometry"):
        gdf = gdf.to_crs("EPSG:4326")
        gdf[lat_col] = gdf.geometry.y
        gdf[long_col] = gdf.geometry.x

    existing = {os.path.splitext(f)[0] for f in os.listdir(output_dir_individual) if f.endswith(".SOL")}
    all_ids = gdf[id_col].astype(str).tolist()
    keep = [str(pid) not in existing for pid in all_ids]
    todo = gdf[keep].reset_index(drop=True)
    print(f"Resume Check: Found {len(all_ids) - len(todo)} existing profiles. Processing {len(todo)} remaining.")
    if todo.empty:
        print("All soil profiles already exist. Skipping iSDAsoil processing.")
        return True

    ids = todo[id_col].astype(str).tolist()
    lats = todo[lat_col].astype(float).to_numpy()
    lons = todo[long_col].astype(float).to_numpy()
    xs, ys = Transformer.from_crs(4326, 3857, always_xy=True).transform(lons, lats)

    props = {}
    for prop in ("clay_content", "sand_content", "carbon_organic", "bulk_density"):
        props[prop] = _sample_property(prop, np.atleast_1d(xs), np.atleast_1d(ys))
    if props["clay_content"] is None or props["sand_content"] is None:
        print("iSDAsoil: clay/sand COGs unreadable; aborting (check network / S3 access).")
        return False

    results, failures = [], []
    csv_header_written = os.path.exists(output_dir_csv)

    for i, (ID, lat, lon) in enumerate(zip(ids, lats, lons)):
        if os.path.exists(os.path.join(output_dir_individual, f"{ID}.SOL")):
            continue
        # Surface (0-20) clay decides coverage: nodata => outside Africa / water.
        if not np.isfinite(props["clay_content"][i, 0]):
            failures.append(_failure(ID, lat, lon,
                "no-coverage: no iSDAsoil value at this location (outside Africa / water / nodata)"))
            continue

        layer_rows = []
        for (bottom, depth_idx) in _LAYERS:
            top = 0 if bottom == 20 else (20 if bottom == 50 else 50)
            clay = props["clay_content"][i, depth_idx]
            sand = props["sand_content"][i, depth_idx]
            oc_gkg = props["carbon_organic"][i, depth_idx] if props["carbon_organic"] is not None else np.nan
            bd = props["bulk_density"][i, depth_idx] if props["bulk_density"] is not None else np.nan
            clay = float(clay) if np.isfinite(clay) else 20.0
            sand = float(sand) if np.isfinite(sand) else 40.0
            om = float(oc_gkg) / 10.0 * 1.724 if np.isfinite(oc_gkg) else 1.0  # g/kg->OC%->OM%
            bd = float(bd) if np.isfinite(bd) else 1.4
            silt = max(0.0, 100.0 - clay - sand)
            SLLL, SDUL, SSAT = _saxton_rawls(sand, clay, om)  # Saxton & Rawls takes OM%
            layer_rows.append({
                "ID": ID, "latitude": lat, "longitude": lon,
                "depth_top": top, "depth_bottom": bottom,
                "clay_pct": clay, "sand_pct": sand, "silt_pct": silt,
                "om_pct": om, "bulk_density": bd,
                "SLLL": SLLL, "SDUL": SDUL, "SSAT": SSAT,
            })
        profile_df = pd.DataFrame(layer_rows)
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
        print(f"[iSDAsoil] {len(fail_df)} of {len(ids)} point(s) produced NO soil profile "
              f"(no-coverage: outside Africa / water). Details -> {failure_log}")

    print("iSDAsoil Processing Complete.")
    return True
