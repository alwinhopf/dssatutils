# Internal helpers for local-raster soil sources.

import glob
import os

import numpy as np
import pandas as pd
from typing import Optional

from .soil_soilgrids_online import _calculate_soil_physics, _format_dssat_sol_file

_DEPTHS = [
    ("0_5", 5, 2.5),
    ("5_15", 15, 10.0),
    ("15_30", 30, 22.5),
    ("30_60", 60, 45.0),
    ("60_100", 100, 80.0),
    ("100_200", 200, 150.0),
]


def _find_raster(root: str, tokens, depth_token: Optional[str] = None) -> Optional[str]:
    tokens = [str(t).lower() for t in tokens]
    depth_bits = [] if depth_token is None else [depth_token.lower(), depth_token.replace("_", "-").lower()]
    for f in sorted(glob.glob(os.path.join(root, "**", "*"), recursive=True)):
        if not f.lower().endswith((".tif", ".tiff", ".vrt")):
            continue
        base = os.path.basename(f).lower()
        if all(t in base for t in tokens) and (not depth_bits or any(d in base for d in depth_bits)):
            return f
    return None


def _sample(path: Optional[str], lats, lons, scale: float = 1.0):
    if not path:
        return np.full(len(lats), np.nan, dtype=float)
    import rasterio
    from pyproj import Transformer

    out = np.full(len(lats), np.nan, dtype=float)
    with rasterio.open(path) as src:
        dst = src.crs.to_string() if src.crs else "EPSG:4326"
        xs, ys = Transformer.from_crs("EPSG:4326", dst, always_xy=True).transform(lons, lats)
        nodata = src.nodata
        for i, cell in enumerate(src.sample(zip(xs, ys), masked=True)):
            v = cell[0]
            if v is np.ma.masked or np.ma.is_masked(v):
                continue
            if nodata is not None and float(v) == float(nodata):
                continue
            out[i] = float(v) * scale
    return out


def _coords(grid_points, id_col, lat_col="LAT", long_col="LONG"):
    if hasattr(grid_points, "geometry") and grid_points.geometry.notna().any():
        g = grid_points.to_crs("EPSG:4326")
        lats = g.geometry.y.values.astype(float)
        lons = g.geometry.x.values.astype(float)
    else:
        g = pd.DataFrame(grid_points)
        lats = g[lat_col].astype(float).values
        lons = g[long_col].astype(float).values
    ids = g[id_col].astype(str).values
    return ids, lats, lons


def _texture_to_pct(cls):
    table = {
        1: (92, 5, 3), 2: (82, 12, 6), 3: (65, 25, 10), 4: (43, 39, 18),
        5: (20, 65, 15), 6: (10, 80, 10), 7: (52, 7, 41), 8: (45, 15, 40),
        9: (32, 34, 34), 10: (20, 40, 40), 11: (10, 34, 56), 12: (22, 20, 58),
    }
    if np.isnan(cls):
        return np.nan, np.nan, np.nan
    return table.get(int(round(cls)), (np.nan, np.nan, np.nan))


def write_mapping(ids, output_csv_path):
    pd.DataFrame({"ID": ids, "SOIL_ID": ids}).to_csv(output_csv_path, index=False)


def write_profiles(df, output_sol_dir, source_name, source_tag):
    os.makedirs(output_sol_dir, exist_ok=True)
    for _, site in df.groupby("ID"):
        _format_dssat_sol_file(site, output_sol_dir, source_name=source_name, source_tag=source_tag)


def add_physics_texture(df):
    phys = df.apply(lambda r: _calculate_soil_physics(r["sand"], r["clay"], r["soc_pct"] * 1.724), axis=1)
    df["SLLL"] = [p["SLLL"] for p in phys]
    df["SDUL"] = [p["SDUL"] for p in phys]
    df["SSAT"] = [p["SSAT"] for p in phys]
    return df


def normalize_bd(values):
    arr = np.asarray(values, dtype=float)
    # SoilGrids-style bdod often arrives in cg/cm3; SLGA/others may already be g/cm3.
    return np.where(arr > 10, arr / 100.0, arr)


DEPTHS = _DEPTHS
find_raster = _find_raster
sample_raster = _sample
coords = _coords
texture_to_pct = _texture_to_pct
