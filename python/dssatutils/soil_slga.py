# Soil source: Soil and Landscape Grid of Australia (SLGA).

import os

import numpy as np
import pandas as pd

from .soil_raster_common import (
    DEPTHS, add_physics_texture, coords, find_raster, normalize_bd,
    sample_raster, write_mapping, write_profiles,
)


def process_soils_slga(
    grid_points,
    slga_raster_dir: str,
    output_csv_path: str,
    output_sol_dir: str,
    id_col: str = "ID",
    lat_col: str = "LAT",
    long_col: str = "LONG",
    depth_specs=None,
    source_name: str = "Soil and Landscape Grid of Australia",
    source_tag: str = "SLGA local rasters",
) -> None:
    """Build DSSAT .SOL files from local per-property/per-depth GeoTIFF/VRT rasters.

    *depth_specs* is a list of ``(depth_token, depth_bottom_cm, depth_center_cm)``
    tuples; it defaults to the SoilGrids standard depths and can be overridden for
    products with different layering (e.g. WISE30sec). Rasters are matched by a
    property token (``sand``/``clay``/``silt``/``bulk``/``organic carbon``) AND the
    depth token appearing in the filename.
    """
    if not os.path.isdir(slga_raster_dir):
        raise FileNotFoundError(f"SLGA raster directory not found: {slga_raster_dir}")
    if depth_specs is None:
        depth_specs = DEPTHS
    ids, lats, lons = coords(grid_points, id_col, lat_col, long_col)
    rows = []
    for dtoken, dbot, dctr in depth_specs:
        sand = sample_raster(find_raster(slga_raster_dir, ["sand"], dtoken), lats, lons)
        clay = sample_raster(find_raster(slga_raster_dir, ["clay"], dtoken), lats, lons)
        silt_path = find_raster(slga_raster_dir, ["silt"], dtoken)
        silt = sample_raster(silt_path, lats, lons) if silt_path else 100.0 - sand - clay
        bd = normalize_bd(sample_raster(find_raster(slga_raster_dir, ["bulk"], dtoken)
                                        or find_raster(slga_raster_dir, ["bd"], dtoken), lats, lons))
        soc = sample_raster(find_raster(slga_raster_dir, ["organic", "carbon"], dtoken)
                            or find_raster(slga_raster_dir, ["soc"], dtoken), lats, lons)
        soc_pct = np.where(soc > 20, soc / 10.0, soc)
        for i, pid in enumerate(ids):
            rows.append({
                "ID": pid, "latitude": lats[i], "longitude": lons[i],
                "depth_bottom": dbot, "depth_center": dctr,
                "sand": sand[i], "clay": clay[i], "silt": silt[i],
                "bdod": bd[i], "soc_pct": soc_pct[i], "cfvo": 0.0,
            })
    df = pd.DataFrame(rows).dropna(subset=["sand", "clay", "bdod", "soc_pct"])
    if df.empty:
        raise RuntimeError("No usable SLGA soil data extracted. Check raster names and coordinates.")
    df = add_physics_texture(df)
    write_mapping(ids, output_csv_path)
    write_profiles(df, output_sol_dir, source_name, source_tag)
