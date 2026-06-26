# Soil source: HiHydroSoil v2.0 global hydraulic-property rasters.

import os

import numpy as np
import pandas as pd

from .soil_raster_common import (
    DEPTHS, coords, find_raster, sample_raster, texture_to_pct, write_mapping,
    write_profiles,
)


def _first_existing(root, token_groups, depth):
    for tokens in token_groups:
        p = find_raster(root, tokens, depth)
        if p:
            return p
    return None


def process_soils_hihydrosoil(
    grid_points,
    hihydrosoil_raster_dir: str,
    output_csv_path: str,
    output_sol_dir: str,
    id_col: str = "ID",
    lat_col: str = "LAT",
    long_col: str = "LONG",
    integer_scale: float = 0.0001,
) -> None:
    """Build DSSAT .SOL files from local HiHydroSoil v2.0 GeoTIFF/VRT rasters.

    HiHydroSoil float layers are often distributed as integers scaled by 10000;
    *integer_scale* defaults to 0.0001 for hydraulic layers.
    """
    if not os.path.isdir(hihydrosoil_raster_dir):
        raise FileNotFoundError(f"HiHydroSoil raster directory not found: {hihydrosoil_raster_dir}")
    ids, lats, lons = coords(grid_points, id_col, lat_col, long_col)
    rows = []
    for dtoken, dbot, dctr in DEPTHS:
        slll = sample_raster(_first_existing(hihydrosoil_raster_dir, [["pf42"], ["pf4.2"], ["pF4.2"]], dtoken), lats, lons, integer_scale)
        sdul_path = _first_existing(hihydrosoil_raster_dir, [["pf25"], ["pf2.5"], ["pF2.5"], ["pf2"], ["pF2"]], dtoken)
        sdul = sample_raster(sdul_path, lats, lons, integer_scale)
        ssat = sample_raster(_first_existing(hihydrosoil_raster_dir, [["thetas"], ["theta_s"], ["sat"]], dtoken), lats, lons, integer_scale)
        ksat = sample_raster(_first_existing(hihydrosoil_raster_dir, [["ksat"], ["conductivity"]], dtoken), lats, lons, integer_scale)
        om = sample_raster(_first_existing(hihydrosoil_raster_dir, [["organic"], ["om"]], dtoken), lats, lons, integer_scale)
        tex = sample_raster(_first_existing(hihydrosoil_raster_dir, [["texture"], ["usda"]], dtoken), lats, lons, 1.0)
        sand = sample_raster(find_raster(hihydrosoil_raster_dir, ["sand"], dtoken), lats, lons)
        clay = sample_raster(find_raster(hihydrosoil_raster_dir, ["clay"], dtoken), lats, lons)
        silt = sample_raster(find_raster(hihydrosoil_raster_dir, ["silt"], dtoken), lats, lons)
        for i, pid in enumerate(ids):
            sa, si, cl = sand[i], silt[i], clay[i]
            if not np.isfinite(sa) or not np.isfinite(si) or not np.isfinite(cl):
                sa, si, cl = texture_to_pct(tex[i])
            bdod = (1.0 - ssat[i]) * 2.65 if np.isfinite(ssat[i]) else np.nan
            rows.append({
                "ID": pid, "latitude": lats[i], "longitude": lons[i],
                "depth_bottom": dbot, "depth_center": dctr,
                "sand": sa, "clay": cl, "silt": si, "bdod": bdod,
                "soc_pct": (om[i] / 1.724) if np.isfinite(om[i]) else 1.0,
                "cfvo": 0.0, "SLLL": slll[i], "SDUL": sdul[i],
                "SSAT": ssat[i],
                # ksat already scaled to cm/day; DSSAT SSKS is cm/h (/24).
                "SSKS": (ksat[i] / 24.0) if np.isfinite(ksat[i]) else np.nan,
            })
    df = pd.DataFrame(rows).dropna(subset=["sand", "clay", "bdod", "SLLL", "SDUL", "SSAT"])
    if df.empty:
        raise RuntimeError("No usable HiHydroSoil data extracted. Check raster names, scaling, and coordinates.")
    df["SSAT"] = np.maximum(df["SSAT"], df["SDUL"] + 0.04)
    df["SDUL"] = np.maximum(df["SDUL"], df["SLLL"] + 0.04)
    write_mapping(ids, output_csv_path)
    write_profiles(df, output_sol_dir, "HiHydroSoil v2.0", "local hydraulic rasters")
