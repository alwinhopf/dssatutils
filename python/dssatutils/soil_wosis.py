# Soil source: processed WoSIS point profiles.

import os
import numpy as np
import pandas as pd

from .soil_raster_common import add_physics_texture, write_mapping, write_profiles


def _nearest(gp_lats, gp_lons, soil_lats, soil_lons):
    R = 6371000.0
    idx = []
    for la, lo in zip(gp_lats, gp_lons):
        dlat = np.radians(soil_lats - la)
        dlon = np.radians(soil_lons - lo)
        a = np.sin(dlat / 2) ** 2 + np.cos(np.radians(la)) * np.cos(np.radians(soil_lats)) * np.sin(dlon / 2) ** 2
        idx.append(int(np.argmin(2 * R * np.arcsin(np.sqrt(a)))))
    return np.array(idx)


def process_soils_wosis(
    grid_points,
    wosis_profile_csv: str,
    output_csv_path: str,
    output_sol_dir: str,
    id_col: str = "ID",
    lat_col: str = "LAT",
    long_col: str = "LONG",
) -> None:
    """Build DSSAT .SOL files from a processed WoSIS layer CSV.

    Required CSV columns: profile_id, latitude, longitude, depth_bottom, sand,
    clay, silt, bdod, soc_pct. Raw WoSIS exports should be harmonized to this
    schema before calling the backend.
    """
    if not os.path.exists(wosis_profile_csv):
        raise FileNotFoundError(wosis_profile_csv)
    src = pd.read_csv(wosis_profile_csv)
    required = ["profile_id", "latitude", "longitude", "depth_bottom", "sand", "clay", "silt", "bdod", "soc_pct"]
    missing = [c for c in required if c not in src.columns]
    if missing:
        raise ValueError(f"WoSIS processed CSV missing columns: {missing}")
    gp = grid_points.to_crs("EPSG:4326") if hasattr(grid_points, "geometry") else pd.DataFrame(grid_points)
    if hasattr(gp, "geometry"):
        gp_lats = gp.geometry.y.values.astype(float); gp_lons = gp.geometry.x.values.astype(float)
    else:
        gp_lats = gp[lat_col].astype(float).values; gp_lons = gp[long_col].astype(float).values
    ids = gp[id_col].astype(str).values
    heads = src.drop_duplicates("profile_id")
    ni = _nearest(gp_lats, gp_lons, heads["latitude"].values, heads["longitude"].values)
    rows = []
    for pid, src_id in zip(ids, heads.iloc[ni]["profile_id"].astype(str).values):
        sub = src[src["profile_id"].astype(str) == src_id].copy()
        sub["ID"] = pid
        sub["depth_center"] = sub.get("depth_center", sub["depth_bottom"] / 2)
        sub["cfvo"] = sub.get("cfvo", 0.0)
        rows.append(sub)
    out = add_physics_texture(pd.concat(rows, ignore_index=True))
    write_mapping(ids, output_csv_path)
    write_profiles(out, output_sol_dir, "WoSIS processed point profiles", "nearest processed profile")
