# Soil source: OpenLandMap / OpenGeoHub global soil COGs (live remote sampling).
#
# Samples cloud-optimized GeoTIFFs over HTTP (no bulk download) via the OpenLandMap
# STAC catalog, so URLs are resolved at runtime and survive version bumps. Texture
# (clay/sand/silt) comes from the global_soil_props product (120 m, 0-30/30-60/
# 60-100 cm); bulk density and organic carbon from the 250 m property COGs at the
# point-depth nearest each layer centre. Water limits are derived via Saxton-Rawls.
#
# Catalog: https://s3.eu-central-1.wasabisys.com/stac/openlandmap/catalog.json

import re

import numpy as np
import pandas as pd
import requests

from .soil_raster_common import (
    add_physics_texture, coords, normalize_bd, write_mapping, write_profiles,
)

_STAC_BASE = "https://s3.eu-central-1.wasabisys.com/stac/openlandmap/"
# property -> STAC collection id
_COLLECTIONS = {
    "clay": "clay.tot_iso.11277.2020.wpct",
    "sand": "sand.tot_iso.11277.2020.wpct",
    "silt": "silt.tot_iso.11277.2020.wpct",
    "bdod": "bulkdens.fineearth_usda.4a1h",
    "oc":   "organic.carbon_usda.6a1c",
}
# DSSAT layers (depth_bottom_cm, depth_center_cm) matching the texture product.
_LAYERS = [(30, 15.0), (60, 45.0), (100, 80.0)]


def _depth_mid(key: str):
    """Parse a depth midpoint (cm) from a STAC asset key like 'b0cm..30cm' or 'b10cm'."""
    m = re.search(r"b(\d+)cm\.\.(\d+)cm", key)
    if m:
        return (int(m.group(1)) + int(m.group(2))) / 2.0
    m = re.search(r"b(\d+)cm", key)
    return float(m.group(1)) if m else None


def _resolve_assets(collection: str):
    """Return [(depth_mid_cm, cog_url)] for the mean COG assets of a collection."""
    col = requests.get(f"{_STAC_BASE}{collection}/collection.json", timeout=60).json()
    item_link = next(l for l in col["links"] if l.get("rel") == "item")
    item_url = requests.compat.urljoin(f"{_STAC_BASE}{collection}/", item_link["href"])
    assets = requests.get(item_url, timeout=60).json()["assets"]
    out = []
    for k, v in assets.items():
        href = v.get("href", "")
        # mean COGs only (key has '_m_'), prefer the finer 120 m where present.
        if "_m_" in k and href.endswith(".tif") and "preview" not in k:
            mid = _depth_mid(k)
            if mid is not None:
                out.append((mid, href, "120m" in k))
    # prefer 120 m assets when both resolutions exist for the same depth
    best = {}
    for mid, href, fine in out:
        if mid not in best or (fine and not best[mid][1]):
            best[mid] = (href, fine)
    return sorted((mid, href) for mid, (href, _) in best.items())


def _sample_cog(url: str, lats, lons):
    import rasterio
    from pyproj import Transformer
    out = np.full(len(lats), np.nan, dtype=float)
    with rasterio.open("/vsicurl/" + url) as src:
        dst = src.crs.to_string() if src.crs else "EPSG:4326"
        xs, ys = Transformer.from_crs("EPSG:4326", dst, always_xy=True).transform(lons, lats)
        nd = src.nodata
        for i, cell in enumerate(src.sample(zip(xs, ys))):
            v = float(cell[0])
            if nd is not None and v == float(nd):
                continue
            out[i] = v
    return out


def _nearest_asset(assets, target_mid):
    return min(assets, key=lambda a: abs(a[0] - target_mid))[1]


def process_soils_openlandmap(
    grid_points, output_csv_path: str, output_sol_dir: str,
    id_col: str = "ID", lat_col: str = "LAT", long_col: str = "LONG",
) -> None:
    """Build DSSAT .SOL files by sampling OpenLandMap COGs over HTTP (live).

    No local data required; resolves COG URLs from the OpenLandMap STAC catalog
    and samples each point. Bulk density is scaled to g/cm3 and organic carbon
    (stored as x5 g/kg) to percent before the Saxton-Rawls water-limit step.
    """
    ids, lats, lons = coords(grid_points, id_col, lat_col, long_col)
    print("--- Starting OpenLandMap COG Extraction (live) ---")
    assets = {p: _resolve_assets(c) for p, c in _COLLECTIONS.items()}

    rows = []
    for dbot, dctr in _LAYERS:
        # texture intervals align with layer; BD/OC use nearest point-depth.
        clay = _sample_cog(_nearest_asset(assets["clay"], dctr), lats, lons)
        sand = _sample_cog(_nearest_asset(assets["sand"], dctr), lats, lons)
        silt = _sample_cog(_nearest_asset(assets["silt"], dctr), lats, lons)
        bd = normalize_bd(_sample_cog(_nearest_asset(assets["bdod"], dctr), lats, lons))
        oc = _sample_cog(_nearest_asset(assets["oc"], dctr), lats, lons)
        soc_pct = oc / 50.0  # x5 g/kg -> g/kg (/5) -> percent (/10)
        for i, pid in enumerate(ids):
            sa, si, cl = sand[i], silt[i], clay[i]
            if np.isfinite(sa) and np.isfinite(cl) and not np.isfinite(si):
                si = max(0.0, 100.0 - sa - cl)
            rows.append({
                "ID": pid, "latitude": lats[i], "longitude": lons[i],
                "depth_bottom": dbot, "depth_center": dctr,
                "sand": sa, "clay": cl, "silt": si, "bdod": bd[i],
                "soc_pct": soc_pct[i], "cfvo": 0.0,
            })
    df = pd.DataFrame(rows).dropna(subset=["sand", "clay", "bdod", "soc_pct"])
    if df.empty:
        raise RuntimeError("No usable OpenLandMap data sampled. Check connectivity / coordinates.")
    df = add_physics_texture(df)
    write_mapping(ids, output_csv_path)
    write_profiles(df, output_sol_dir, "OpenLandMap (OpenGeoHub)", "STAC COG sampling")
