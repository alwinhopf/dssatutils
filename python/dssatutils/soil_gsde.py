# Soil source: GSDE (Global Soil Dataset for Earth System Modeling), local rasters.
#
# GSDE (Shangguan et al. 2014) is a global 1 km product with 8 vertical layers to
# ~2.3 m. Reuses the generic per-property/per-depth raster reader with GSDE's own
# depth table so the texture/BD/OC -> Saxton-Rawls water-limit pipeline is shared.
# Expected input: one GeoTIFF/VRT per property per layer with a property token
# (sand/clay/silt/bulk/organic carbon) AND a layer token (l1..l8) in the filename.

from .soil_slga import process_soils_slga

# GSDE 8 layers (cm bottoms): 4.5, 9.1, 16.6, 28.9, 49.3, 82.9, 138.3, 229.6.
_GSDE_DEPTHS = [
    ("l1", 5, 2.0), ("l2", 9, 7.0), ("l3", 17, 13.0), ("l4", 29, 23.0),
    ("l5", 49, 39.0), ("l6", 83, 66.0), ("l7", 138, 110.0), ("l8", 230, 184.0),
]


def process_soils_gsde(
    grid_points, gsde_raster_dir: str, output_csv_path: str, output_sol_dir: str,
    id_col: str = "ID", lat_col: str = "LAT", long_col: str = "LONG",
    depth_specs=None,
) -> None:
    """Build DSSAT .SOL files from local GSDE per-property/per-layer rasters."""
    return process_soils_slga(
        grid_points, gsde_raster_dir, output_csv_path, output_sol_dir,
        id_col=id_col, lat_col=lat_col, long_col=long_col,
        depth_specs=depth_specs or _GSDE_DEPTHS,
        source_name="GSDE (Global Soil Dataset for ESM)", source_tag="GSDE local rasters")
