# Soil source: China Soil Characteristics Dataset (BNU, Shangguan et al. 2013).
#
# National-grade 1 km soil for China (8 vertical layers), pairing with the CMFD
# weather backend. Reuses the generic per-property/per-depth raster reader with
# the 8-layer depth table. Expected input: one GeoTIFF/VRT per property per layer
# with a property token (sand/clay/silt/bulk/organic carbon) AND a layer token
# (l1..l8) in the filename.

from .soil_slga import process_soils_slga

# Same 8-layer scheme as GSDE (cm bottoms): 5, 9, 17, 29, 49, 83, 138, 230.
_CHINA_DEPTHS = [
    ("l1", 5, 2.0), ("l2", 9, 7.0), ("l3", 17, 13.0), ("l4", 29, 23.0),
    ("l5", 49, 39.0), ("l6", 83, 66.0), ("l7", 138, 110.0), ("l8", 230, 184.0),
]


def process_soils_china(
    grid_points, china_raster_dir: str, output_csv_path: str, output_sol_dir: str,
    id_col: str = "ID", lat_col: str = "LAT", long_col: str = "LONG",
    depth_specs=None,
) -> None:
    """Build DSSAT .SOL files from local China (BNU) per-property/per-layer rasters."""
    return process_soils_slga(
        grid_points, china_raster_dir, output_csv_path, output_sol_dir,
        id_col=id_col, lat_col=lat_col, long_col=long_col,
        depth_specs=depth_specs or _CHINA_DEPTHS,
        source_name="China Soil Dataset (BNU)", source_tag="BNU local rasters")
