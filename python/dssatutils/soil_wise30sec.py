# Soil source: WISE30sec derived soil properties (ISRIC), local rasters.
#
# WISE30sec (Batjes 2016) reports properties on 7 fixed depth layers, NOT the
# SoilGrids 0-5/5-15/... intervals. This backend reuses the generic per-property/
# per-depth raster reader but with WISE30sec's own depth table, so the texture/BD/
# OC -> Saxton-Rawls water-limit pipeline is shared while the depths stay correct.
#
# Expected input: one GeoTIFF/VRT per property per depth layer, with a property
# token (sand/clay/silt/bulk/organic carbon) AND a depth token (d1..d7) in the
# filename. WISE30sec is natively distributed as a map-unit raster plus an
# attribute table; rasterize that join to per-property/per-depth grids first.

from .soil_slga import process_soils_slga

# WISE30sec layers D1..D7: 0-20, 20-40, 40-60, 60-80, 80-100, 100-150, 150-200 cm.
_WISE_DEPTHS = [
    ("d1", 20, 10.0),
    ("d2", 40, 30.0),
    ("d3", 60, 50.0),
    ("d4", 80, 70.0),
    ("d5", 100, 90.0),
    ("d6", 150, 125.0),
    ("d7", 200, 175.0),
]


def process_soils_wise30sec(
    grid_points,
    wise30sec_raster_dir: str,
    output_csv_path: str,
    output_sol_dir: str,
    id_col: str = "ID",
    lat_col: str = "LAT",
    long_col: str = "LONG",
    depth_specs=None,
) -> None:
    """Build DSSAT .SOL files from local WISE30sec per-property/per-depth rasters.

    Filenames must carry a property token and a depth token (d1..d7). Pass
    *depth_specs* to override the default WISE30sec depth table.
    """
    return process_soils_slga(
        grid_points, wise30sec_raster_dir, output_csv_path, output_sol_dir,
        id_col=id_col, lat_col=lat_col, long_col=long_col,
        depth_specs=depth_specs or _WISE_DEPTHS,
        source_name="ISRIC WISE30sec", source_tag="WISE30sec local rasters")
