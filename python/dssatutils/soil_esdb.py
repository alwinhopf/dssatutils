# Soil source: European Soil Database / ESDAC (ESDB), local rasters.
#
# ESDB gives full-profile European soil (topsoil + subsoil), complementing LUCAS
# which is topsoil-only (0-20 cm). Reuses the generic per-property/per-depth raster
# reader. Expected input: GeoTIFF/VRT per property per layer with a property token
# (sand/clay/silt/bulk/organic carbon) AND a layer token (top/sub) in the filename.

from .soil_slga import process_soils_slga

_ESDB_DEPTHS = [("top", 30, 15.0), ("sub", 100, 65.0)]


def process_soils_esdb(
    grid_points, esdb_raster_dir: str, output_csv_path: str, output_sol_dir: str,
    id_col: str = "ID", lat_col: str = "LAT", long_col: str = "LONG",
    depth_specs=None,
) -> None:
    """Build DSSAT .SOL files from local ESDB (Europe) per-property/per-layer rasters."""
    return process_soils_slga(
        grid_points, esdb_raster_dir, output_csv_path, output_sol_dir,
        id_col=id_col, lat_col=lat_col, long_col=long_col,
        depth_specs=depth_specs or _ESDB_DEPTHS,
        source_name="European Soil Database (ESDAC)", source_tag="ESDB local rasters")
