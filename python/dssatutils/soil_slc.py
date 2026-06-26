# Soil source: Soil Landscapes of Canada (SLC), local rasters.
#
# SLC (Agriculture & Agri-Food Canada) pairs with the ANUSPLIN weather backend to
# give a Canada stack. SLC is natively polygon map units + component attribute
# tables; rasterize that join (per property, top/sub layers) before use. Expected
# input: GeoTIFF/VRT per property per layer with a property token (sand/clay/silt/
# bulk/organic carbon) AND a layer token (top/sub) in the filename.

from .soil_slga import process_soils_slga

# SLC components are typically summarised as topsoil/subsoil horizons.
_SLC_DEPTHS = [("top", 30, 15.0), ("sub", 100, 65.0)]


def process_soils_slc(
    grid_points, slc_raster_dir: str, output_csv_path: str, output_sol_dir: str,
    id_col: str = "ID", lat_col: str = "LAT", long_col: str = "LONG",
    depth_specs=None,
) -> None:
    """Build DSSAT .SOL files from local SLC (Canada) per-property/per-layer rasters."""
    return process_soils_slga(
        grid_points, slc_raster_dir, output_csv_path, output_sol_dir,
        id_col=id_col, lat_col=lat_col, long_col=long_col,
        depth_specs=depth_specs or _SLC_DEPTHS,
        source_name="Soil Landscapes of Canada", source_tag="SLC local rasters")
