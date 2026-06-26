# Soil source: European Soil Database / ESDAC (ESDB), local rasters.
# Full-profile Europe (topsoil + subsoil), complementing topsoil-only LUCAS.
# Rasters need a property token (sand/clay/silt/bulk/organic carbon) AND a layer
# token (top/sub) in the filename.

ESDB_DEPTHS <- data.frame(
  token = c("top", "sub"), bottom = c(30, 100), center = c(15, 65)
)

process_soils_esdb <- function(grid_points, esdb_raster_dir, output_csv_path,
                               output_sol_dir, id_col = "ID", lat_col = "LAT",
                               long_col = "LONG", depth_specs = NULL) {
  process_soils_slga(grid_points, esdb_raster_dir, output_csv_path, output_sol_dir,
                     id_col = id_col, lat_col = lat_col, long_col = long_col,
                     depth_specs = if (is.null(depth_specs)) ESDB_DEPTHS else depth_specs,
                     source_name = "European Soil Database (ESDAC)",
                     source_tag = "ESDB local rasters")
}
