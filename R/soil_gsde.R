# Soil source: GSDE (Global Soil Dataset for Earth System Modeling), local rasters.
# Global 1 km, 8 layers to ~2.3 m. Rasters need a property token (sand/clay/silt/
# bulk/organic carbon) AND a layer token (l1..l8) in the filename.

GSDE_DEPTHS <- data.frame(
  token = c("l1", "l2", "l3", "l4", "l5", "l6", "l7", "l8"),
  bottom = c(5, 9, 17, 29, 49, 83, 138, 230),
  center = c(2, 7, 13, 23, 39, 66, 110, 184)
)

process_soils_gsde <- function(grid_points, gsde_raster_dir, output_csv_path,
                               output_sol_dir, id_col = "ID", lat_col = "LAT",
                               long_col = "LONG", depth_specs = NULL) {
  process_soils_slga(grid_points, gsde_raster_dir, output_csv_path, output_sol_dir,
                     id_col = id_col, lat_col = lat_col, long_col = long_col,
                     depth_specs = if (is.null(depth_specs)) GSDE_DEPTHS else depth_specs,
                     source_name = "GSDE (Global Soil Dataset for ESM)",
                     source_tag = "GSDE local rasters")
}
