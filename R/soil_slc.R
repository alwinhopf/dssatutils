# Soil source: Soil Landscapes of Canada (SLC), local rasters.
# Pairs with ANUSPLIN. SLC is natively polygon map units + component tables;
# rasterize that join (top/sub layers) before use. Rasters need a property token
# (sand/clay/silt/bulk/organic carbon) AND a layer token (top/sub) in the filename.

SLC_DEPTHS <- data.frame(
  token = c("top", "sub"), bottom = c(30, 100), center = c(15, 65)
)

process_soils_slc <- function(grid_points, slc_raster_dir, output_csv_path,
                              output_sol_dir, id_col = "ID", lat_col = "LAT",
                              long_col = "LONG", depth_specs = NULL) {
  process_soils_slga(grid_points, slc_raster_dir, output_csv_path, output_sol_dir,
                     id_col = id_col, lat_col = lat_col, long_col = long_col,
                     depth_specs = if (is.null(depth_specs)) SLC_DEPTHS else depth_specs,
                     source_name = "Soil Landscapes of Canada",
                     source_tag = "SLC local rasters")
}
