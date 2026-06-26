# Soil source: WISE30sec derived soil properties (ISRIC), local rasters.
#
# WISE30sec (Batjes 2016) reports properties on 7 fixed depth layers, NOT the
# SoilGrids 0-5/5-15/... intervals. Reuse the generic per-property/per-depth
# raster reader but with WISE30sec's own depth table. Expected input: one
# GeoTIFF/VRT per property per depth layer, with a property token (sand/clay/
# silt/bulk/organic carbon) AND a depth token (d1..d7) in the filename.

# WISE30sec layers D1..D7: 0-20, 20-40, 40-60, 60-80, 80-100, 100-150, 150-200 cm.
WISE30SEC_DEPTHS <- data.frame(
  token = c("d1", "d2", "d3", "d4", "d5", "d6", "d7"),
  bottom = c(20, 40, 60, 80, 100, 150, 200),
  center = c(10, 30, 50, 70, 90, 125, 175)
)

process_soils_wise30sec <- function(grid_points, wise30sec_raster_dir,
                                    output_csv_path, output_sol_dir,
                                    id_col = "ID", lat_col = "LAT",
                                    long_col = "LONG", depth_specs = NULL) {
  process_soils_slga(grid_points, wise30sec_raster_dir, output_csv_path,
                     output_sol_dir, id_col = id_col,
                     lat_col = lat_col, long_col = long_col,
                     depth_specs = if (is.null(depth_specs)) WISE30SEC_DEPTHS else depth_specs,
                     source_name = "ISRIC WISE30sec",
                     source_tag = "WISE30sec local rasters")
}
