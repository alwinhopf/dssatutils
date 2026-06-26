process_soils_slga <- function(grid_points, slga_raster_dir, output_csv_path,
                               output_sol_dir, id_col = "ID",
                               lat_col = "LAT", long_col = "LONG",
                               depth_specs = NULL,
                               source_name = "Soil and Landscape Grid of Australia",
                               source_tag = "SLGA local rasters") {
  if (!dir.exists(slga_raster_dir)) stop(sprintf("SLGA raster directory not found: %s", slga_raster_dir))
  if (is.null(depth_specs)) depth_specs <- SOIL_RASTER_DEPTHS
  pts <- sf::st_transform(grid_points, 4326)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  xy <- sf::st_coordinates(pts); lats <- xy[, 2]; lons <- xy[, 1]
  pts_vect <- terra::vect(pts)
  rows <- list()
  for (j in seq_len(nrow(depth_specs))) {
    d <- depth_specs[j, ]
    sand <- soil_sample_raster(soil_find_raster(slga_raster_dir, "sand", d$token), pts_vect)
    clay <- soil_sample_raster(soil_find_raster(slga_raster_dir, "clay", d$token), pts_vect)
    silt_p <- soil_find_raster(slga_raster_dir, "silt", d$token)
    silt <- if (!is.na(silt_p)) soil_sample_raster(silt_p, pts_vect) else 100 - sand - clay
    bd <- soil_sample_raster(soil_find_raster(slga_raster_dir, c("bulk"), d$token), pts_vect)
    bd <- ifelse(bd > 10, bd / 100, bd)
    soc <- soil_sample_raster(soil_find_raster(slga_raster_dir, c("organic", "carbon"), d$token), pts_vect)
    soc_pct <- ifelse(soc > 20, soc / 10, soc)
    rows[[j]] <- data.frame(ID = ids, latitude = lats, longitude = lons,
                            depth_bottom = d$bottom, depth_center = d$center,
                            sand = sand, clay = clay, silt = silt,
                            bdod = bd, soc_pct = soc_pct, cfvo = 0)
  }
  df <- do.call(rbind, rows)
  df <- df[complete.cases(df[, c("sand", "clay", "bdod", "soc_pct")]), ]
  if (!nrow(df)) stop("No usable SLGA soil data extracted. Check raster names and coordinates.")
  df <- soil_add_physics(df)
  soil_write_mapping(ids, output_csv_path)
  soil_write_profiles(df, output_sol_dir, source_name, source_tag)
}
