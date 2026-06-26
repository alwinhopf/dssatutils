process_soils_wosis <- function(grid_points, wosis_profile_csv, output_csv_path,
                                output_sol_dir, id_col = "ID",
                                lat_col = "LAT", long_col = "LONG") {
  if (!file.exists(wosis_profile_csv)) stop(wosis_profile_csv)
  src <- utils::read.csv(wosis_profile_csv, stringsAsFactors = FALSE)
  req <- c("profile_id", "latitude", "longitude", "depth_bottom", "sand", "clay", "silt", "bdod", "soc_pct")
  miss <- setdiff(req, names(src))
  if (length(miss)) stop(sprintf("WoSIS processed CSV missing columns: %s", paste(miss, collapse = ", ")))
  pts <- sf::st_transform(grid_points, 4326)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  heads <- src[!duplicated(src$profile_id), ]
  src_sf <- sf::st_as_sf(heads, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  nearest <- sf::st_nearest_feature(pts, src_sf)
  rows <- list()
  for (i in seq_along(ids)) {
    sid <- heads$profile_id[nearest[i]]
    sub <- src[src$profile_id == sid, , drop = FALSE]
    sub$ID <- ids[i]
    if (!"depth_center" %in% names(sub)) sub$depth_center <- sub$depth_bottom / 2
    if (!"cfvo" %in% names(sub)) sub$cfvo <- 0
    rows[[i]] <- sub
  }
  df <- do.call(rbind, rows)
  df <- soil_add_physics(df)
  soil_write_mapping(ids, output_csv_path)
  soil_write_profiles(df, output_sol_dir, "WoSIS processed point profiles", "nearest processed profile")
}
