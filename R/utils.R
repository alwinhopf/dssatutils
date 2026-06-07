# Shared package-private helper functions for dssatutils

# Extract IDs, latitudes, and longitudes from sf object or data.frame
.extract_coords <- function(shapefile, id_col, lat_col, lon_col) {
  ids <- as.character(shapefile[[id_col]])
  if (inherits(shapefile, "sf")) {
    shapefile_wgs84 <- sf::st_transform(shapefile, 4326)
    coords <- sf::st_coordinates(shapefile_wgs84)
    lons <- as.numeric(coords[, 1])
    lats <- as.numeric(coords[, 2])
  } else {
    lats <- as.numeric(shapefile[[lat_col]])
    lons <- as.numeric(shapefile[[lon_col]])
  }
  list(ids = ids, lats = lats, lons = lons)
}
