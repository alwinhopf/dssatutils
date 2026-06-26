SILO_WEATHER_VARS <- list(
  TMAX = list(tokens = c("max_temp", "tmax", "maxt"), kind = "temp", required = TRUE),
  TMIN = list(tokens = c("min_temp", "tmin", "mint"), kind = "temp", required = TRUE),
  RAIN = list(tokens = c("rain", "precip", "ppt"), kind = "rain", required = TRUE),
  SRAD = list(tokens = c("radiation", "srad", "solar"), kind = "srad", required = TRUE),
  # SILO distributes vapour pressure ("vp", hPa), not dewpoint: convert via kind="vp".
  TDEW = list(tokens = c("vp"), kind = "vp", required = FALSE),
  WIND = list(tokens = c("wind"), kind = "wind", required = FALSE)
)

process_weather_silo <- function(shapefile, start_year, end_year, output_dir,
                                 id_col, lat_col, lon_col, n_cores, log_file,
                                 silo_nc_dir) {
  message(sprintf("--- Starting SILO Processing (Years: %d-%d) ---", start_year, end_year))
  written <- process_local_netcdf_weather(shapefile, start_year, end_year, output_dir,
                                          id_col, lat_col, lon_col, log_file,
                                          silo_nc_dir, SILO_WEATHER_VARS,
                                          "SILO Australia", "SILO")
  message(sprintf("\nSILO processing complete: %d point(s) written.\n", written))
}
