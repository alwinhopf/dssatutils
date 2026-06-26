CRUJRA_WEATHER_VARS <- list(
  TMAX = list(tokens = c("tmax", "tasmax"), kind = "temp", required = TRUE),
  TMIN = list(tokens = c("tmin", "tasmin"), kind = "temp", required = TRUE),
  TMEAN = list(tokens = c("tas", "tmean"), kind = "temp", required = FALSE),
  RAIN = list(tokens = c("pr", "precip"), kind = "rain", required = TRUE),
  SRAD = list(tokens = c("rsds", "srad"), kind = "srad", required = TRUE),
  WIND = list(tokens = c("wind", "sfcwind"), kind = "wind", required = FALSE),
  RH2M = list(tokens = c("rh", "hurs"), kind = "rh", required = FALSE)
)

process_weather_crujra <- function(shapefile, start_year, end_year, output_dir,
                                   id_col, lat_col, lon_col, n_cores, log_file,
                                   crujra_nc_dir) {
  message(sprintf("--- Starting CRU-JRA Processing (Years: %d-%d) ---", start_year, end_year))
  written <- process_local_netcdf_weather(shapefile, start_year, end_year, output_dir,
                                          id_col, lat_col, lon_col, log_file,
                                          crujra_nc_dir, CRUJRA_WEATHER_VARS,
                                          "CRU-JRA", "CRUJ")
  message(sprintf("\nCRU-JRA processing complete: %d point(s) written.\n", written))
}
