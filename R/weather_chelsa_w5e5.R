CHELSA_W5E5_VARS <- list(
  TMAX = list(tokens = c("tasmax"), kind = "temp", required = TRUE),
  TMIN = list(tokens = c("tasmin"), kind = "temp", required = TRUE),
  TMEAN = list(tokens = c("tas_", "tas."), kind = "temp", required = FALSE),
  RAIN = list(tokens = c("pr"), kind = "rain", required = TRUE),
  SRAD = list(tokens = c("rsds"), kind = "srad", required = TRUE)
)

process_weather_chelsa_w5e5 <- function(shapefile, start_year, end_year, output_dir,
                                         id_col, lat_col, lon_col, n_cores, log_file,
                                         chelsa_nc_dir) {
  end_year <- min(as.integer(end_year), 2016)
  message(sprintf("--- Starting CHELSA-W5E5 Processing (Years: %d-%d) ---", start_year, end_year))
  written <- process_local_netcdf_weather(shapefile, start_year, end_year, output_dir,
                                          id_col, lat_col, lon_col, log_file,
                                          chelsa_nc_dir, CHELSA_W5E5_VARS,
                                          "CHELSA-W5E5", "CHW5")
  message(sprintf("\nCHELSA-W5E5 processing complete: %d point(s) written.\n", written))
}
