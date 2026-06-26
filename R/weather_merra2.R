# Weather source: NASA MERRA-2 reanalysis daily aggregates (local NetCDF).
# Global ~0.5 x 0.625 degree, 1980-present. Temperatures in Kelvin auto-convert.

MERRA2_WEATHER_VARS <- list(
  TMAX = list(tokens = c("t2mmax", "tmax", "tasmax"), kind = "temp", required = TRUE),
  TMIN = list(tokens = c("t2mmin", "tmin", "tasmin"), kind = "temp", required = TRUE),
  TMEAN = list(tokens = c("t2mmean", "t2m", "tas"), kind = "temp", required = FALSE),
  RAIN = list(tokens = c("prectot", "precip", "pr"), kind = "rain", required = TRUE),
  SRAD = list(tokens = c("swgdn", "swgnt", "rsds", "srad"), kind = "srad", required = TRUE),
  WIND = list(tokens = c("speed", "wind", "sfcwind"), kind = "wind", required = FALSE),
  RH2M = list(tokens = c("rh2m", "rh", "hurs"), kind = "rh", required = FALSE)
)

process_weather_merra2 <- function(shapefile, start_year, end_year, output_dir,
                                   id_col, lat_col, lon_col, n_cores, log_file,
                                   merra2_nc_dir) {
  start_year <- max(as.integer(start_year), 1980)
  message(sprintf("--- Starting MERRA-2 Processing (Years: %d-%d) ---", start_year, end_year))
  written <- process_local_netcdf_weather(shapefile, start_year, end_year, output_dir,
                                          id_col, lat_col, lon_col, log_file,
                                          merra2_nc_dir, MERRA2_WEATHER_VARS,
                                          "NASA MERRA-2", "MER2")
  message(sprintf("\nMERRA-2 processing complete: %d point(s) written.\n", written))
}
