# Weather source: ANUSPLIN (Agriculture & Agri-Food Canada) daily 10 km grids.
# All of Canada, 1950-2015 (maxt/mint/pcp); SRAD/RH/wind written -99.

ANUSPLIN_WEATHER_VARS <- list(
  TMAX = list(tokens = c("maxt", "tmax", "tasmax"), kind = "temp", required = TRUE),
  TMIN = list(tokens = c("mint", "tmin", "tasmin"), kind = "temp", required = TRUE),
  RAIN = list(tokens = c("pcp", "precip", "pr", "rain"), kind = "rain", required = TRUE)
)

process_weather_anusplin <- function(shapefile, start_year, end_year, output_dir,
                                     id_col, lat_col, lon_col, n_cores, log_file,
                                     anusplin_nc_dir) {
  end_year <- min(as.integer(end_year), 2015)
  message(sprintf("--- Starting ANUSPLIN Processing (Years: %d-%d) ---", start_year, end_year))
  written <- process_local_netcdf_weather(shapefile, start_year, end_year, output_dir,
                                          id_col, lat_col, lon_col, log_file,
                                          anusplin_nc_dir, ANUSPLIN_WEATHER_VARS,
                                          "ANUSPLIN Canada", "ANUS")
  message(sprintf("\nANUSPLIN processing complete: %d point(s) written.\n", written))
}
