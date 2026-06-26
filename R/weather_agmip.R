AGMIP_WEATHER_VARS <- list(
  TMAX = list(tokens = c("tmax", "tasmax"), kind = "temp", required = TRUE),
  TMIN = list(tokens = c("tmin", "tasmin"), kind = "temp", required = TRUE),
  TMEAN = list(tokens = c("tmean", "tas"), kind = "temp", required = FALSE),
  RAIN = list(tokens = c("prate", "precip", "pr"), kind = "rain", required = TRUE),
  SRAD = list(tokens = c("srad", "rsds", "swdown"), kind = "srad", required = TRUE),
  WIND = list(tokens = c("wind", "sfcwind", "wnd"), kind = "wind", required = FALSE),
  TDEW = list(tokens = c("tdew", "dew"), kind = "temp", required = FALSE),
  RH2M = list(tokens = c("rh", "rhum", "hurs"), kind = "rh", required = FALSE)
)

.process_weather_agmip <- function(shapefile, start_year, end_year, output_dir,
                                   id_col, lat_col, lon_col, log_file,
                                   nc_dir, product, insi) {
  start_year <- max(as.integer(start_year), 1980)
  end_year <- min(as.integer(end_year), 2010)
  message(sprintf("--- Starting %s Processing (Years: %d-%d) ---", product, start_year, end_year))
  written <- process_local_netcdf_weather(shapefile, start_year, end_year, output_dir,
                                          id_col, lat_col, lon_col, log_file,
                                          nc_dir, AGMIP_WEATHER_VARS, product, insi)
  message(sprintf("\n%s processing complete: %d point(s) written.\n", product, written))
}

process_weather_agmerra <- function(shapefile, start_year, end_year, output_dir,
                                    id_col, lat_col, lon_col, n_cores, log_file,
                                    agmerra_nc_dir) {
  .process_weather_agmip(shapefile, start_year, end_year, output_dir,
                         id_col, lat_col, lon_col, log_file,
                         agmerra_nc_dir, "AgMERRA", "AGMR")
}

process_weather_agcfsr <- function(shapefile, start_year, end_year, output_dir,
                                   id_col, lat_col, lon_col, n_cores, log_file,
                                   agcfsr_nc_dir) {
  .process_weather_agmip(shapefile, start_year, end_year, output_dir,
                         id_col, lat_col, lon_col, log_file,
                         agcfsr_nc_dir, "AgCFSR", "AGCF")
}
