# Weather source: Princeton Global Forcing (PGF) daily reanalysis (local NetCDF).
# Global 0.25 degree; full-variable alternative to AgMERRA/CRU-JRA.

PGF_WEATHER_VARS <- list(
  TMAX = list(tokens = c("tmax", "tasmax"), kind = "temp", required = TRUE),
  TMIN = list(tokens = c("tmin", "tasmin"), kind = "temp", required = TRUE),
  TMEAN = list(tokens = c("tas", "tmean"), kind = "temp", required = FALSE),
  RAIN = list(tokens = c("prcp", "precip", "pr"), kind = "rain", required = TRUE),
  SRAD = list(tokens = c("dswrf", "rsds", "srad", "swdown"), kind = "srad", required = TRUE),
  WIND = list(tokens = c("wind", "wnd", "sfcwind"), kind = "wind", required = FALSE),
  RH2M = list(tokens = c("rh", "hurs"), kind = "rh", required = FALSE)
)

process_weather_pgf <- function(shapefile, start_year, end_year, output_dir,
                                id_col, lat_col, lon_col, n_cores, log_file,
                                pgf_nc_dir) {
  message(sprintf("--- Starting PGF Processing (Years: %d-%d) ---", start_year, end_year))
  written <- process_local_netcdf_weather(shapefile, start_year, end_year, output_dir,
                                          id_col, lat_col, lon_col, log_file,
                                          pgf_nc_dir, PGF_WEATHER_VARS,
                                          "Princeton Global Forcing", "PGF")
  message(sprintf("\nPGF processing complete: %d point(s) written.\n", written))
}
