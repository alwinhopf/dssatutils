# Weather source: TerraClimate monthly climate, disaggregated to daily DSSAT input.
#
# TerraClimate is MONTHLY (~4 km); it cannot be written straight to a daily .WTH.
# This backend expands each month into continuous daily records: Tmax/Tmin/SRAD/
# wind are held at the monthly value and the monthly precip TOTAL is spread evenly
# across the days of the month. Runnable but deliberately smooth -- screening /
# climatology only, NOT for day-to-day variability (intensity, dry spells, stress).

# DSSAT var -> (filename tokens, unit-conversion kind, required?)
TERRACLIMATE_WEATHER_VARS <- list(
  TMAX = list(tokens = c("tmmx", "tmax"), kind = "temp", required = TRUE),
  TMIN = list(tokens = c("tmmn", "tmin"), kind = "temp", required = TRUE),
  RAIN = list(tokens = c("ppt", "precip", "pr"), kind = "rain", required = TRUE),
  SRAD = list(tokens = c("srad", "rsds"), kind = "srad", required = FALSE),
  WIND = list(tokens = c("wind", "ws"), kind = "wind", required = FALSE)
)

.terraclimate_expand_daily <- function(monthly_by_var, pid) {
  tmax_m <- monthly_by_var[["TMAX"]][[pid]]
  tmin_m <- monthly_by_var[["TMIN"]][[pid]]
  if (is.null(tmax_m) || is.null(tmin_m)) return(NULL)
  codes <- intersect(names(tmax_m), names(tmin_m))
  if (!length(codes)) return(NULL)
  codes <- sort(codes)
  get1 <- function(v, code) {
    s <- monthly_by_var[[v]][[pid]]
    if (is.null(s) || is.na(s[code])) NA_real_ else as.numeric(s[code])
  }
  frames <- list()
  for (code in codes) {
    year <- as.integer(substr(code, 1, 4))
    doy <- as.integer(substr(code, 5, 7))
    month <- as.integer(format(as.Date(sprintf("%d-01-01", year)) + (doy - 1), "%m"))
    ndays <- as.integer(format(
      as.Date(sprintf("%d-%02d-01", ifelse(month == 12, year + 1, year),
                      ifelse(month == 12, 1, month + 1))) - 1, "%d"))
    srad <- get1("SRAD", code); wind <- get1("WIND", code)
    ppt <- get1("RAIN", code)
    rain_daily <- if (is.finite(ppt)) ppt / ndays else -99
    days <- as.Date(sprintf("%d-%02d-01", year, month)) + (seq_len(ndays) - 1)
    frames[[code]] <- data.frame(
      DATE = sprintf("%d%03d", year, as.integer(format(days, "%j"))),
      YEAR = year, MM = month,
      SRAD = ifelse(is.finite(srad), srad, -99),
      TMAX = as.numeric(tmax_m[code]), TMIN = as.numeric(tmin_m[code]),
      RAIN = rain_daily, TDEW = -99, RH2M = -99,
      WIND = ifelse(is.finite(wind), wind, -99))
  }
  do.call(rbind, frames)
}

process_weather_terraclimate <- function(shapefile, start_year, end_year, output_dir,
                                         id_col, lat_col, lon_col, n_cores, log_file,
                                         terraclimate_nc_dir) {
  if (!nzchar(terraclimate_nc_dir) || !dir.exists(terraclimate_nc_dir))
    stop(sprintf("TerraClimate needs a local NetCDF directory: %s", terraclimate_nc_dir))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  warning("TerraClimate is monthly; disaggregated to daily for screening/climatology only.")
  message(sprintf("--- Starting TerraClimate Processing (Years: %d-%d) ---", start_year, end_year))

  pts <- sf::st_transform(shapefile, 4326)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  pts_vect <- terra::vect(pts)
  xy <- sf::st_coordinates(pts); lats <- xy[, 2]; lons <- xy[, 1]

  per_var <- list()
  for (v in names(TERRACLIMATE_WEATHER_VARS)) {
    spec <- TERRACLIMATE_WEATHER_VARS[[v]]
    path <- weather_find_nc_file(terraclimate_nc_dir, spec$tokens)
    if (is.na(path)) {
      if (isTRUE(spec$required))
        stop(sprintf("TerraClimate required variable %s not found in %s", v, terraclimate_nc_dir))
      next
    }
    per_var[[v]] <- weather_extract_netcdf_series(path, ids, pts_vect,
                                                  start_year, end_year, spec$kind)
  }

  written <- 0
  for (k in seq_along(ids)) {
    pid <- ids[k]
    tryCatch({
      df <- .terraclimate_expand_daily(per_var, pid)
      if (is.null(df) || !nrow(df)) stop("No overlapping monthly TMAX/TMIN extracted for this point.")
      weather_write_wth(df, pid, lats[k], lons[k], output_dir, "TerraClimate monthly->daily", "TCLM")
      written <- written + 1
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\nTerraClimate point %s: %s\n", pid, conditionMessage(e))
      cat(msg); write(msg, file = log_file, append = TRUE)
    })
  }
  message(sprintf("\nTerraClimate processing complete: %d point(s) written.\n", written))
}
