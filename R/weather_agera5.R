# File: weather_agera5.R
# ---------------------------------------------------------------------------
# Weather source: AgERA5 (ECMWF agrometeorological reanalysis) -> DSSAT .WTH.
#
# AgERA5 is ERA5 reprocessed for agriculture: global, 0.1deg (~10 km), daily,
# 1979-present, with the daily statistics crop models need (24h max/min/mean
# temperature, solar radiation flux, precipitation flux, RH, wind, dewpoint).
# Covers the poles (unlike CHIRPS) and is higher-res than NASA POWER.
#
# ACCESS (requires a free key, NOT keyless):
#   1. Register at the Copernicus CDS: https://cds.climate.copernicus.eu/
#   2. Store your key (ecmwfr::wf_set_key) or via ~/.cdsapirc.
#   3. install.packages(c("ecmwfr","terra"))
#   Dataset: "sis-agrometeorological-indicators"
#
# Requests are queued by the CDS, so first runs can be slow; downloads are
# cached under `agera5_cache_dir`. Mirrors the Python weather_agera5.py.
#
# NOTE: provided for parity; validate once against your CDS account. Requires
# ecmwfr + terra. Same .WTH format as the NASA POWER / Open-Meteo writers.
# ---------------------------------------------------------------------------

library(lubridate)
library(dplyr)
library(terra)

# AgERA5 CDS variable -> (variable, selector kind+value) and DSSAT unit handling.
# 2m_relative_humidity uses a fixed-hour `time` selector (NOT a 24-hour
# statistic); fluxes take no selector. sel_kind is "statistic", "time", or NA.
.agera5_vars <- list(
  TMAX = list(var = "2m_temperature",          sel_kind = "statistic", sel = "24_hour_maximum"),  # K->C
  TMIN = list(var = "2m_temperature",          sel_kind = "statistic", sel = "24_hour_minimum"),  # K->C
  SRAD = list(var = "solar_radiation_flux",    sel_kind = NA,          sel = NA),                 # J/m2->MJ
  RAIN = list(var = "precipitation_flux",      sel_kind = NA,          sel = NA),                 # mm/day
  TDEW = list(var = "2m_dewpoint_temperature", sel_kind = "statistic", sel = "24_hour_mean"),     # K->C
  RH2M = list(var = "2m_relative_humidity",    sel_kind = "time",      sel = "15_00"),            # %  mid-afternoon
  WIND = list(var = "10m_wind_speed",          sel_kind = "statistic", sel = "24_hour_mean")      # m/s
)

process_weather_agera5 <- function(shapefile, start_year, end_year, output_dir,
                                   id_col, lat_col, lon_col, n_cores, log_file,
                                   agera5_cache_dir) {
  if (!requireNamespace("ecmwfr", quietly = TRUE))
    stop("AgERA5 needs the 'ecmwfr' package + a Copernicus CDS key. install.packages('ecmwfr')")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(agera5_cache_dir)) dir.create(agera5_cache_dir, recursive = TRUE)

  ids  <- as.character(shapefile[[id_col]])
  lats <- as.numeric(shapefile[[lat_col]]); lons <- as.numeric(shapefile[[lon_col]])
  end_year <- min(end_year, lubridate::year(Sys.Date()))

  message(sprintf("--- Starting AgERA5 Download (Years: %d-%d) ---", start_year, end_year))
  message("  NOTE: AgERA5 requires a Copernicus CDS API key and queues requests; first run can be slow.")

  pad <- 0.2
  area <- c(max(lats) + pad, min(lons) - pad, min(lats) - pad, max(lons) + pad)  # N,W,S,E
  pts <- terra::vect(data.frame(lon = lons, lat = lats), geom = c("lon", "lat"),
                     crs = "EPSG:4326")

  # point_series[[pid]][[VAR]] is a named numeric vector keyed by "YYYYDOY".
  point_series <- setNames(lapply(ids, function(x)
    setNames(vector("list", length(.agera5_vars)), names(.agera5_vars))), ids)

  for (yr in start_year:end_year) {
    for (vname in names(.agera5_vars)) {
      spec <- .agera5_vars[[vname]]
      tag <- sprintf("%s_%s_%d", spec$var, ifelse(is.na(spec$sel), "na", spec$sel), yr)
      dest <- file.path(agera5_cache_dir, sprintf("agera5_%s.nc", tag))
      if (!file.exists(dest)) {
        req <- list(dataset_short_name = "sis-agrometeorological-indicators",
                    variable = spec$var, year = as.character(yr),
                    month = sprintf("%02d", 1:12), day = sprintf("%02d", 1:31),
                    area = area, version = "2_0",  # AgERA5 v2 (v1.1 deprecated 2026-06)
                    target = basename(dest))
        if (!is.na(spec$sel_kind)) req[[spec$sel_kind]] <- spec$sel  # statistic OR time
        tryCatch(
          ecmwfr::wf_request(request = req, path = agera5_cache_dir),
          error = function(e) message(sprintf("  AgERA5 download failed (%s): %s",
                                              tag, conditionMessage(e))))
      }
      if (!file.exists(dest)) next
      tryCatch({
        r <- terra::rast(dest)
        tvals <- terra::time(r)
        date_codes <- sprintf("%d%03d", lubridate::year(tvals), lubridate::yday(tvals))
        ex <- terra::extract(r, pts, ID = FALSE)
        for (j in seq_along(ids)) {
          v <- as.numeric(ex[j, ])
          if (vname %in% c("TMAX", "TMIN", "TDEW")) v <- v - 273.15
          if (vname == "SRAD") v <- v * 1e-6
          names(v) <- date_codes
          point_series[[ids[j]]][[vname]] <- c(point_series[[ids[j]]][[vname]], v)
        }
      }, error = function(e)
        message(sprintf("  AgERA5 extract failed (%s): %s", tag, conditionMessage(e))))
    }
  }

  written <- 0
  for (i in seq_along(ids)) {
    pid <- ids[i]
    tryCatch({
      ps <- point_series[[pid]]
      if (is.null(ps$TMAX) || length(ps$TMAX) == 0) stop("No AgERA5 data for point.")
      dates <- names(ps$TMAX)
      get <- function(var) as.numeric(ps[[var]][dates])
      wd <- data.frame(DATE = dates, SRAD = get("SRAD"), TMAX = get("TMAX"),
                       TMIN = get("TMIN"), RAIN = get("RAIN"), TDEW = get("TDEW"),
                       RH2M = get("RH2M"), WIND = get("WIND"))
      wd[is.na(wd)] <- -99
      wd$YEAR <- as.integer(substr(wd$DATE, 1, 4))
      wd$MM <- lubridate::month(as.Date(wd$DATE, format = "%Y%j"))
      wd$TAVG <- (wd$TMAX + wd$TMIN) / 2
      tav <- mean(wd$TAVG, na.rm = TRUE)
      monthly <- wd %>% group_by(YEAR, MM) %>% summarise(m = mean(TAVG), .groups = "drop")
      amp <- mean((monthly %>% group_by(YEAR) %>% summarise(a = max(m) - min(m)))$a, na.rm = TRUE)

      hdr <- sprintf(
        "$WEATHER DATA: AgERA5 (Point ID: %s)\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  AGE5 %8.4f %8.4f   -99 %5.1f %5.1f   2.0  10.0\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
        pid, lats[i], lons[i], tav, amp)
      lines <- with(wd, sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                                DATE, SRAD, TMAX, TMIN, RAIN, TDEW, RH2M, WIND))
      lines <- gsub("-99.0", "  -99", lines, fixed = TRUE)
      writeLines(c(hdr, lines), con = file.path(output_dir, sprintf("%s.WTH", pid)))
      written <- written + 1
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\nAgERA5 point %s: %s\n", pid, conditionMessage(e))
      cat(msg); write(msg, file = log_file, append = TRUE)
    })
  }
  message(sprintf("\nAgERA5 processing complete: %d/%d points written to '%s'.\n",
                  written, length(ids), output_dir))
}
