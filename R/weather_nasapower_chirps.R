# File: weather_nasapower_chirps.R
# ---------------------------------------------------------------------------
# HYBRID weather source: NASA POWER (all variables) + CHIRPS (rainfall).
#
# NASA POWER is global and provides every DSSAT variable, but its rainfall is
# coarse (~0.5deg). CHIRPS is daily, ~0.05deg, 1981-present, station-blended,
# and far better for precipitation (especially tropics / semi-arid: Africa,
# India). This fetches NASA POWER per point and REPLACES the RAIN column with
# CHIRPS within its 50S-50N coverage, falling back to NASA POWER rain outside
# that band or over CHIRPS no-data cells, so output stays global.
#
# CHIRPS netCDF: https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/
# Cite: Funk et al. 2015, Sci. Data 2:150066.
#
# Mirrors process_weather_nasapower() but takes an extra `chirps_cache_dir`
# argument (like GridMET's cache dir). Resolution via global CHIRPS_RESOLUTION
# ("p05" default or "p25").
# ---------------------------------------------------------------------------

library(nasapower)
library(lubridate)
library(dplyr)
library(terra)

process_weather_nasapower_chirps <- function(shapefile, start_year, end_year,
                                             output_dir, id_col, lat_col, lon_col,
                                             n_cores, log_file, chirps_cache_dir) {

  res <- if (exists("CHIRPS_RESOLUTION")) CHIRPS_RESOLUTION else "p05"
  chirps_lat_limit <- 50.0
  chirps_nodata <- -9999

  message(sprintf("--- Starting NASA-POWER + CHIRPS(%s) Hybrid (Years: %d-%d) ---",
                  res, start_year, end_year))

  start_date_str <- paste0(start_year, "-01-01")
  current_year <- lubridate::year(Sys.Date())
  if (end_year == current_year) {
    end_date_str <- as.character(Sys.Date() - 2)
    message(sprintf("End year is current year. Fetching up to: %s", end_date_str))
  } else {
    end_date_str <- paste0(end_year, "-12-31")
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(chirps_cache_dir, showWarnings = FALSE, recursive = TRUE)

  ids  <- as.character(shapefile[[id_col]])
  lats <- as.numeric(shapefile[[lat_col]])
  lons <- as.numeric(shapefile[[lon_col]])

  # --- 1. CHIRPS: download yearly netCDFs and extract per-point daily rain ---
  # chirps_rain[[point_id]] is a named numeric vector keyed by "YYYYDOY".
  chirps_rain <- setNames(vector("list", length(ids)), ids)
  if (any(abs(lats) <= chirps_lat_limit)) {
    base_url <- "https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf"
    pts <- terra::vect(data.frame(lon = lons, lat = lats),
                       geom = c("lon", "lat"), crs = "EPSG:4326")
    for (yr in start_year:end_year) {
      fname <- sprintf("chirps-v2.0.%d.days_%s.nc", yr, res)
      dest  <- file.path(chirps_cache_dir, fname)
      if (!file.exists(dest) || file.info(dest)$size == 0) {
        url <- sprintf("%s/%s/%s", base_url, res, fname)
        message(sprintf("  Downloading CHIRPS %d (%s)... (large; cached after first run)", yr, res))
        tryCatch(
          utils::download.file(url, dest, mode = "wb", quiet = TRUE),
          error = function(e) message(sprintf("  CHIRPS %d download failed: %s",
                                              yr, conditionMessage(e))))
      }
      if (!file.exists(dest) || file.info(dest)$size == 0) next
      tryCatch({
        r <- terra::rast(dest)                      # one layer per day
        tvals <- terra::time(r)
        if (all(is.na(tvals))) {
          # Fall back to layer names if time() is unset.
          tvals <- as.Date(gsub(".*?(\\d{4}\\.\\d{2}\\.\\d{2}).*", "\\1",
                                names(r)), format = "%Y.%m.%d")
        }
        date_codes <- sprintf("%d%03d", lubridate::year(tvals),
                              lubridate::yday(tvals))
        ex <- terra::extract(r, pts, ID = FALSE)     # rows=points, cols=days
        ex <- as.matrix(ex)
        ex[ex <= chirps_nodata] <- NA
        for (j in seq_along(ids)) {
          vals <- ex[j, ]
          keep <- !is.na(vals)
          if (any(keep)) {
            v <- setNames(as.numeric(vals[keep]), date_codes[keep])
            chirps_rain[[ids[j]]] <- c(chirps_rain[[ids[j]]], v)
          }
        }
      }, error = function(e)
        message(sprintf("  CHIRPS %d extraction failed: %s", yr, conditionMessage(e))))
    }
  } else {
    message("  All points outside CHIRPS coverage (|lat| > 50); using NASA-POWER rain.")
  }

  nasa_params <- c("T2M_MAX", "T2M_MIN", "ALLSKY_SFC_SW_DWN", "PRECTOTCORR",
                   "T2MDEW", "RH2M", "WS2M")

  message(sprintf("Processing %d point(s) (NASA-POWER + CHIRPS merge)...", length(ids)))

  # --- 2. NASA POWER per point + CHIRPS rain merge (serial; get_power is the
  #        network bottleneck and nasapower self-throttles) ---
  for (i in seq_len(nrow(shapefile))) {
    latitude  <- lats[i]; longitude <- lons[i]; point_id <- ids[i]
    output_file <- file.path(output_dir, sprintf("%s.WTH", point_id))
    if (file.exists(output_file)) next

    tryCatch({
      power_data <- nasapower::get_power(
        community = "AG", lonlat = c(longitude, latitude),
        pars = nasa_params, dates = c(start_date_str, end_date_str),
        temporal_api = "DAILY")
      if (nrow(power_data) == 0) stop("No data returned from NASA-POWER.")

      weather_data <- power_data %>%
        dplyr::rename(SRAD = ALLSKY_SFC_SW_DWN, TMAX = T2M_MAX, TMIN = T2M_MIN,
                      RAIN = PRECTOTCORR, TDEW = T2MDEW, RH2M = RH2M, WIND = WS2M) %>%
        dplyr::mutate(DATE = sprintf("%d%03d", YEAR, DOY)) %>%
        dplyr::mutate(across(where(is.numeric), ~ ifelse(. == -999, -99, .)))

      # Merge CHIRPS rainfall over NASA-POWER rain (within coverage band).
      n_chirps <- 0
      cr <- chirps_rain[[point_id]]
      if (abs(latitude) <= chirps_lat_limit && !is.null(cr) && length(cr) > 0) {
        idx <- match(weather_data$DATE, names(cr))
        hit <- !is.na(idx)
        weather_data$RAIN[hit] <- as.numeric(cr[idx[hit]])
        n_chirps <- sum(hit)
      }
      rain_src <- if (n_chirps > 0)
        sprintf("CHIRPS(%s) where available, %d days; NASA-POWER otherwise", res, n_chirps)
      else "NASA-POWER (CHIRPS unavailable here)"

      weather_data$TAVG <- (weather_data$TMAX + weather_data$TMIN) / 2
      tav <- mean(weather_data$TAVG, na.rm = TRUE)
      monthly_temps <- weather_data %>% dplyr::group_by(YEAR, MM) %>%
        dplyr::summarise(TAVG_MON = mean(TAVG, na.rm = TRUE), .groups = "drop")
      annual_amps <- monthly_temps %>% dplyr::group_by(YEAR) %>%
        dplyr::summarise(AMP_YR = max(TAVG_MON) - min(TAVG_MON), .groups = "drop")
      amp <- mean(annual_amps$AMP_YR, na.rm = TRUE)

      wth_header <- sprintf(
        "$WEATHER DATA: NASA-POWER + CHIRPS rain (Point ID: %s) [%s]\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  NAPC %8.4f %8.4f   -99 %5.1f %5.1f   2.0   2.0\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
        point_id, rain_src, latitude, longitude, tav, amp)

      weather_lines <- with(weather_data, sprintf(
        "%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
        DATE, SRAD, TMAX, TMIN, RAIN, TDEW, RH2M, WIND))
      weather_lines <- gsub("-99.0", "  -99", weather_lines, fixed = TRUE)
      writeLines(c(wth_header, weather_lines), con = output_file)

    }, error = function(e) {
      error_message <- sprintf(
        "\n--- ERROR on task %d ---\nFailed point ID: %s\nCoords: Lat: %.3f, Lon: %.3f\nError: %s\n",
        i, point_id, latitude, longitude, conditionMessage(e))
      cat(error_message); write(error_message, file = log_file, append = TRUE)
    })
  }

  message(sprintf("\nNASA-POWER + CHIRPS processing complete. Check '%s'.\n", output_dir))
}
