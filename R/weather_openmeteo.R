# File: weather_openmeteo.R
# ---------------------------------------------------------------------------
# Open-Meteo Historical Weather (ERA5 / ERA5-Land reanalysis) -> DSSAT .WTH.
#
# WHY THIS SOURCE: truly GLOBAL daily coverage (Europe, Asia, Africa, Oceania,
# South America) from 1940 onward, NO API KEY, no registration. Complements
# DAYMET (North America) and GRIDMET (US). A higher-resolution (~9 km, ERA5-Land)
# alternative to NASA POWER for non-US regions.
#
# API:     https://open-meteo.com/en/docs/historical-weather-api
# License: ERA5 data is CC-BY 4.0 (Copernicus/ECMWF) — cite when publishing.
#
# Same argument signature as process_weather_nasapower() — a drop-in source.
# ---------------------------------------------------------------------------

# Log-wind-profile factor: 10 m -> 2 m wind (FAO-56), ~0.748.
.OMET_WIND_10M_TO_2M <- 0.748

process_weather_openmeteo <- function(shapefile, start_year, end_year, output_dir,
                                      id_col, lat_col, lon_col, n_cores, log_file) {

  message(sprintf("--- Starting Open-Meteo (ERA5) Download (Years: %d-%d) ---",
                  start_year, end_year))

  start_date_str <- paste0(start_year, "-01-01")
  current_year <- lubridate::year(Sys.Date())
  if (end_year >= current_year) {
    end_date_str <- as.character(Sys.Date() - 6)  # ERA5 archive lags ~5 days
    message(sprintf("End year is current/future. Fetching data up to: %s", end_date_str))
  } else {
    end_date_str <- paste0(end_year, "-12-31")
  }

  daily_vars <- paste(c("temperature_2m_max", "temperature_2m_min",
                        "precipitation_sum", "shortwave_radiation_sum",
                        "wind_speed_10m_max"), collapse = ",")

  if (n_cores > 1) {
    message(sprintf(
      "Open-Meteo is rate-limited; using 1 weather core instead of %d to avoid 429 failures.",
      n_cores
    ))
    n_cores <- 1L
  }
  request_delay <- .dssatutils_config_number(
    "weather.openmeteo.request_delay_seconds",
    20
  )
  if (is.na(request_delay) || request_delay < 0) request_delay <- 20
  minutely_sleep <- .dssatutils_config_number(
    "weather.openmeteo.minutely_limit_sleep_seconds",
    75
  )
  if (is.na(minutely_sleep) || minutely_sleep < 1) minutely_sleep <- 75
  hourly_sleep <- .dssatutils_config_number(
    "weather.openmeteo.hourly_limit_sleep_seconds",
    3700
  )
  if (is.na(hourly_sleep) || hourly_sleep < 1) hourly_sleep <- 3700
  max_attempts <- as.integer(.dssatutils_config_number(
    "weather.openmeteo.max_attempts",
    12
  ))
  if (is.na(max_attempts) || max_attempts < 1L) max_attempts <- 12L
  archive_url <- .dssatutils_config_get(
    "weather.openmeteo.archive_url",
    "https://archive-api.open-meteo.com/v1/archive"
  )

  cl <- parallel::makeCluster(n_cores)
  # Safety net: release the cluster even if the download errors out before the
  # explicit stopCluster() below. try() keeps it harmless on the normal path.
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  doParallel::registerDoParallel(cl)
  message(sprintf("Registered %d cores for parallel Open-Meteo download.", n_cores))

  # Extract coordinates and IDs robustly
  coords_list <- .extract_coords(shapefile, id_col, lat_col, lon_col)
  ids <- coords_list$ids
  lats <- coords_list$lats
  lons <- coords_list$lons

  foreach(i = 1:nrow(shapefile),
          .packages = c("httr", "jsonlite", "lubridate", "dplyr"),
          .export = c(".OMET_WIND_10M_TO_2M")) %dopar% {

    latitude  <- lats[i]
    longitude <- lons[i]
    point_id  <- ids[i]
    output_file <- file.path(output_dir, sprintf("%s.WTH", point_id))

    if (file.exists(output_file)) return(NULL)

    tryCatch({
      if (request_delay > 0 && i > 1) Sys.sleep(request_delay)

      # --- Fetch with back-off that respects Open-Meteo's rate-limit messages ---
      resp <- NULL
      last_status <- NA_integer_
      last_reason <- NA_character_
      for (attempt in seq_len(max_attempts)) {
        r <- httr::GET(archive_url,
                       query = list(latitude = latitude, longitude = longitude,
                                    start_date = start_date_str, end_date = end_date_str,
                                    daily = daily_vars, windspeed_unit = "ms",
                                    timezone = "UTC"),
                       httr::timeout(180))
        last_status <- httr::status_code(r)
        body_txt <- httr::content(r, as = "text", encoding = "UTF-8")
        if (last_status == 200) {
          resp <- r
          break
        }

        last_reason <- tryCatch({
          parsed <- jsonlite::fromJSON(body_txt)
          if (!is.null(parsed$reason)) as.character(parsed$reason) else substr(body_txt, 1, 250)
        }, error = function(e) substr(body_txt, 1, 250))

        wait_seconds <- 10 * attempt
        if (last_status == 429 && grepl("hour", last_reason, ignore.case = TRUE)) {
          wait_seconds <- hourly_sleep
        } else if (last_status == 429 && grepl("minute", last_reason, ignore.case = TRUE)) {
          wait_seconds <- minutely_sleep
        }
        retry_message <- sprintf(
          "Open-Meteo %s attempt %d/%d returned HTTP %s: %s; waiting %.0f sec",
          point_id, attempt, max_attempts, last_status, last_reason, wait_seconds
        )
        cat(retry_message, "\n")
        write(retry_message, file = log_file, append = TRUE)
        Sys.sleep(wait_seconds)
      }
      if (is.null(resp)) {
        stop(sprintf("Open-Meteo request failed after %d attempt(s). Last HTTP %s: %s",
                     max_attempts, last_status, last_reason))
      }

      d     <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"))
      daily <- d$daily
      if (is.null(daily) || length(daily$time) == 0) stop("No data returned from Open-Meteo.")

      dts <- as.Date(daily$time)
      weather_data <- data.frame(
        YEAR = lubridate::year(dts),
        MM   = lubridate::month(dts),
        DOY  = lubridate::yday(dts),
        SRAD = daily$shortwave_radiation_sum,
        TMAX = daily$temperature_2m_max,
        TMIN = daily$temperature_2m_min,
        RAIN = daily$precipitation_sum,
        WIND = daily$wind_speed_10m_max * .OMET_WIND_10M_TO_2M,
        TDEW = -99,   # Open-Meteo has no daily dewpoint
        RH2M = -99    # ... or daily RH
      )
      weather_data$DATE <- sprintf("%d%03d", weather_data$YEAR, weather_data$DOY)
      num_cols <- c("SRAD", "TMAX", "TMIN", "RAIN", "WIND")
      weather_data[num_cols] <- lapply(weather_data[num_cols],
                                       function(x) ifelse(is.na(x), -99, x))

      weather_data$TAVG <- (weather_data$TMAX + weather_data$TMIN) / 2
      tav <- mean(weather_data$TAVG, na.rm = TRUE)
      monthly <- weather_data %>% dplyr::group_by(YEAR, MM) %>%
        dplyr::summarise(TAVG_MON = mean(TAVG, na.rm = TRUE), .groups = "drop")
      annual <- monthly %>% dplyr::group_by(YEAR) %>%
        dplyr::summarise(AMP_YR = max(TAVG_MON) - min(TAVG_MON), .groups = "drop")
      amp <- mean(annual$AMP_YR, na.rm = TRUE)

      wth_header <- sprintf(
        "$WEATHER DATA: OPEN-METEO ERA5 (Point ID: %s)\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  OMET %8.4f %8.4f   -99 %5.1f %5.1f   2.0   2.0\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
        point_id, latitude, longitude, tav, amp)

      # Guard against values that would overflow a %6.1f field and shift every
      # downstream column (see weather_nasapower.R). Local so it is visible in
      # the parallel worker; corrupt readings become the DSSAT missing value.
      clamp_wth <- function(x) ifelse(!is.na(x) & (x >= 9999.95 | x <= -999.95), -99, x)
      weather_lines <- with(weather_data, sprintf(
        "%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
        DATE, clamp_wth(SRAD), clamp_wth(TMAX), clamp_wth(TMIN),
        clamp_wth(RAIN), clamp_wth(TDEW), clamp_wth(RH2M), clamp_wth(WIND)))
      weather_lines <- gsub("-99.0", "  -99", weather_lines, fixed = TRUE)

      writeLines(c(wth_header, weather_lines), con = output_file)

    }, error = function(e) {
      error_message <- sprintf(
        "\n--- ERROR on task %d ---\nFailed point ID: %s\nCoords: Lat %.3f, Lon %.3f\nError: %s\n",
        i, point_id, latitude, longitude, conditionMessage(e))
      cat(error_message)
      write(error_message, file = log_file, append = TRUE)
    })
  }

  parallel::stopCluster(cl)
  message(sprintf("\nOpen-Meteo processing complete. Check the '%s' directory.\n", output_dir))
}
