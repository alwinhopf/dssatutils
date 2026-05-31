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

library(httr)
library(jsonlite)
library(lubridate)
library(foreach)
library(doParallel)
library(dplyr)

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

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  message(sprintf("Registered %d cores for parallel Open-Meteo download.", n_cores))

  foreach(i = 1:nrow(shapefile),
          .packages = c("httr", "jsonlite", "lubridate", "dplyr"),
          .export = c(".OMET_WIND_10M_TO_2M")) %dopar% {

    latitude  <- shapefile[[lat_col]][i]
    longitude <- shapefile[[lon_col]][i]
    point_id  <- shapefile[[id_col]][i]
    output_file <- file.path(output_dir, sprintf("%s.WTH", point_id))

    if (file.exists(output_file)) return(NULL)

    tryCatch({
      # --- Fetch with simple exponential back-off (handles 429) ---
      resp <- NULL
      for (attempt in 1:4) {
        r <- httr::GET("https://archive-api.open-meteo.com/v1/archive",
                       query = list(latitude = latitude, longitude = longitude,
                                    start_date = start_date_str, end_date = end_date_str,
                                    daily = daily_vars, windspeed_unit = "ms",
                                    timezone = "UTC"),
                       httr::timeout(180))
        if (httr::status_code(r) == 200) { resp <- r; break }
        Sys.sleep(5 * attempt)
      }
      if (is.null(resp)) stop("Open-Meteo request failed after retries.")

      d <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"))
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

      weather_lines <- with(weather_data, sprintf(
        "%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
        DATE, SRAD, TMAX, TMIN, RAIN, TDEW, RH2M, WIND))
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

  stopCluster(cl)
  message(sprintf("\nOpen-Meteo processing complete. Check the '%s' directory.\n", output_dir))
}
