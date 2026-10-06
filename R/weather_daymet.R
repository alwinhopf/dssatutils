# File: weather_daymet.R


# Wrap all your existing Daymet logic in a function
process_weather_daymet <- function(shapefile, start_year, end_year, output_dir, 
                                   id_col, lat_col, lon_col, n_cores, log_file) {
  
  message(sprintf("--- Starting DAYMET Download (Years: %d-%d) ---", start_year, end_year))
  
  # Check for current year, which Daymet does NOT support
  current_year <- lubridate::year(Sys.Date())
  if(end_year >= current_year) {
    message(sprintf("WARNING: Daymet data is not available for the current year (%d).", current_year))
    message(sprintf("Adjusting end year to %d.", current_year - 1))
    end_year <- current_year - 1
  }
  
  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores < 1) {
    stop("n_cores must be a positive integer.")
  }

  if (n_cores == 1) {
    # Keep single-core runs in the current R process so local mocks and test
    # fixtures are visible. This also avoids unnecessary worker startup.
    foreach::registerDoSEQ()
    message("Registered sequential backend for Daymet download.")
  } else {
    cl <- parallel::makeCluster(n_cores)
    # Safety net: release the cluster even if the download errors out before the
    # explicit stopCluster() below. try() keeps it harmless on the normal path.
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    doParallel::registerDoParallel(cl)
    message(sprintf("Registered %d cores for parallel Daymet download.", n_cores))
  }
  
  # Extract coordinates and IDs robustly
  coords_list <- .extract_coords(shapefile, id_col, lat_col, lon_col)
  ids <- coords_list$ids
  lats <- coords_list$lats
  lons <- coords_list$lons
  
  leap_years <- seq(start_year, end_year)[lubridate::leap_year(seq(start_year, end_year))]
  
  # This is your exact code from STEP 2, just inside a function
  foreach(i = 1:nrow(shapefile), .export = c(".format_wth_value", ".weather_wind_run"), .packages = c("daymetr", "lubridate")) %dopar% {
    
    latitude <- lats[i]
    longitude <- lons[i]
    point_id <- ids[i]
    output_file <- file.path(output_dir, sprintf("%s.WTH", point_id))
    
    if (file.exists(output_file)) {
      return(NULL)
    }
    
    tryCatch({
      daymet_data <- dssatutils:::.provider_retry(function() daymetr::download_daymet(
        lat = latitude,
        lon = longitude,
        start = start_year,
        end = end_year,
        internal = TRUE,
        silent = TRUE
      ))
      
      weather_data <- daymet_data$data
      
      # Calculate TAV and AMP
      weather_data$t_avg <- (weather_data$`tmax..deg.c.` + weather_data$`tmin..deg.c.`) / 2
      tav <- mean(weather_data$t_avg)
      
      weather_data$month <- lubridate::month(as.Date(weather_data$yday - 1, origin = paste0(weather_data$year, "-01-01")))
      monthly_temps <- aggregate(t_avg ~ year + month, data = weather_data, FUN = mean)
      annual_amps <- aggregate(t_avg ~ year, data = monthly_temps, FUN = function(x) max(x) - min(x))
      amp <- mean(annual_amps$t_avg)
      
      # Calculate other weather variables
      weather_data$srad_mj <- (weather_data$`srad..W.m.2.` * weather_data$`dayl..s.`) / 1000000
      vp_Pa <- weather_data$`vp..Pa.`
      weather_data$t_dew <- (237.3 * log(vp_Pa / 611.2)) / (17.27 - log(vp_Pa / 611.2))
      es_Pa <- 611.2 * exp((17.67 * weather_data$t_avg) / (weather_data$t_avg + 243.5))
      weather_data$rh_2m <- 100 * (vp_Pa / es_Pa)
      weather_data$rh_2m[weather_data$rh_2m > 100] <- 100
      weather_data$wind <- -99
      weather_data$DATE <- sprintf("%d%03d", weather_data$year, weather_data$yday)
      
      # Handle Leap Years
      for (year in leap_years) {
        if(year %in% weather_data$year) {
          leap_day <- weather_data[weather_data$year == year & weather_data$yday == 365, ]
          if (nrow(leap_day) > 0) {
            leap_day$yday <- 366
            leap_day$DATE <- sprintf("%d366", year)
            weather_data <- rbind(weather_data, leap_day)
          }
        }
      }
      
      weather_data <- weather_data[order(weather_data$year, weather_data$yday), ]
      
      # Format and write the .WTH file
      wth_header <- sprintf(
        "$WEATHER DATA: DayMet Data (Point ID: %s)\n@ INSI      LAT      LONG  ELEV   TAV   AMP REFHT WNDHT\n DMET  %8.4f %9.4f   -99 %5.1f %5.1f   -99   -99\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
        point_id, latitude, longitude, tav, amp
      )

      # Shared formatting preserves the DSSAT separator and checks width.

      weather_lines <- with(weather_data, {
        sprintf(
          "%7s%s%s%s%s%s%s%s",
          DATE, .format_wth_value(srad_mj), .format_wth_value(`tmax..deg.c.`), .format_wth_value(`tmin..deg.c.`),
          .format_wth_value(`prcp..mm.day.`), .format_wth_value(t_dew), .format_wth_value(rh_2m), .format_wth_value(wind))
      })
      
      weather_lines <- gsub("-99.0", "  -99", weather_lines, fixed = TRUE)
      wth_content <- c(wth_header, weather_lines)
      writeLines(wth_content, con = output_file)
      
    },  
    error = function(e) {
      error_message <- sprintf(
        "\n--- ERROR on task %d ---\nFailed to process point ID: %s\nCoords: Lat: %.3f, Lon: %.3f\nOriginal error: %s\n",
        i, point_id, latitude, longitude, conditionMessage(e)
      )
      cat(error_message)
      write(error_message, file = log_file, append = TRUE)
    })
  } # End foreach loop
  
  if (exists("cl", inherits = FALSE)) {
    parallel::stopCluster(cl)
  }
  message(sprintf("\nDaymet processing complete. Check the '%s' directory.\n", output_dir))
}
