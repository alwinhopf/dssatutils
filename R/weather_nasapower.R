# File: weather_nasapower.R

# Load libraries needed *only* for this task
library(nasapower)
library(lubridate)
library(foreach)
library(doParallel)
library(dplyr)
library(sf)

# Define the function with the *exact same arguments* as the Daymet one
process_weather_nasapower <- function(shapefile, start_year, end_year, output_dir,
                                      id_col, lat_col, lon_col, n_cores, log_file) {
  
  message(sprintf("--- Starting NASA-POWER Download (Years: %d-%d) ---", start_year, end_year))
  
  # --- Handle Date Ranges ---
  start_date_str <- paste0(start_year, "-01-01")
  
  # Check if end_year is the current year. If so, use today's date.
  current_year <- lubridate::year(Sys.Date())
  if (end_year == current_year) {
    # NASA-POWER data is usually available up to a day or two ago.
    end_date_str <- as.character(Sys.Date() - 2)
    message(sprintf("End year is current year. Fetching data up to: %s", end_date_str))
  } else {
    end_date_str <- paste0(end_year, "-12-31")
  }  # ← FIX 1: else block now properly closed before cluster setup
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  message(sprintf("Registered %d cores for parallel NASA-POWER download.", n_cores))
  
  # Define the DSSAT variables from NASA-POWER
  nasa_params <- c(
    "T2M_MAX",           # Max Temp (C)
    "T2M_MIN",           # Min Temp (C)
    "ALLSKY_SFC_SW_DWN", # Solar Radiation (MJ/m^2/day)
    "PRECTOTCORR",       # Precipitation (mm/day)
    "T2MDEW",            # Dewpoint Temp (C)
    "RH2M",              # Relative Humidity (%)
    "WS2M"               # Wind speed at 2 m height (m/s)
  )  # ← FIX 2: nasa_params vector now properly closed
  
  foreach(i = 1:nrow(shapefile), .packages = c("nasapower", "lubridate", "dplyr")) %dopar% {
    
    latitude  <- shapefile[[lat_col]][i]
    longitude <- shapefile[[lon_col]][i]
    point_id  <- shapefile[[id_col]][i]
    output_file <- file.path(output_dir, sprintf("%s.WTH", point_id))
    
    if (file.exists(output_file)) {
      return(NULL)
    }  # ← FIX 3: if block properly closed
    
    tryCatch({
      
      # Download NASA-POWER data
      power_data <- get_power(
        community    = "AG",
        lonlat       = c(longitude, latitude),
        pars         = nasa_params,
        dates        = c(start_date_str, end_date_str),
        temporal_api = "DAILY"
      )  # ← FIX 4: get_power() call properly closed
      
      # Check if data is empty
      if (nrow(power_data) == 0) {
        stop("No data returned from NASA-POWER.")
      }
      
      # Rename and format
      weather_data <- power_data %>%
        dplyr::rename(
          SRAD = ALLSKY_SFC_SW_DWN,
          TMAX = T2M_MAX,
          TMIN = T2M_MIN,
          RAIN = PRECTOTCORR,
          TDEW = T2MDEW,
          RH2M = RH2M,
          WIND = WS2M
        ) %>%
        dplyr::mutate(
          DATE = sprintf("%d%03d", YEAR, DOY)
        ) %>%
        # NASA-POWER uses -999 for missing values
        dplyr::mutate(across(where(is.numeric), ~ ifelse(. == -999, -99, .)))
      
      # Calculate TAV and AMP
      weather_data$TAVG <- (weather_data$TMAX + weather_data$TMIN) / 2
      tav <- mean(weather_data$TAVG, na.rm = TRUE)
      
      monthly_temps <- weather_data %>%
        dplyr::group_by(YEAR, MM) %>%
        dplyr::summarise(TAVG_MON = mean(TAVG, na.rm = TRUE), .groups = "drop")
      
      annual_amps <- monthly_temps %>%
        dplyr::group_by(YEAR) %>%
        dplyr::summarise(AMP_YR = max(TAVG_MON) - min(TAVG_MON), .groups = "drop")
      
      amp <- mean(annual_amps$AMP_YR, na.rm = TRUE)
      
      # Format and write the .WTH file
      wth_header <- sprintf(
        #"$WEATHER DATA: NASA-POWER (Point ID: %s)\n@ INSI      LAT     LONG  ELEV  TAV  AMP REFHT WNDHT\n  NASA %8.4f %8.4f  -99 %5.1f %5.1f   2.0   2.0\n@ DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
        "$WEATHER DATA: NASA-POWER (Point ID: %s)\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  NASA %8.4f %8.4f   -99 %5.1f %5.1f   2.0   2.0\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
        point_id, latitude, longitude, tav, amp
      )
      
      weather_lines <- with(weather_data, {
        sprintf(
          "%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
          DATE, SRAD, TMAX, TMIN, RAIN, TDEW, RH2M, WIND
        )
      })
      
      weather_lines <- gsub("-99.0", "  -99", weather_lines, fixed = TRUE)
      wth_content   <- c(wth_header, weather_lines)
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
    
  }  # End foreach loop
  
  stopCluster(cl)
  message(sprintf("\nNASA-POWER processing complete. Check the '%s' directory.\n", output_dir))
  
}