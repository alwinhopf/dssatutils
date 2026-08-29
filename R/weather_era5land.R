# File: R/weather_era5land.R
# Centralized ERA5-Land weather downloader and parser for dssatutils.

.log_worker_message <- function(log_file, level = "INFO", point_id = NULL, msg = "") {
  if (is.null(log_file) || !nzchar(log_file)) return(invisible(NULL))
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  id_part <- if (!is.null(point_id) && !is.na(point_id) && nzchar(as.character(point_id))) {
    sprintf(" [ID=%s]", as.character(point_id))
  } else {
    ""
  }
  line <- sprintf(
    "[%s] [%s] [WEATHER_ERA5_LAND]%s %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    level,
    id_part,
    gsub("[\r\n\t]+", " ", paste(as.character(msg), collapse = " "))
  )
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(NULL)
}

.fill_na_with_neighbor_mean <- function(x, window = 2, min_neighbors = 2, max_iter = 10) {
  if (!is.numeric(x)) return(x)

  out <- as.numeric(x)
  if (!anyNA(out)) return(out)

  n <- length(out)

  for (iter in seq_len(max_iter)) {
    missing_idx <- which(is.na(out))
    if (length(missing_idx) == 0) break

    prev_out <- out

    for (idx in missing_idx) {
      neighbor_idx <- setdiff(seq.int(max(1, idx - window), min(n, idx + window)), idx)
      neighbor_vals <- prev_out[neighbor_idx]
      neighbor_vals <- neighbor_vals[!is.na(neighbor_vals)]
      if (length(neighbor_vals) >= min_neighbors) {
        out[idx] <- mean(neighbor_vals)
      }
    }

    if (identical(which(is.na(out)), which(is.na(prev_out)))) break
  }

  out
}

.find_first_matching_column <- function(nm, patterns) {
  if (length(nm) == 0) return(NA_character_)
  for (pat in patterns) {
    idx <- grep(pat, nm, ignore.case = TRUE, perl = TRUE)
    if (length(idx) > 0) return(nm[idx[1]])
  }
  NA_character_
}

.parse_hourly_datetime <- function(df) {
  nm <- names(df)

  date_col <- .find_first_matching_column(nm, c("^date$", "date", "valid_date"))
  time_col <- .find_first_matching_column(nm, c("^time$", "hour", "valid_time"))
  dt_col <- .find_first_matching_column(nm, c("datetime", "date_time", "valid_datetime", "^timestamp$", "^time$"))

  if (!is.na(dt_col)) {
    raw <- df[[dt_col]]
    parsed <- suppressWarnings(as.POSIXct(raw, tz = "UTC"))
    if (all(is.na(parsed))) parsed <- suppressWarnings(as.POSIXct(raw, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))
    if (all(is.na(parsed))) parsed <- suppressWarnings(as.POSIXct(raw, tz = "UTC", format = "%Y-%m-%dT%H:%M:%S"))
    if (!all(is.na(parsed))) return(parsed)
  }

  if (!is.na(date_col) && !is.na(time_col)) {
    raw <- paste(df[[date_col]], df[[time_col]])
    parsed <- suppressWarnings(as.POSIXct(raw, tz = "UTC"))
    if (all(is.na(parsed))) parsed <- suppressWarnings(as.POSIXct(raw, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))
    if (all(is.na(parsed))) parsed <- suppressWarnings(as.POSIXct(raw, tz = "UTC", format = "%Y-%m-%d %H:%M"))
    if (!all(is.na(parsed))) return(parsed)
  }

  if (!is.na(date_col)) {
    parsed <- suppressWarnings(as.POSIXct(df[[date_col]], tz = "UTC"))
    if (!all(is.na(parsed))) return(parsed)
  }

  stop("Could not identify a datetime column in the ERA5-Land CSV download.")
}

.calc_rh_from_temp_dew <- function(temp_c, dew_c) {
  rh <- 100 * exp((17.625 * dew_c) / (243.04 + dew_c) - (17.625 * temp_c) / (243.04 + temp_c))
  pmin(pmax(rh, 0), 100)
}

.ensure_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

.cap_end_date_era5_land <- function(end_year, lag_days = 5) {
  requested_end <- as.Date(sprintf("%d-12-31", end_year))
  min(requested_end, Sys.Date() - lag_days)
}

.build_era5_land_request <- function(latitude,
                                     longitude,
                                     start_date,
                                     end_date,
                                     target_file,
                                     data_format = "csv",
                                     variables = c(
                                       "2m_temperature",
                                       "2m_dewpoint_temperature",
                                       "total_precipitation",
                                       "surface_solar_radiation_downwards",
                                       "10m_u_component_of_wind",
                                       "10m_v_component_of_wind"
                                     )) {
  if (!nzchar(target_file)) stop("target_file must be supplied.")
  if (!(data_format %in% c("csv", "netcdf"))) stop("data_format must be 'csv' or 'netcdf'.")

  list(
    dataset_short_name = "reanalysis-era5-land-timeseries",
    variable = variables,
    location = list(
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude)
    ),
    date = list(
      start = as.character(start_date),
      end = as.character(end_date)
    ),
    data_format = data_format,
    target = basename(target_file)
  )
}

.download_era5_land_point_csv <- function(latitude,
                                          longitude,
                                          start_date,
                                          end_date,
                                          target_file,
                                          cds_user = "ecmwfr",
                                          verbose = FALSE) {
  dir.create(dirname(target_file), recursive = TRUE, showWarnings = FALSE)
  req <- .build_era5_land_request(
    latitude = latitude,
    longitude = longitude,
    start_date = start_date,
    end_date = end_date,
    target_file = target_file,
    data_format = "csv"
  )
  ecmwfr::wf_request(
    request = req,
    user = cds_user,
    transfer = TRUE,
    path = dirname(target_file),
    verbose = verbose
  )
}

.read_era5_land_hourly_csv <- function(csv_file) {
  df <- readr::read_csv(csv_file, show_col_types = FALSE, progress = FALSE, name_repair = "minimal")
  if (nrow(df) == 0) stop(sprintf("Downloaded ERA5-Land file is empty: %s", csv_file))

  nm <- names(df)
  time_values <- .parse_hourly_datetime(df)

  temp_col <- .find_first_matching_column(nm, c("^2m_temperature$", "2m_temperature"))
  dew_col  <- .find_first_matching_column(nm, c("^2m_dewpoint_temperature$", "2m_dewpoint_temperature"))
  tp_col   <- .find_first_matching_column(nm, c("^total_precipitation$", "total_precipitation"))
  ssrd_col <- .find_first_matching_column(nm, c("^surface_solar_radiation_downwards$", "surface_solar_radiation_downwards"))
  u10_col  <- .find_first_matching_column(nm, c("^10m_u_component_of_wind$", "10m_u_component_of_wind"))
  v10_col  <- .find_first_matching_column(nm, c("^10m_v_component_of_wind$", "10m_v_component_of_wind"))

  needed <- c(temp_col, dew_col, tp_col, ssrd_col, u10_col, v10_col)
  if (any(is.na(needed))) {
    stop(sprintf(
      "Could not find all required ERA5-Land columns in %s. Found names: %s",
      csv_file,
      paste(nm, collapse = ", ")
    ))
  }

  data.frame(
    DATETIME_UTC = time_values,
    T2M_K = .ensure_numeric(df[[temp_col]]),
    DEW_K = .ensure_numeric(df[[dew_col]]),
    TP_M = .ensure_numeric(df[[tp_col]]),
    SSRD_J = .ensure_numeric(df[[ssrd_col]]),
    U10 = .ensure_numeric(df[[u10_col]]),
    V10 = .ensure_numeric(df[[v10_col]]),
    stringsAsFactors = FALSE
  )
}

.aggregate_era5_land_to_daily <- function(hourly_df,
                                          start_date,
                                          end_date,
                                          utc_offset_hours = NULL) {
  if (!is.null(utc_offset_hours) && !is.na(utc_offset_hours)) {
    hourly_df$DATETIME_LOCAL <- hourly_df$DATETIME_UTC + lubridate::hours(as.numeric(utc_offset_hours))
  } else {
    hourly_df$DATETIME_LOCAL <- hourly_df$DATETIME_UTC
  }

  hourly_df$DATE_obj <- as.Date(hourly_df$DATETIME_LOCAL)
  hourly_df$T2M_C <- hourly_df$T2M_K - 273.15
  hourly_df$DEW_C <- hourly_df$DEW_K - 273.15
  hourly_df$RAIN_MM <- pmax(hourly_df$TP_M, 0) * 1000
  hourly_df$SRAD_MJ <- pmax(hourly_df$SSRD_J, 0) / 1e6
  hourly_df$WIND_MS <- sqrt(hourly_df$U10^2 + hourly_df$V10^2)
  hourly_df$RH2M <- .calc_rh_from_temp_dew(hourly_df$T2M_C, hourly_df$DEW_C)

  # Avoid using %>% inside parallel workers if namespace isn't fully attached
  daily <- hourly_df %>%
    dplyr::group_by(DATE_obj) %>%
    dplyr::summarise(
      SRAD = sum(SRAD_MJ, na.rm = TRUE),
      TMAX = max(T2M_C, na.rm = TRUE),
      TMIN = min(T2M_C, na.rm = TRUE),
      RAIN = sum(RAIN_MM, na.rm = TRUE),
      TDEW = mean(DEW_C, na.rm = TRUE),
      RH2M = mean(RH2M, na.rm = TRUE),
      WIND = mean(WIND_MS, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(DATE_obj)

  full_calendar <- data.frame(DATE_obj = seq(as.Date(start_date), as.Date(end_date), by = "day"))
  weather_data <- full_calendar %>%
    dplyr::left_join(daily, by = "DATE_obj") %>%
    dplyr::mutate(
      YEAR = lubridate::year(DATE_obj),
      MM = lubridate::month(DATE_obj),
      DOY = lubridate::yday(DATE_obj),
      DATE = sprintf("%d%03d", YEAR, DOY)
    )

  vars_to_repair <- c("SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND")
  missing_before <- vapply(weather_data[vars_to_repair], function(x) sum(is.na(x)), integer(1))

  for (v in vars_to_repair) {
    weather_data[[v]] <- .fill_na_with_neighbor_mean(weather_data[[v]], window = 2, min_neighbors = 2)
  }

  weather_data$RAIN[!is.na(weather_data$RAIN) & weather_data$RAIN < 0] <- 0
  weather_data$RH2M[!is.na(weather_data$RH2M)] <- pmin(pmax(weather_data$RH2M[!is.na(weather_data$RH2M)], 0), 100)

  bad_temp_order_idx <- which(!is.na(weather_data$TMAX) & !is.na(weather_data$TMIN) & weather_data$TMAX < weather_data$TMIN)
  if (length(bad_temp_order_idx) > 0) {
    tmp <- weather_data$TMAX[bad_temp_order_idx]
    weather_data$TMAX[bad_temp_order_idx] <- weather_data$TMIN[bad_temp_order_idx]
    weather_data$TMIN[bad_temp_order_idx] <- tmp
  }

  missing_after <- vapply(weather_data[vars_to_repair], function(x) sum(is.na(x)), integer(1))
  attr(weather_data, "missing_before") <- missing_before
  attr(weather_data, "missing_after") <- missing_after
  weather_data
}

.write_dssat_weather_file <- function(weather_data,
                                      latitude,
                                      longitude,
                                      output_file,
                                      point_id) {
  weather_data$TAVG <- (weather_data$TMAX + weather_data$TMIN) / 2
  tav <- mean(weather_data$TAVG, na.rm = TRUE)

  monthly_temps <- weather_data %>%
    dplyr::group_by(YEAR, MM) %>%
    dplyr::summarise(TAVG_MON = mean(TAVG, na.rm = TRUE), .groups = "drop")

  annual_amps <- monthly_temps %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarise(AMP_YR = max(TAVG_MON, na.rm = TRUE) - min(TAVG_MON, na.rm = TRUE), .groups = "drop")

  amp <- mean(annual_amps$AMP_YR, na.rm = TRUE)

  weather_data_out <- weather_data %>%
    dplyr::mutate(
      dplyr::across(c(SRAD, TMAX, TMIN, RAIN, TDEW, RH2M, WIND), ~ ifelse(is.na(.), -99, .)),
      DATE = sprintf("%d%03d", YEAR, DOY)
    )

  lines <- c(
    sprintf("$WEATHER DATA: ERA5-LAND  (Point ID: %s)", point_id),
    "@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT",
    sprintf("E5LD    %8.4f %8.4f   -99 %5.1f %5.1f   -99   -99", latitude, longitude, tav, amp),
    "@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND"
  )

  body <- with(
    weather_data_out,
    sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f", DATE, SRAD, TMAX, TMIN, RAIN, TDEW, RH2M, WIND)
  )

  writeLines(c(lines, body), con = output_file)
  invisible(output_file)
}

#' Download ERA5-Land weather and write DSSAT .WTH files for all points
#'
#' @param shapefile sf object or data frame with point coordinates.
#' @param start_year Start year.
#' @param end_year End year.
#' @param output_dir Output directory for .WTH files.
#' @param id_col Point id column name.
#' @param lat_col Latitude column name.
#' @param lon_col Longitude column name.
#' @param n_cores Number of parallel workers.
#' @param log_file Log file path.
#' @param cds_user Keyring user entry for ecmwfr.
#' @param utc_offset_hours Optional fixed offset applied before daily aggregation.
#' @param cache_dir Directory for raw ERA5 downloads.
#' @param keep_raw_downloads Keep raw hourly CSVs.
#' @param availability_lag_days Conservative lag for latest available data. Defaults to 5.
#' @export
process_weather_era5_land <- function(shapefile,
                                      start_year,
                                      end_year,
                                      output_dir,
                                      id_col,
                                      lat_col,
                                      lon_col,
                                      n_cores,
                                      log_file,
                                      cds_user = "ecmwfr",
                                      utc_offset_hours = NULL,
                                      cache_dir = file.path(output_dir, "_era5_cache"),
                                      keep_raw_downloads = FALSE,
                                      availability_lag_days = 5) {
  .dssatutils_ensure_cds_credentials(user = cds_user, require_ecmwfr = TRUE)
  
  message(sprintf("--- Starting ERA5-Land Download (Years: %d-%d) ---", start_year, end_year))

  start_date_str <- paste0(start_year, "-01-01")
  end_date <- .cap_end_date_era5_land(end_year, lag_days = availability_lag_days)
  end_date_str <- as.character(end_date)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  coords_list <- .extract_coords(shapefile, id_col, lat_col, lon_col)
  ids <- coords_list$ids
  lats <- coords_list$lats
  lons <- coords_list$lons

  .process_one <- function(i) {
    latitude <- lats[i]
    longitude <- lons[i]
    point_id <- ids[i]
    output_file <- file.path(output_dir, sprintf("%s.WTH", point_id))
    raw_csv <- file.path(cache_dir, sprintf("%s_era5land_hourly.csv", point_id))

    if (file.exists(output_file)) return(NULL)

    tryCatch({
      .download_era5_land_point_csv(
        latitude = latitude,
        longitude = longitude,
        start_date = start_date_str,
        end_date = end_date_str,
        target_file = raw_csv,
        cds_user = cds_user,
        verbose = FALSE
      )

      hourly <- .read_era5_land_hourly_csv(raw_csv)
      daily <- .aggregate_era5_land_to_daily(
        hourly_df = hourly,
        start_date = start_date_str,
        end_date = end_date_str,
        utc_offset_hours = utc_offset_hours
      )

      missing_before <- attr(daily, "missing_before")
      missing_after <- attr(daily, "missing_after")
      if (any(missing_before > 0)) {
        repair_msg <- paste(sprintf("%s %d->%d", names(missing_before), missing_before, missing_after), collapse = "; ")
        .log_worker_message(log_file, "WARN", point_id,
                            sprintf("Missing ERA5-Land daily values detected and repaired where possible: %s", repair_msg))
      }

      .write_dssat_weather_file(
        weather_data = daily,
        latitude = latitude,
        longitude = longitude,
        output_file = output_file,
        point_id = point_id
      )

      if (!keep_raw_downloads && file.exists(raw_csv)) unlink(raw_csv)
      .log_worker_message(log_file, "INFO", point_id,
                          sprintf("Successfully created ERA5-Land weather file: %s", basename(output_file)))
      NULL
    }, error = function(e) {
      .log_worker_message(log_file, "ERROR", point_id, conditionMessage(e))
      NULL
    })
  }

  requested_cores <- suppressWarnings(as.integer(n_cores))
  if (is.na(requested_cores) || requested_cores < 1L) requested_cores <- 1L
  if (requested_cores == 1L || length(ids) <= 1L) {
    # Keep one-core operation in the current process. Besides avoiding PSOCK
    # startup overhead, this makes debugging and mocked/offline tests reliable.
    lapply(seq_along(ids), .process_one)
  } else {
    cl <- parallel::makeCluster(min(requested_cores, length(ids)))
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export required internal functions and package names to parallel workers.
    parallel::clusterExport(
      cl,
      varlist = c(
        ".download_era5_land_point_csv", ".build_era5_land_request",
        ".read_era5_land_hourly_csv", ".parse_hourly_datetime",
        ".find_first_matching_column", ".ensure_numeric",
        ".aggregate_era5_land_to_daily", ".fill_na_with_neighbor_mean",
        ".calc_rh_from_temp_dew", ".write_dssat_weather_file",
        ".log_worker_message", ".dssatutils_cds_default_url",
        ".dssatutils_cds_rc_candidates", ".dssatutils_read_cdsapirc",
        ".dssatutils_prompt_secret", "setup_cds_credentials",
        ".dssatutils_ensure_cds_credentials", "cds_user"
      ),
      envir = environment()
    )

    foreach::foreach(
      i = seq_along(ids),
      .packages = c("ecmwfr", "readr", "dplyr", "lubridate"),
      .export = ".process_one"
    ) %dopar% .process_one(i)
  }

  invisible(TRUE)
}
