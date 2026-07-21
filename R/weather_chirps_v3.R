# File: weather_chirps_v3.R
# ---------------------------------------------------------------------------
# CHIRPS v3 daily rainfall extraction and NASA POWER hybrid.
#
# extract_chirps_v3_rainfall() returns per-point named daily rainfall vectors
# keyed by DSSAT DATE (YYYYDOY). process_weather_nasapower_chirps_v3() uses
# those vectors to replace NASA POWER rainfall while keeping all other weather
# variables from NASA POWER.
# ---------------------------------------------------------------------------

.chirps_v3_config <- function(key, default) {
  .dssatutils_config_get(paste0("weather.chirps_v3.", key), default)
}

.chirps_v3_lat_limit <- function() {
  .dssatutils_config_number("weather.chirps_v3.latitude_limit", 60)
}

.chirps_v3_nodata <- function() {
  .dssatutils_config_number("weather.chirps_v3.nodata", -9999)
}

.chirps_v3_options <- function(product = NULL, stream = NULL,
                               fetch_mode = NULL, resolution = NULL) {
  product <- tolower(if (is.null(product)) .chirps_v3_config("product", "rnl") else product)
  stream <- tolower(if (is.null(stream)) .chirps_v3_config("stream", "final") else stream)
  fetch_mode <- tolower(if (is.null(fetch_mode)) .chirps_v3_config("fetch_mode", "monthly_netcdf") else fetch_mode)
  resolution <- tolower(if (is.null(resolution)) .chirps_v3_config("resolution", "p05") else resolution)
  if (!product %in% c("rnl", "sat")) stop("chirps_product must be 'rnl' or 'sat'.", call. = FALSE)
  if (!stream %in% c("final", "prelim")) stop("chirps_stream must be 'final' or 'prelim'.", call. = FALSE)
  if (!fetch_mode %in% c("monthly_netcdf", "yearly_netcdf", "gee", "remote_cog")) {
    stop("chirps_fetch_mode must be 'monthly_netcdf', 'yearly_netcdf', 'gee', or 'remote_cog'.", call. = FALSE)
  }
  if (resolution != "p05") {
    stop("CHIRPS v3 daily data is currently supported only at resolution 'p05'.", call. = FALSE)
  }
  if (stream == "prelim" && product != "sat") {
    stop("CHIRPS v3 preliminary daily data is available only for product 'sat'.", call. = FALSE)
  }
  if (stream == "prelim" && fetch_mode %in% c("monthly_netcdf", "gee", "remote_cog")) {
    stop("CHIRPS v3 preliminary daily NetCDF is currently exposed only byYear.", call. = FALSE)
  }
  list(product = product, stream = stream, fetch_mode = fetch_mode, resolution = resolution)
}

.chirps_v3_file_info <- function(year, month = NULL, product = NULL,
                                 stream = NULL, fetch_mode = NULL,
                                 resolution = NULL) {
  opt <- .chirps_v3_options(product, stream, fetch_mode, resolution)
  if (opt$fetch_mode == "monthly_netcdf") {
    if (is.null(month)) stop("month is required for monthly_netcdf.", call. = FALSE)
    fname <- sprintf("chirps-v3.0.%d.%02d.days_%s.nc", as.integer(year),
                     as.integer(month), opt$resolution)
    base_url <- .chirps_v3_config("base_url", "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/daily")
    url <- sprintf("%s/%s/%s/netcdf/byMonth/%s", base_url,
                   opt$stream, opt$product, fname)
  } else {
    fname <- sprintf("chirps-v3.0.%s.%d.days_%s.nc", opt$product,
                     as.integer(year), opt$resolution)
    base_url <- .chirps_v3_config("base_url", "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/daily")
    url <- sprintf("%s/%s/%s/netcdf/byYear/%s", base_url,
                   opt$stream, opt$product, fname)
  }
  list(fname = fname, url = url, options = opt)
}

.chirps_v3_cache_path <- function(cache_dir, fname, product, stream, fetch_mode) {
  out_dir <- file.path(cache_dir, sprintf("v3_%s_%s_%s", stream, product, fetch_mode))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(out_dir, fname)
}

.chirps_v3_valid_netcdf <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) return(FALSE)
  tryCatch({
    r <- suppressWarnings(terra::rast(path))
    if (terra::nlyr(r) < 1) return(FALSE)
    e <- terra::ext(r)
    xy <- matrix(c((e[1] + e[2]) / 2, (e[3] + e[4]) / 2), ncol = 2)
    invisible(terra::extract(r[[1]], xy))
    TRUE
  }, error = function(e) {
    message(sprintf("  Invalid/corrupt CHIRPS v3 NetCDF %s: %s",
                    basename(path), conditionMessage(e)))
    FALSE
  })
}

.get_chirps_download_timeout <- function(timeout = NULL, default_val = NULL) {
  if (!is.null(timeout) && !is.na(timeout)) {
    return(as.integer(timeout))
  }
  cfg_val <- .dssatutils_config_number("weather.chirps_v3.download_timeout_seconds", NA)
  if (!is.na(cfg_val)) return(as.integer(cfg_val))
  # Legacy fallback for older scripts.
  opt_val <- getOption("dssatutils.timeout")
  if (!is.null(opt_val)) {
    val <- suppressWarnings(as.integer(opt_val))
    if (!is.na(val)) return(val)
  }
  if (is.null(default_val) || is.na(default_val)) default_val <- 14400
  default_val
}

.download_chirps_v3_file <- function(year, month, cache_dir, product, stream,
                                     fetch_mode, resolution, timeout = NULL) {
  info <- .chirps_v3_file_info(year, month, product, stream, fetch_mode, resolution)
  dest <- .chirps_v3_cache_path(cache_dir, info$fname, product, stream, fetch_mode)
  if (file.exists(dest) && .chirps_v3_valid_netcdf(dest)) return(dest)
  if (file.exists(dest)) {
    message(sprintf("  Removing corrupt cached CHIRPS v3 file: %s", basename(dest)))
    unlink(dest)
  }

  tmp <- paste0(dest, ".part-", Sys.getpid())
  unlink(Sys.glob(paste0(dest, ".part-*")), force = TRUE)
  old_timeout <- getOption("timeout")
  chirps_timeout <- .get_chirps_download_timeout(timeout, default_val = 14400)
  options(timeout = max(chirps_timeout, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  label <- if (!is.null(month)) sprintf("%d-%02d", year, month) else as.character(year)
  message(sprintf("  Downloading CHIRPS v3 %s/%s %s (%s)...",
                  stream, product, label, resolution))
  ok <- tryCatch({
    utils::download.file(info$url, tmp, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("  CHIRPS v3 download failed: ", conditionMessage(e))
    FALSE
  })
  if (!ok || !.chirps_v3_valid_netcdf(tmp)) {
    unlink(tmp)
    return(NULL)
  }
  file.rename(tmp, dest)
  dest
}

.chirps_v3_months_for_range <- function(start_year, end_year, months = NULL) {
  today <- Sys.Date()
  if (!is.null(months)) {
    months <- unique(as.integer(months))
    if (any(is.na(months)) || any(months < 1L | months > 12L)) {
      stop("chirps_months values must be calendar months 1..12.", call. = FALSE)
    }
  }
  rows <- list()
  for (yr in seq.int(as.integer(start_year), as.integer(end_year))) {
    if (yr > lubridate::year(today)) next
    last_month <- if (yr == lubridate::year(today)) lubridate::month(today) else 12L
    for (mo in seq_len(last_month)) {
      if (!is.null(months) && !mo %in% months) next
      rows[[length(rows) + 1L]] <- c(year = yr, month = mo)
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

.extract_chirps_v3_rain_remote_cog_native <- function(shapefile, start_year, end_year, id_col, lat_col, lon_col,
                                                      chirps_cache_dir, product, stream, months) {
  coords_list <- .extract_coords(shapefile, id_col, lat_col, lon_col)
  ids <- coords_list$ids
  lats <- coords_list$lats
  lons <- coords_list$lons

  out <- setNames(vector("list", length(ids)), ids)

  `%dopar%` <- foreach::`%dopar%`

  # Auto-register parallel backend if not registered or serial
  if (!foreach::getDoParRegistered() || foreach::getDoParWorkers() <= 1) {
    system_cores <- parallel::detectCores()
    cores_to_use <- if (!is.na(system_cores)) min(8, max(1, system_cores - 1)) else 2
    cl <- parallel::makeCluster(cores_to_use)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  for (yr in start_year:end_year) {
    for (i in seq_along(ids)) {
      pid <- ids[i]; lat <- lats[i]; lon <- lons[i]
      d_start <- as.Date(sprintf("%d-01-01", yr))
      d_end <- as.Date(sprintf("%d-12-31", yr))
      dates_seq <- seq(d_start, d_end, by = "1 day")
      if (!is.null(months)) {
        dates_seq <- dates_seq[lubridate::month(dates_seq) %in% as.integer(months)]
      }
      date_codes <- sprintf("%d%03d", lubridate::year(dates_seq), lubridate::yday(dates_seq))
      month_tag <- if (is.null(months)) "all" else {
        paste0("m", paste(sprintf("%02d", sort(unique(as.integer(months)))), collapse = "-"))
      }
      cache_fn <- sprintf("cog_cache_%.5f_%.5f_%s_%s_%d_%s.csv",
                          lat, lon, product, stream, yr, month_tag)
      cache_path <- file.path(chirps_cache_dir, cache_fn)
      if (file.exists(cache_path)) {
        tryCatch({
          df_cache <- utils::read.csv(cache_path, stringsAsFactors = FALSE)
          cache_dates <- as.character(df_cache$DATE)
          if (identical(cache_dates, date_codes) && !anyDuplicated(cache_dates)) {
            v <- setNames(as.numeric(df_cache$RAIN), cache_dates)
            v <- v[!is.na(v)]
            out[[pid]] <- c(out[[pid]], v)
            next
          } else {
            message(sprintf("  Incomplete CHIRPS COG cache for point %s year %d; re-fetching.",
                            pid, yr))
          }
        }, error = function(e) {})
      }

      if (length(dates_seq) == 0) next

      base_url <- .chirps_v3_config("base_url", "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/daily")
      urls <- sprintf("/vsicurl/%s/%s/%s/cogs/%d/chirps-v3.0.%s.%d.%02d.%02d.cog",
                      base_url, stream, product, yr, product, yr, lubridate::month(dates_seq), lubridate::day(dates_seq))

      fetched <- foreach::foreach(url = urls, .combine = rbind) %dopar% {
        tryCatch({
          r <- suppressWarnings(terra::rast(url))
          c(value = as.numeric(suppressWarnings(terra::extract(r, matrix(c(lon, lat), ncol = 2))[1, 1])),
            ok = 1)
        }, error = function(e) c(value = NA_real_, ok = 0))
      }

      vals <- fetched[, "value"]
      vals[is.na(vals) | vals <= .chirps_v3_nodata()] <- NA_real_

      if (all(fetched[, "ok"] == 1)) {
        df_cache <- data.frame(DATE = date_codes, RAIN = vals)
        utils::write.csv(df_cache, cache_path, row.names = FALSE)
      } else {
        message(sprintf(
          "  Warning: %d CHIRPS COG request(s) failed for point %s, year %d; incomplete results were not cached.",
          sum(fetched[, "ok"] != 1), pid, yr))
      }

      keep <- !is.na(vals)
      if (any(keep)) {
        v <- setNames(as.numeric(vals[keep]), date_codes[keep])
        out[[pid]] <- c(out[[pid]], v)
      }
    }
  }
  out
}


.run_python_gee_extraction <- function(shapefile, start_year, end_year, id_col, lat_col, lon_col,
                                       chirps_cache_dir, product, stream, fetch_mode, gee_project = NULL) {
  temp_csv <- tempfile(fileext = ".csv")
  temp_out <- tempfile(fileext = ".json")
  on.exit(unlink(c(temp_csv, temp_out)), add = TRUE)

  df_pts <- data.frame(
    id = as.character(shapefile[[id_col]]),
    lat = as.numeric(shapefile[[lat_col]]),
    lon = as.numeric(shapefile[[lon_col]]),
    stringsAsFactors = FALSE
  )
  utils::write.csv(df_pts, temp_csv, row.names = FALSE)

  project_arg <- if (is.null(gee_project) || gee_project == "") "None" else sprintf("'%s'", gee_project)
  py_code <- sprintf(
    "import pandas as pd; from dssatutils.weather_chirps_v3 import _extract_chirps_v3_rain_cli; _extract_chirps_v3_rain_cli('%s', %d, %d, '%s', '%s', '%s', '%s', '%s', gee_project=%s)",
    gsub("\\\\", "/", temp_csv), start_year, end_year, product, stream, fetch_mode, gsub("\\\\", "/", chirps_cache_dir), gsub("\\\\", "/", temp_out), project_arg
  )

  py_exe <- Sys.which("python")
  if (!nzchar(py_exe)) py_exe <- "python"

  res <- system2(py_exe, c("-c", shQuote(py_code)), stdout = TRUE, stderr = TRUE)
  status <- attr(res, "status")
  if (!is.null(status) && status != 0) {
    stop(sprintf("Python Earth Engine/COG extraction failed with status %s:\n%s",
                 status, paste(res, collapse = "\n")), call. = FALSE)
  }

  if (!file.exists(temp_out)) {
    stop("Python extraction completed but output file was not created.", call. = FALSE)
  }

  data <- jsonlite::fromJSON(temp_out)
  out <- setNames(vector("list", length(shapefile[[id_col]])), as.character(shapefile[[id_col]]))
  for (pid in names(data)) {
    vals <- unlist(data[[pid]])
    if (length(vals) > 0) {
      out[[pid]] <- setNames(as.numeric(vals), names(vals))
    }
  }
  out
}


extract_chirps_v3_rainfall <- function(shapefile, start_year, end_year,
                                       id_col, lat_col, lon_col,
                                       chirps_cache_dir,
                                       product = NULL, stream = NULL,
                                       fetch_mode = NULL, resolution = NULL,
                                       months = NULL, gee_project = NULL,
                                       timeout = NULL) {
  opt <- .chirps_v3_options(product, stream, fetch_mode, resolution)
  if (!is.null(months)) {
    months <- sort(unique(as.integer(months)))
    if (any(is.na(months)) || any(months < 1L | months > 12L)) {
      stop("chirps_months/months values must be calendar months 1..12.", call. = FALSE)
    }
  }
  dir.create(chirps_cache_dir, showWarnings = FALSE, recursive = TRUE)
  coords_list <- .extract_coords(shapefile, id_col, lat_col, lon_col)
  ids <- coords_list$ids
  lats <- coords_list$lats
  lons <- coords_list$lons
  out <- setNames(vector("list", length(ids)), ids)
  in_band <- abs(lats) <= .chirps_v3_lat_limit()
  if (!any(in_band)) {
    message("  All points outside CHIRPS v3 coverage (|lat| > 60); using fallback rainfall.")
    return(out)
  }

  extraction_shape <- shapefile[in_band, , drop = FALSE]

  if (opt$fetch_mode == "remote_cog") {
    message("  Extracting CHIRPS v3 rainfall remotely from COGs in parallel...")
    extracted <- .extract_chirps_v3_rain_remote_cog_native(extraction_shape, start_year, end_year,
                                                     id_col, lat_col, lon_col,
                                                     chirps_cache_dir, opt$product, opt$stream, months)
    out[names(extracted)] <- extracted
    return(out)
  }

  if (opt$fetch_mode == "gee") {
    message("  Extracting CHIRPS v3 rainfall via Google Earth Engine...")
    extracted <- .run_python_gee_extraction(extraction_shape, start_year, end_year,
                                      id_col, lat_col, lon_col,
                                      chirps_cache_dir, opt$product, opt$stream, opt$fetch_mode,
                                      gee_project = gee_project)
    out[names(extracted)] <- extracted
    return(out)
  }

  files <- character()
  if (opt$fetch_mode == "monthly_netcdf") {
    month_rows <- .chirps_v3_months_for_range(start_year, end_year, months = months)
    if (!is.null(month_rows) && nrow(month_rows) > 0) {
      for (i in seq_len(nrow(month_rows))) {
        p <- .download_chirps_v3_file(month_rows[i, "year"], month_rows[i, "month"],
                                      chirps_cache_dir, opt$product, opt$stream,
                                      opt$fetch_mode, opt$resolution, timeout = timeout)
        if (!is.null(p) && file.exists(p)) files <- c(files, p)
      }
    }
  } else {
    if (!is.null(months)) {
      stop("chirps_months/months is only supported with monthly_netcdf.", call. = FALSE)
    }
    current_year <- lubridate::year(Sys.Date())
    for (yr in seq.int(as.integer(start_year), min(as.integer(end_year), current_year))) {
      p <- .download_chirps_v3_file(yr, NULL, chirps_cache_dir, opt$product,
                                    opt$stream, opt$fetch_mode, opt$resolution, timeout = timeout)
      if (!is.null(p) && file.exists(p)) files <- c(files, p)
    }
  }
  if (!length(files)) {
    message("  No CHIRPS v3 files available; using fallback rainfall.")
    return(out)
  }

  extract_ids <- ids[in_band]
  pts <- terra::vect(data.frame(lon = lons[in_band], lat = lats[in_band]),
                     geom = c("lon", "lat"), crs = "EPSG:4326")
  message(sprintf("  Extracting CHIRPS v3 %s/%s rainfall for %d point(s) from %d file(s)...",
                  opt$stream, opt$product, length(extract_ids), length(files)))
  for (path in files) {
    tryCatch({
      r <- suppressWarnings(terra::rast(path))
      tvals <- terra::time(r)
      if (all(is.na(tvals))) {
        tvals <- as.Date(gsub(".*?(\\d{4}\\.\\d{2}\\.\\d{2}).*", "\\1",
                               names(r)), format = "%Y.%m.%d")
      }
      date_codes <- sprintf("%d%03d", lubridate::year(tvals), lubridate::yday(tvals))
      ex <- as.matrix(suppressWarnings(terra::extract(r, pts, ID = FALSE)))
      ex[ex <= .chirps_v3_nodata()] <- NA
      for (j in seq_along(extract_ids)) {
        vals <- ex[j, ]
        keep <- !is.na(vals)
        if (any(keep)) {
          out[[extract_ids[j]]] <- c(out[[extract_ids[j]]], setNames(as.numeric(vals[keep]), date_codes[keep]))
        }
      }
    }, error = function(e)
      message(sprintf("  CHIRPS v3 extraction failed (%s): %s",
                      basename(path), conditionMessage(e))))
  }
  out
}

process_weather_nasapower_chirps_v3 <- function(shapefile, start_year, end_year,
                                                output_dir, id_col, lat_col, lon_col,
                                                n_cores, log_file, chirps_cache_dir,
                                                chirps_product = NULL,
                                                chirps_stream = NULL,
                                                chirps_fetch_mode = NULL,
                                                chirps_resolution = NULL,
                                                chirps_months = NULL,
                                                chirps_gee_project = NULL,
                                                timeout = NULL) {
  opt <- .chirps_v3_options(chirps_product, chirps_stream,
                            chirps_fetch_mode, chirps_resolution)
  message(sprintf("--- Starting NASA-POWER + CHIRPS-v3 %s/%s (%s, Years: %d-%d) ---",
                  opt$stream, opt$product, opt$fetch_mode, start_year, end_year))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  start_date_str <- paste0(start_year, "-01-01")
  current_year <- lubridate::year(Sys.Date())
  if (end_year == current_year) {
    end_date_str <- as.character(Sys.Date() - 2)
    message(sprintf("End year is current year. Fetching NASA POWER up to: %s", end_date_str))
  } else {
    end_date_str <- paste0(end_year, "-12-31")
  }

  coords_list <- .extract_coords(shapefile, id_col, lat_col, lon_col)
  ids <- coords_list$ids
  lats <- coords_list$lats
  lons <- coords_list$lons
  chirps_rain <- extract_chirps_v3_rainfall(
    shapefile, start_year, end_year, id_col, lat_col, lon_col, chirps_cache_dir,
    product = opt$product, stream = opt$stream,
    fetch_mode = opt$fetch_mode, resolution = opt$resolution,
    months = chirps_months, gee_project = chirps_gee_project,
    timeout = timeout
  )

  nasa_params <- c("T2M_MAX", "T2M_MIN", "ALLSKY_SFC_SW_DWN", "PRECTOTCORR",
                   "T2MDEW", "RH2M", "WS2M")

  for (i in seq_along(ids)) {
    latitude <- lats[i]; longitude <- lons[i]; point_id <- ids[i]
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
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ ifelse(. == -999, -99, .)))

      merge <- merge_rainfall_into_weather(weather_data, chirps_rain[[point_id]])
      weather_data <- merge$weather_data
      n_chirps <- if (abs(latitude) <= .chirps_v3_lat_limit()) merge$n_replaced else 0L
      rain_src <- if (n_chirps > 0)
        sprintf("CHIRPS-v3 %s/%s where available, %d days; NASA-POWER otherwise",
                opt$stream, opt$product, n_chirps)
      else "NASA-POWER (CHIRPS-v3 unavailable here)"

      weather_data <- .normalize_weather_missing_values(
        weather_data,
        point_id = point_id,
        output_file = output_file,
        output_dir = output_dir,
        source_label = "NASA_POWER_CHIRPS_V3"
      )

      weather_data$TAVG <- (weather_data$TMAX + weather_data$TMIN) / 2
      tav <- mean(weather_data$TAVG, na.rm = TRUE)
      monthly_temps <- weather_data %>% dplyr::group_by(YEAR, MM) %>%
        dplyr::summarise(TAVG_MON = mean(TAVG, na.rm = TRUE), .groups = "drop")
      annual_amps <- monthly_temps %>% dplyr::group_by(YEAR) %>%
        dplyr::summarise(AMP_YR = max(TAVG_MON) - min(TAVG_MON), .groups = "drop")
      amp <- mean(annual_amps$AMP_YR, na.rm = TRUE)

      wth_header <- sprintf(
        "$WEATHER DATA: NASA-POWER + CHIRPS-v3 rain (Point ID: %s) [%s]\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  NCV3 %8.4f %8.4f   -99 %5.1f %5.1f   2.0   2.0\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
        point_id, rain_src, latitude, longitude, tav, amp)

      clamp_wth <- function(x) ifelse(!is.na(x) & (x >= 9999.95 | x <= -999.95), -99, x)
      weather_lines <- with(weather_data, sprintf(
        "%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
        DATE, clamp_wth(SRAD), clamp_wth(TMAX), clamp_wth(TMIN),
        clamp_wth(RAIN), clamp_wth(TDEW), clamp_wth(RH2M), clamp_wth(WIND)))
      weather_lines <- gsub("-99.0", "  -99", weather_lines, fixed = TRUE)
      writeLines(c(wth_header, weather_lines), con = output_file)
    }, error = function(e) {
      error_message <- sprintf(
        "\n--- ERROR on task %d ---\nFailed point ID: %s\nCoords: Lat: %.3f, Lon: %.3f\nError: %s\n",
        i, point_id, latitude, longitude, conditionMessage(e))
      cat(error_message); write(error_message, file = log_file, append = TRUE)
    })
  }

  message(sprintf("\nNASA-POWER + CHIRPS-v3 processing complete. Check '%s'.\n", output_dir))
}
