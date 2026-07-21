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
#   2. Store your key with ecmwfr::wf_set_key(), or provide a Python-style
#      .cdsapirc and this module will import it into ecmwfr's keyring entry.
#   3. install.packages(c("ecmwfr","terra"))
#   Dataset: "sis-agrometeorological-indicators"
#
# Requests are queued by the CDS, so first runs can be slow; downloads are
# cached under `agera5_cache_dir`. Mirrors the Python weather_agera5.py.
#
# NOTE: provided for parity; validate once against your CDS account. Requires
# ecmwfr + terra. Same .WTH format as the NASA POWER / Open-Meteo writers.
# ---------------------------------------------------------------------------

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

.agera5_timeseries_vars <- list(
  TMAX = list(var = "2m_temperature_24_hour_maximum", col = "Temperature_Air_2m_Max_24h"),
  TMIN = list(var = "2m_temperature_24_hour_minimum", col = "Temperature_Air_2m_Min_24h"),
  SRAD = list(var = "solar_radiation_flux", col = "Solar_Radiation_Flux"),
  RAIN = list(var = "precipitation_flux", col = "Precipitation_Flux"),
  TDEW = list(var = "2m_dewpoint_temperature_24_hour_mean", col = "Dew_Point_Temperature_2m_Mean_24h"),
  RH2M = list(var = "2m_relative_humidity_at_15_00", col = "Relative_Humidity_2m_15h"),
  WIND = list(var = "10m_wind_speed_24_hour_mean", col = "Wind_Speed_10m_Mean_24h")
)

AGERA5_TIMESERIES_MAX_EXTENT_DEG <- 5.0
AGERA5_TIMESERIES_DEFAULT_CHUNK_DEG <- 4.5
AGERA5_TIMESERIES_PAD_DEG <- 0.2

.agera5_data_files <- function(nc_path, zip_path, unzip_dir) {
  if (file.exists(nc_path)) return(nc_path)

  if (file.exists(zip_path)) {
    dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
    nc_files <- list.files(unzip_dir, pattern = "\\.nc$", full.names = TRUE)
    if (!length(nc_files)) {
      utils::unzip(zip_path, exdir = unzip_dir)
      nc_files <- list.files(unzip_dir, pattern = "\\.nc$", full.names = TRUE)
    }
    return(sort(nc_files))
  }

  nc_files <- list.files(unzip_dir, pattern = "\\.nc$", full.names = TRUE)
  sort(nc_files)
}

AGERA5_CDS_REQUEST_CAP <- 4L

.agera5_cds_rc_candidates <- function() {
  .dssatutils_cds_rc_candidates()
}

.agera5_read_cdsapirc <- function(paths = .agera5_cds_rc_candidates()) {
  rc <- .dssatutils_read_cdsapirc(paths)
  if (is.null(rc)) return(NULL)
  list(token = rc$key, url = rc$url, path = rc$path)
}

.agera5_ensure_ecmwfr_key <- function(user = "ecmwfr", quiet = FALSE) {
  .dssatutils_ensure_cds_credentials(
    user = user,
    prompt = interactive(),
    quiet = quiet,
    require_ecmwfr = TRUE
  )
  invisible(TRUE)
}

.agera5_job <- function(vname, yr, spec, area, agera5_cache_dir) {
  tag <- sprintf("%s_%s_%d", spec$var, ifelse(is.na(spec$sel), "na", spec$sel), yr)
  list(
    vname = vname,
    yr = yr,
    spec = spec,
    tag = tag,
    area = area,
    nc_dest = file.path(agera5_cache_dir, sprintf("agera5_%s.nc", tag)),
    zip_dest = file.path(agera5_cache_dir, sprintf("agera5_%s.zip", tag)),
    unzip_dir = file.path(agera5_cache_dir, sprintf("agera5_%s_nc", tag)),
    cache_dir = agera5_cache_dir
  )
}

.agera5_download_job <- function(job) {
  data_files <- .agera5_data_files(job$nc_dest, job$zip_dest, job$unzip_dir)
  if (length(data_files)) {
    return(list(ok = TRUE, cached = TRUE, job = job, data_files = data_files,
                message = sprintf("  AgERA5 cache hit (%s)", job$tag)))
  }

  .agera5_ensure_ecmwfr_key(quiet = TRUE)

  req <- list(dataset_short_name = "sis-agrometeorological-indicators",
              variable = job$spec$var, year = as.character(job$yr),
              month = sprintf("%02d", 1:12), day = sprintf("%02d", 1:31),
              area = job$area, version = "2_0",
              target = basename(job$zip_dest))
  if (!is.na(job$spec$sel_kind)) req[[job$spec$sel_kind]] <- job$spec$sel

  err <- NULL
  tryCatch(
    ecmwfr::wf_request(request = req, path = job$cache_dir),
    error = function(e) err <<- conditionMessage(e)
  )
  data_files <- .agera5_data_files(job$nc_dest, job$zip_dest, job$unzip_dir)
  if (length(data_files)) {
    return(list(ok = TRUE, cached = FALSE, job = job, data_files = data_files,
                message = sprintf("  AgERA5 downloaded (%s)", job$tag)))
  }

  list(ok = FALSE, cached = FALSE, job = job, data_files = character(),
       message = sprintf("  AgERA5 download failed (%s): %s",
                         job$tag, if (is.null(err)) "no data file returned" else err))
}

.agera5_date_bounds_for_year <- function(year) {
  start <- as.Date(sprintf("%d-01-01", as.integer(year)))
  end <- as.Date(sprintf("%d-12-31", as.integer(year)))
  latest_safe <- Sys.Date() - 10
  if (end > latest_safe) end <- latest_safe
  if (end < start) return(NULL)
  c(as.character(start), as.character(end))
}

.agera5_split_timeseries_chunks <- function(lats, lons,
                                            chunk_degrees = AGERA5_TIMESERIES_DEFAULT_CHUNK_DEG,
                                            pad = AGERA5_TIMESERIES_PAD_DEG) {
  chunk_degrees <- as.numeric(chunk_degrees)
  pad <- as.numeric(pad)
  if (!is.finite(chunk_degrees) || chunk_degrees <= 0) {
    stop("agera5_timeseries_chunk_degrees must be a positive number.", call. = FALSE)
  }
  if (!is.finite(pad) || pad < 0) stop("AgERA5 chunk padding must be non-negative.", call. = FALSE)
  max_raw <- max(0.1, AGERA5_TIMESERIES_MAX_EXTENT_DEG - 2 * pad)
  step <- min(max(0.1, chunk_degrees), max_raw)
  lat_min <- min(lats); lat_max <- max(lats)
  lon_min <- min(lons); lon_max <- max(lons)
  chunks <- list()
  south <- lat_min
  repeat {
    north <- min(lat_max, south + step)
    west <- lon_min
    repeat {
      east <- min(lon_max, west + step)
      idx <- which(lats >= south - 1e-12 & lats <= north + 1e-12 &
                   lons >= west - 1e-12 & lons <= east + 1e-12)
      if (length(idx)) {
        chunks[[length(chunks) + 1L]] <- list(
          idx = idx,
          bounds = c(south = south, west = west, north = north, east = east),
          area = c(min(90, north + pad), max(-180, west - pad),
                   max(-90, south - pad), min(179.9, east + pad))
        )
      }
      if (east >= lon_max) break
      west <- east
    }
    if (north >= lat_max) break
    south <- north
  }
  chunks
}

.agera5_slug_float <- function(value) {
  gsub("\\.", "p", gsub("-", "m", sprintf("%.4f", as.numeric(value)), fixed = TRUE))
}

.agera5_timeseries_cache_path <- function(cache_dir, year, area, data_format = "csv") {
  ext <- if (tolower(data_format) == "csv") "csv" else "nc"
  tag <- paste(vapply(area, .agera5_slug_float, character(1)), collapse = "_")
  file.path(cache_dir, sprintf("agera5_timeseries_%d_%s.%s", as.integer(year), tag, ext))
}

.agera5_download_timeseries_job <- function(job) {
  data_format <- tolower(job$data_format)
  if (data_format != "csv") {
    stop("AgERA5 time-series backend currently supports data_format='csv'.", call. = FALSE)
  }
  bounds <- .agera5_date_bounds_for_year(job$year)
  if (is.null(bounds)) return(NULL)
  dest <- .agera5_timeseries_cache_path(job$cache_dir, job$year, job$area, data_format)
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)

  req <- list(
    dataset_short_name = "sis-agrometeorological-indicators-timeseries",
    variable = vapply(.agera5_timeseries_vars, `[[`, character(1), "var"),
    date = unname(bounds),
    data_format = data_format,
    area = as.numeric(job$area),
    target = basename(dest)
  )
  err <- NULL
  tryCatch(
    ecmwfr::wf_request(request = req, path = job$cache_dir),
    error = function(e) err <<- conditionMessage(e)
  )
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)
  message(sprintf("  AgERA5 time-series download failed (%d, area=%s): %s",
                  job$year, paste(job$area, collapse = ","),
                  if (is.null(err)) "no data file returned" else err))
  NULL
}

.agera5_find_timeseries_column <- function(df, expected) {
  hit <- which(tolower(names(df)) == tolower(expected))
  if (!length(hit)) {
    stop(sprintf("AgERA5 time-series CSV missing expected column '%s'.", expected), call. = FALSE)
  }
  names(df)[hit[1]]
}

.agera5_read_timeseries_csv <- function(path) {
  raw <- readr::read_csv(path, show_col_types = FALSE)
  date_col <- .agera5_find_timeseries_column(raw, "valid_time")
  lat_col <- .agera5_find_timeseries_column(raw, "latitude")
  lon_col <- .agera5_find_timeseries_column(raw, "longitude")
  valid_time <- as.POSIXct(raw[[date_col]], tz = "UTC")
  if (any(is.na(valid_time))) stop("AgERA5 time-series CSV contains invalid valid_time values.", call. = FALSE)
  out <- data.frame(
    valid_time = valid_time,
    latitude = as.numeric(raw[[lat_col]]),
    longitude = as.numeric(raw[[lon_col]])
  )
  for (vname in names(.agera5_timeseries_vars)) {
    col <- .agera5_find_timeseries_column(raw, .agera5_timeseries_vars[[vname]]$col)
    values <- as.numeric(raw[[col]])
    if (vname %in% c("TMAX", "TMIN", "TDEW")) values <- values - 273.15
    if (vname == "SRAD") values <- values * 1e-6
    out[[vname]] <- values
  }
  out$DATE <- sprintf("%d%03d", lubridate::year(out$valid_time), lubridate::yday(out$valid_time))
  out
}

.agera5_write_wth <- function(wd, pid, lat, lon, output_dir) {
  temp_ok <- wd$TMAX > -90 & wd$TMIN > -90
  tavg <- ifelse(temp_ok, (wd$TMAX + wd$TMIN) / 2, NA_real_)
  tav <- mean(tavg, na.rm = TRUE)
  monthly <- stats::aggregate(tavg, list(YEAR = wd$YEAR, MM = wd$MM), mean, na.rm = TRUE)
  amp_by_year <- stats::aggregate(monthly$x, list(YEAR = monthly$YEAR),
                                  function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  amp <- mean(amp_by_year$x, na.rm = TRUE)
  if (!is.finite(tav) || !is.finite(amp)) stop("No valid AgERA5 temperature climatology for point.")

  hdr <- sprintf(
    "$WEATHER DATA: AgERA5 (Point ID: %s)\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  AGE5 %8.4f %8.4f   -99 %5.1f %5.1f   2.0  10.0\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
    pid, lat, lon, tav, amp)
  clamp_wth <- function(x) {
    x[is.na(x) | x >= 9999.95 | x <= -999.95] <- -99
    x
  }
  lines <- with(wd, sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                            DATE, clamp_wth(SRAD), clamp_wth(TMAX), clamp_wth(TMIN),
                            clamp_wth(RAIN), clamp_wth(TDEW), clamp_wth(RH2M), clamp_wth(WIND)))
  lines <- gsub("-99.0", "  -99", lines, fixed = TRUE)
  out <- file.path(output_dir, sprintf("%s.WTH", pid))
  writeLines(c(hdr, lines), con = out)
  out
}

.process_weather_agera5_timeseries <- function(shapefile, start_year, end_year, output_dir,
                                                id_col, lat_col, lon_col, n_cores, log_file,
                                                agera5_cache_dir, agera5_data_format = "csv",
                                                agera5_timeseries_chunk_degrees = AGERA5_TIMESERIES_DEFAULT_CHUNK_DEG) {
  if (tolower(agera5_data_format) != "csv") {
    stop("AgERA5 time-series backend currently supports data_format='csv'.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(agera5_cache_dir, recursive = TRUE, showWarnings = FALSE)
  coords <- .extract_coords(shapefile, id_col, lat_col, lon_col)
  ids <- coords$ids; lats <- coords$lats; lons <- coords$lons
  end_year <- min(as.integer(end_year), lubridate::year(Sys.Date()))
  chunks <- .agera5_split_timeseries_chunks(lats, lons, agera5_timeseries_chunk_degrees)
  message(sprintf("--- Starting AgERA5 Time-Series Download (Years: %d-%d) ---", start_year, end_year))
  message(sprintf("  Backend: sis-agrometeorological-indicators-timeseries (%d area chunk(s), format=%s).",
                  length(chunks), agera5_data_format))

  jobs <- list()
  for (yr in seq.int(as.integer(start_year), end_year)) {
    if (is.null(.agera5_date_bounds_for_year(yr))) next
    for (chunk in chunks) {
      jobs[[length(jobs) + 1L]] <- list(year = yr, area = chunk$area,
                                       cache_dir = agera5_cache_dir,
                                       data_format = agera5_data_format,
                                       chunk = chunk)
    }
  }
  paths <- lapply(jobs, .agera5_download_timeseries_job)
  point_series <- setNames(lapply(ids, function(id)
    setNames(vector("list", length(.agera5_timeseries_vars)), names(.agera5_timeseries_vars))), ids)

  for (k in seq_along(jobs)) {
    path <- paths[[k]]
    if (is.null(path) || !file.exists(path)) next
    tryCatch({
      df <- .agera5_read_timeseries_csv(path)
      grids <- unique(df[c("latitude", "longitude")])
      for (j in jobs[[k]]$chunk$idx) {
        dist <- (grids$latitude - lats[j])^2 + (grids$longitude - lons[j])^2
        nearest <- grids[which.min(dist), ]
        sub <- df[abs(df$latitude - nearest$latitude) < 1e-9 &
                  abs(df$longitude - nearest$longitude) < 1e-9, ]
        sub <- sub[order(sub$valid_time), ]
        for (vname in names(.agera5_timeseries_vars)) {
          good <- is.finite(sub[[vname]])
          point_series[[ids[j]]][[vname]] <- c(
            point_series[[ids[j]]][[vname]],
            setNames(sub[[vname]][good], sub$DATE[good]))
        }
      }
    }, error = function(e) {
      msg <- sprintf("  AgERA5 time-series parse failed (%s): %s", path, conditionMessage(e))
      message(msg)
      if (!is.null(log_file)) write(msg, file = log_file, append = TRUE)
    })
  }

  written <- 0L
  for (i in seq_along(ids)) {
    tryCatch({
      ps <- point_series[[ids[i]]]
      if (!length(ps$TMAX)) stop("No AgERA5 time-series data extracted for this point.")
      dates <- sort(unique(unlist(lapply(ps, names), use.names = FALSE)))
      get_values <- function(vname) {
        values <- ps[[vname]]
        values <- values[!duplicated(names(values), fromLast = TRUE)]
        as.numeric(values[dates])
      }
      wd <- data.frame(DATE = dates, SRAD = get_values("SRAD"), TMAX = get_values("TMAX"),
                       TMIN = get_values("TMIN"), RAIN = get_values("RAIN"), TDEW = get_values("TDEW"),
                       RH2M = get_values("RH2M"), WIND = get_values("WIND"))
      parsed <- as.Date(wd$DATE, format = "%Y%j")
      wd$YEAR <- lubridate::year(parsed); wd$MM <- lubridate::month(parsed)
      .agera5_write_wth(wd, ids[i], lats[i], lons[i], output_dir)
      written <- written + 1L
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\nAgERA5 time-series point %s: %s\n", ids[i], conditionMessage(e))
      cat(msg)
      if (!is.null(log_file)) write(msg, file = log_file, append = TRUE)
    })
  }
  message(sprintf("\nAgERA5 time-series processing complete: %d/%d points written to '%s'.\n",
                  written, length(ids), output_dir))
  invisible(NULL)
}

process_weather_agera5 <- function(shapefile, start_year, end_year, output_dir,
                                   id_col, lat_col, lon_col, n_cores, log_file,
                                   agera5_cache_dir, agera5_backend = "gridded",
                                   agera5_data_format = "csv",
                                   agera5_timeseries_chunk_degrees = AGERA5_TIMESERIES_DEFAULT_CHUNK_DEG) {
  backend <- gsub("-", "_", tolower(if (is.null(agera5_backend)) "gridded" else agera5_backend), fixed = TRUE)
  if (!backend %in% c("gridded", "grid", "classic", "timeseries", "time_series", "ts")) {
    stop("agera5_backend must be 'gridded' or 'timeseries'.", call. = FALSE)
  }
  if (!requireNamespace("ecmwfr", quietly = TRUE))
    stop("AgERA5 needs the 'ecmwfr' package + a Copernicus CDS key. install.packages('ecmwfr')")
  .agera5_ensure_ecmwfr_key()
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(agera5_cache_dir)) dir.create(agera5_cache_dir, recursive = TRUE)
  if (backend %in% c("timeseries", "time_series", "ts")) {
    return(.process_weather_agera5_timeseries(
      shapefile, start_year, end_year, output_dir, id_col, lat_col, lon_col,
      n_cores, log_file, agera5_cache_dir, agera5_data_format,
      agera5_timeseries_chunk_degrees
    ))
  }

  # Extract coordinates and IDs robustly
  coords_list <- .extract_coords(shapefile, id_col, lat_col, lon_col)
  ids <- coords_list$ids
  lats <- coords_list$lats
  lons <- coords_list$lons
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

  jobs <- list()
  for (yr in start_year:end_year) {
    for (vname in names(.agera5_vars)) {
      jobs[[length(jobs) + 1L]] <- .agera5_job(vname, yr, .agera5_vars[[vname]], area, agera5_cache_dir)
    }
  }

  requested_cores <- suppressWarnings(as.integer(n_cores))
  if (is.na(requested_cores) || requested_cores < 1L) requested_cores <- 1L
  workers <- min(requested_cores, AGERA5_CDS_REQUEST_CAP, length(jobs))
  missing_jobs <- vapply(jobs, function(job) {
    !length(.agera5_data_files(job$nc_dest, job$zip_dest, job$unzip_dir))
  }, logical(1))
  if (!any(missing_jobs)) workers <- 1L
  message(sprintf(
    "  AgERA5 cache/download phase: %d variable-year job(s), %d missing; using %d concurrent CDS request(s) (cap=%d).",
    length(jobs), sum(missing_jobs), workers, AGERA5_CDS_REQUEST_CAP
  ))

  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      c(".agera5_data_files", ".agera5_download_job", ".agera5_cds_rc_candidates",
        ".agera5_read_cdsapirc", ".agera5_ensure_ecmwfr_key",
        ".dssatutils_cds_default_url", ".dssatutils_cds_rc_candidates",
        ".dssatutils_read_cdsapirc", ".dssatutils_prompt_secret",
        "setup_cds_credentials", ".dssatutils_ensure_cds_credentials"),
      envir = parent.env(environment())
    )
    parallel::clusterEvalQ(cl, { library(ecmwfr); NULL })
    download_results <- parallel::parLapply(cl, jobs, .agera5_download_job)
  } else {
    download_results <- lapply(jobs, .agera5_download_job)
  }

  for (res in download_results) {
    if (!isTRUE(res$cached) || !isTRUE(res$ok)) message(res$message)
    if (!isTRUE(res$ok) && !is.null(log_file)) {
      write(res$message, file = log_file, append = TRUE)
    }
  }

  for (res in download_results) {
    if (!isTRUE(res$ok) || !length(res$data_files)) next
    vname <- res$job$vname
    tag <- res$job$tag
      tryCatch({
        r <- suppressWarnings(terra::rast(res$data_files))
        tvals <- terra::time(r)
        date_codes <- sprintf("%d%03d", lubridate::year(tvals), lubridate::yday(tvals))
        ex <- suppressWarnings(terra::extract(r, pts, ID = FALSE))
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
      wd$YEAR <- as.integer(substr(wd$DATE, 1, 4))
      wd$MM <- lubridate::month(as.Date(wd$DATE, format = "%Y%j"))
      .agera5_write_wth(wd, pid, lats[i], lons[i], output_dir)
      written <- written + 1
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\nAgERA5 point %s: %s\n", pid, conditionMessage(e))
      cat(msg)
      if (!is.null(log_file)) write(msg, file = log_file, append = TRUE)
    })
  }
  message(sprintf("\nAgERA5 processing complete: %d/%d points written to '%s'.\n",
                  written, length(ids), output_dir))
}
