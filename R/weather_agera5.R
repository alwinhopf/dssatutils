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
AGERA5_TIMESERIES_PAD_DEG <- 0.0
AGERA5_GRID_DEG <- 0.1

.agera5_data_files <- function(nc_path, zip_path, unzip_dir) {
  if (file.exists(nc_path) && file.info(nc_path)$size > 0) return(nc_path)

  if (file.exists(zip_path)) {
    listing <- try(utils::unzip(zip_path, list = TRUE), silent = TRUE)
    valid <- !inherits(listing, "try-error") && nrow(listing) > 0 &&
      any(grepl("\\.nc$", listing$Name, ignore.case = TRUE) & listing$Length > 0)
    if (!valid) return(character())
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
  area_tag <- paste(vapply(area, .agera5_slug_float, character(1)), collapse = "_")
  tag <- sprintf("%s_%s_%d_%s", spec$var, ifelse(is.na(spec$sel), "na", spec$sel), yr, area_tag)
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

  partial <- paste0(job$zip_dest, ".partial")
  if (file.exists(partial)) unlink(partial)
  req <- list(dataset_short_name = "sis-agrometeorological-indicators",
              variable = job$spec$var, year = as.character(job$yr),
              month = sprintf("%02d", 1:12), day = sprintf("%02d", 1:31),
              area = job$area, version = "2_0",
              target = basename(partial))
  if (!is.na(job$spec$sel_kind)) req[[job$spec$sel_kind]] <- job$spec$sel

  err <- NULL
  tryCatch(
    ecmwfr::wf_request(request = req, path = job$cache_dir),
    error = function(e) err <<- conditionMessage(e)
  )
  partial_files <- .agera5_data_files("", partial, paste0(job$unzip_dir, ".partial"))
  if (length(partial_files)) {
    if (!file.rename(partial, job$zip_dest)) {
      err <- "download completed but atomic cache rename failed"
    }
    unlink(paste0(job$unzip_dir, ".partial"), recursive = TRUE)
  }
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
  lats <- as.numeric(lats)
  lons <- as.numeric(lons)
  chunk_degrees <- as.numeric(chunk_degrees)
  pad <- as.numeric(pad)
  if (!length(lats) || length(lats) != length(lons) ||
      any(!is.finite(lats)) || any(!is.finite(lons)) ||
      any(lats < -90 | lats > 90) || any(lons < -180 | lons > 180)) {
    stop("AgERA5 coordinates must be finite paired latitude/longitude values.", call. = FALSE)
  }
  if (!is.finite(chunk_degrees) || chunk_degrees <= 0) {
    stop("agera5_timeseries_chunk_degrees must be a positive number.", call. = FALSE)
  }
  if (!is.finite(pad) || pad < 0) stop("AgERA5 chunk padding must be non-negative.", call. = FALSE)
  max_raw <- max(AGERA5_GRID_DEG, AGERA5_TIMESERIES_MAX_EXTENT_DEG - 2 * pad)
  cells_per_chunk <- max(1L, as.integer(round(chunk_degrees / AGERA5_GRID_DEG)))
  cells_per_chunk <- min(cells_per_chunk,
                         max(1L, as.integer(floor(max_raw / AGERA5_GRID_DEG))))

  # Snap to AgERA5's fixed global 0.1-degree lattice and use globally anchored
  # chunks. Cache paths then remain stable across point subsets and resolutions.
  lat_cell <- pmin(90, pmax(-90, round(lats / AGERA5_GRID_DEG) * AGERA5_GRID_DEG))
  lon_cell <- pmin(179.9, pmax(-180, round(lons / AGERA5_GRID_DEG) * AGERA5_GRID_DEG))
  lat_idx <- as.integer(round((lat_cell + 90) / AGERA5_GRID_DEG))
  lon_idx <- as.integer(round((lon_cell + 180) / AGERA5_GRID_DEG))
  lat_chunk <- lat_idx %/% cells_per_chunk
  lon_chunk <- lon_idx %/% cells_per_chunk
  keys <- paste(lat_chunk, lon_chunk, sep = ":")

  chunks <- list()
  for (key in sort(unique(keys))) {
    idx <- which(keys == key)
    south_idx <- min(lat_chunk[idx]) * cells_per_chunk
    west_idx <- min(lon_chunk[idx]) * cells_per_chunk
    north_idx <- min(1800L, south_idx + cells_per_chunk - 1L)
    east_idx <- min(3599L, west_idx + cells_per_chunk - 1L)
    south <- -90 + south_idx * AGERA5_GRID_DEG
    north <- -90 + north_idx * AGERA5_GRID_DEG
    west <- -180 + west_idx * AGERA5_GRID_DEG
    east <- -180 + east_idx * AGERA5_GRID_DEG
    chunks[[length(chunks) + 1L]] <- list(
      idx = idx,
      bounds = c(south = south, west = west, north = north, east = east),
      # CDS rejects a degenerate area (north == south or west == east).  Treat
      # the snapped values as grid-cell centres and request their half-cell
      # envelope.  A 0.1-degree chunk therefore remains exactly one canonical
      # AgERA5 cell while satisfying the API's area geometry.
      area = c(min(90, north + AGERA5_GRID_DEG / 2 + pad),
               max(-180, west - AGERA5_GRID_DEG / 2 - pad),
               max(-90, south - AGERA5_GRID_DEG / 2 - pad),
               min(179.9, east + AGERA5_GRID_DEG / 2 + pad))
    )
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

.agera5_acquire_lock <- function(lock_path, max_wait_sec = 120, stale_sec = NULL) {
  # Kernel locks survive long requests and are released automatically on exit.
  # Never unlink lock files: concurrent callers must lock the same inode.
  if (dir.exists(lock_path)) return(NULL) # legacy directory: do not steal it
  filelock::lock(lock_path, timeout = max_wait_sec * 1000)
}

.agera5_release_lock <- function(lock) {
  if (!is.null(lock)) filelock::unlock(lock)
}

.agera5_validate_timeseries_csv <- function(path, expected_year = NULL, bounds = NULL, area = NULL) {
  if (is.null(path) || !is.character(path) || !length(path) || !file.exists(path)) return(FALSE)
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) return(FALSE)

  df <- tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0L) return(FALSE)

  csv_cols_lower <- tolower(names(df))
  req_vars_lower <- tolower(c(
    "Solar_Radiation_Flux", "Temperature_Air_2m_Max_24h", "Temperature_Air_2m_Min_24h",
    "Precipitation_Flux", "Dew_Point_Temperature_2m_Mean_24h", "Relative_Humidity_2m_15h",
    "Wind_Speed_10m_Mean_24h"
  ))
  if (!all(req_vars_lower %in% csv_cols_lower)) return(FALSE)

  date_col <- intersect("valid_time", csv_cols_lower)
  if (!length(date_col)) return(FALSE)
  lat_col <- intersect("latitude", csv_cols_lower)
  if (!length(lat_col)) return(FALSE)
  lon_col <- intersect("longitude", csv_cols_lower)
  if (!length(lon_col)) return(FALSE)

  date_col_name <- names(df)[which(csv_cols_lower == date_col[1])[1]]
  lat_col_name <- names(df)[which(csv_cols_lower == lat_col[1])[1]]
  lon_col_name <- names(df)[which(csv_cols_lower == lon_col[1])[1]]

  lat <- suppressWarnings(as.numeric(df[[lat_col_name]]))
  lon <- suppressWarnings(as.numeric(df[[lon_col_name]]))
  if (any(!is.finite(lat) | !is.finite(lon) | abs(lat) > 90 | abs(lon) > 180)) return(FALSE)
  if (!is.null(area) && any(lat > area[1] + 1e-6 | lat < area[3] - 1e-6 |
                           lon < area[2] - 1e-6 | lon > area[4] + 1e-6)) return(FALSE)
  dates <- tryCatch(as.Date(df[[date_col_name]]), error = function(e) as.Date(NA))
  if (any(is.na(dates))) return(FALSE)

  for (rv in req_vars_lower) {
    col_name <- names(df)[which(csv_cols_lower == rv)[1]]
    num_vals <- suppressWarnings(as.numeric(df[[col_name]]))
    if (any(!is.finite(num_vals) | num_vals == -99)) return(FALSE)
  }

  if (!is.null(expected_year) && is.finite(as.numeric(expected_year))) {
    yr <- as.integer(expected_year)
    if (is.null(bounds)) bounds <- .agera5_date_bounds_for_year(yr)
    if (is.null(bounds)) return(FALSE)
    exp_start <- as.Date(bounds[1])
    exp_end <- as.Date(bounds[2])
    expected_days <- as.integer(exp_end - exp_start) + 1L

    cell_keys <- paste(df[[lat_col_name]], df[[lon_col_name]], sep = "_")
    cells <- unique(cell_keys)
    for (cell in cells) {
      idx <- which(cell_keys == cell)
      cell_dates <- sort(dates[idx])
      if (length(cell_dates) != expected_days) return(FALSE)
      if (anyDuplicated(cell_dates)) return(FALSE)
      if (cell_dates[1] != exp_start || cell_dates[length(cell_dates)] != exp_end) return(FALSE)
      if (length(cell_dates) > 1L && any(diff(cell_dates) != 1)) return(FALSE)
    }
  }

  TRUE
}

.agera5_download_timeseries_impl <- function(job) {
  data_format <- tolower(job$data_format)
  if (data_format != "csv") stop("AgERA5 time-series backend requires CSV.", call. = FALSE)
  bounds <- job$bounds
  if (is.null(bounds)) bounds <- .agera5_date_bounds_for_year(job$year)
  if (is.null(bounds)) return(NULL)
  dest <- .agera5_timeseries_cache_path(job$cache_dir, job$year, job$area, data_format)
  valid <- function(path) .agera5_validate_timeseries_csv(path, job$year, bounds, job$area)
  lock <- .agera5_acquire_lock(paste0(dest, ".lock"))
  if (is.null(lock)) {
    message("AgERA5 cache busy (including possible legacy lock): ", basename(dest))
    return(NULL)
  }
  on.exit(.agera5_release_lock(lock), add = TRUE)
  if (valid(dest)) return(dest)
  # No mutations before acquiring ownership. Keep invalid canonical data until
  # a verified replacement is ready; cache-only reads never delete evidence.
  promote <- function(source) {
    if (!valid(source)) return(FALSE)
    stage <- tempfile(".publish-", tmpdir = dirname(dest), fileext = ".csv")
    on.exit(unlink(stage), add = TRUE)
    if (!file.copy(source, stage) || !valid(stage)) return(FALSE)
    if (file.exists(dest)) {
      evidence <- tempfile(paste0(basename(dest), ".invalid-"), tmpdir = dirname(dest))
      if (!file.copy(dest, evidence)) stop("Cannot preserve invalid cache evidence")
    }
    # Same-filesystem rename; fail closed rather than copying onto a live file.
    if (!file.rename(stage, dest)) stop("Cannot atomically publish AgERA5 CSV")
    TRUE
  }
  for (candidate in c(paste0(dest, ".csv"), paste0(dest, ".partial.csv"), paste0(dest, ".partial"))) {
    if (promote(candidate)) return(dest)
  }
  if (isTRUE(job$cache_only)) return(NULL)
  .agera5_ensure_ecmwfr_key(quiet = TRUE)
  stage_dir <- tempfile(".agera5-request-", tmpdir = job$cache_dir)
  dir.create(stage_dir)
  on.exit(unlink(stage_dir, recursive = TRUE), add = TRUE)
  req <- list(dataset_short_name = "sis-agrometeorological-indicators-timeseries",
              variable = vapply(.agera5_timeseries_vars, `[[`, character(1), "var"),
              date = unname(bounds), data_format = data_format,
              area = as.numeric(job$area), target = basename(dest))
  err <- NULL
  returned <- tryCatch(dssatutils:::.provider_retry(function() ecmwfr::wf_request(request = req, path = stage_dir)),
                       error = function(e) { err <<- e; NULL })
  canonical_stage <- file.path(stage_dir, basename(dest))
  candidates <- c(canonical_stage, paste0(canonical_stage, ".csv"),
                  paste0(canonical_stage, ".partial.csv"))
  if (is.character(returned)) candidates <- unique(c(candidates, returned[!is.na(returned)]))
  for (candidate in candidates) if (promote(candidate)) return(dest)
  if (inherits(err, "dssat_connectivity_error")) stop(err)
  message(sprintf("AgERA5 time-series download failed (%d): %s", job$year,
                  if (is.null(err)) "no complete CSV returned" else conditionMessage(err)))
  NULL
}

.agera5_download_timeseries_job <- function(job) .agera5_download_timeseries_impl(job)

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

.agera5_write_wth <- function(wd, pid, lat, lon, output_dir, filename = sprintf("%s.WTH", pid)) {
  # Provider writers serialize their source values and defer physical-quality
  # decisions to the shared engine-level is_wth_valid() gate. This keeps
  # AgERA5 consistent with the other weather adapters and preserves the raw
  # daily values for diagnostics.
  tavg <- (wd$TMAX + wd$TMIN) / 2
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
  out <- file.path(output_dir, filename)
  writeLines(c(hdr, lines), con = out)
  out
}

.process_weather_agera5_timeseries <- function(shapefile, start_year, end_year, output_dir,
                                                id_col, lat_col, lon_col, n_cores, log_file,
                                                agera5_cache_dir, agera5_data_format = "csv",
                                                agera5_timeseries_chunk_degrees = AGERA5_TIMESERIES_DEFAULT_CHUNK_DEG,
                                                cache_only = FALSE) {
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
  message(sprintf("  Backend: sis-agrometeorological-indicators-timeseries (%d area chunk(s), format=%s, cache_only=%s).",
                  length(chunks), agera5_data_format, cache_only))

  effective_end <- min(as.Date(sprintf("%04d-12-31", end_year)), Sys.Date() - 10)
  jobs <- list()
  for (yr in seq.int(as.integer(start_year), end_year)) {
    if (is.null(.agera5_date_bounds_for_year(yr))) next
    for (chunk in chunks) {
      jobs[[length(jobs) + 1L]] <- list(year = yr, area = chunk$area,
                                       cache_dir = agera5_cache_dir,
                                       data_format = agera5_data_format,
                                       chunk = chunk,
                                       cache_only = isTRUE(cache_only),
                                       bounds = c(sprintf("%04d-01-01", yr), as.character(min(as.Date(sprintf("%04d-12-31", yr)), effective_end))))
    }
  }
  requested_cores <- suppressWarnings(as.integer(n_cores))
  if (is.na(requested_cores) || requested_cores < 1L) requested_cores <- 1L
  workers <- min(requested_cores, AGERA5_CDS_REQUEST_CAP, max(1L, length(jobs)))
  message(sprintf(
    "  AgERA5 time-series cache/download phase: %d year-area job(s); using %d concurrent CDS request(s) (cap=%d).",
    length(jobs), workers, AGERA5_CDS_REQUEST_CAP
  ))
  if (workers > 1L && length(jobs) > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      c(".agera5_download_timeseries_job", ".agera5_date_bounds_for_year",
        ".agera5_timeseries_cache_path", ".agera5_slug_float",
        ".agera5_timeseries_vars", ".agera5_ensure_ecmwfr_key",
        ".agera5_cds_rc_candidates", ".agera5_read_cdsapirc",
        ".dssatutils_cds_default_url", ".dssatutils_cds_rc_candidates",
        ".dssatutils_read_cdsapirc", ".dssatutils_prompt_secret",
        "setup_cds_credentials", ".dssatutils_ensure_cds_credentials",
        ".agera5_download_timeseries_impl", ".agera5_validate_timeseries_csv", ".agera5_acquire_lock",
        ".agera5_release_lock"),
      envir = parent.env(environment())
    )
    parallel::clusterEvalQ(cl, {
      # Exported closures must resolve these helpers from the worker, not a
      # potentially older installed namespace (e.g. during pkgload development).
      for (name in c(".agera5_download_timeseries_job", ".agera5_download_timeseries_impl",
                     ".agera5_validate_timeseries_csv", ".agera5_timeseries_cache_path",
                     ".agera5_date_bounds_for_year", ".agera5_acquire_lock", ".agera5_release_lock")) {
        fun <- get(name, envir = .GlobalEnv)
        environment(fun) <- .GlobalEnv
        assign(name, fun, envir = .GlobalEnv)
      }
      library(ecmwfr)
      NULL
    })
    paths <- parallel::parLapply(cl, jobs, .agera5_download_timeseries_job)
  } else {
    paths <- lapply(jobs, .agera5_download_timeseries_job)
  }
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

  expected_start <- as.Date(sprintf("%04d-01-01", as.integer(start_year)))
  expected_end <- effective_end
  all_expected_dates <- seq.Date(expected_start, expected_end, by = "day")
  all_expected_codes <- sprintf("%d%03d", lubridate::year(all_expected_dates), lubridate::yday(all_expected_dates))

  written <- 0L
  for (i in seq_along(ids)) {
    pid <- ids[i]
    tryCatch({
      ps <- point_series[[pid]]
      if (is.null(ps) || !length(ps$TMAX)) {
        stop("No AgERA5 time-series data extracted for point.")
      }
      extracted_dates <- names(ps$TMAX)
      missing_dates <- setdiff(all_expected_codes, extracted_dates)
      if (length(missing_dates) > 0L) {
        stop(sprintf(
          "Incomplete time series: missing %d day(s) between %s and %s (first missing: %s, last missing: %s).",
          length(missing_dates), as.character(expected_start), as.character(expected_end),
          missing_dates[1], missing_dates[length(missing_dates)]
        ))
      }
      dates <- all_expected_codes
      get_values <- function(vname) {
        values <- ps[[vname]]
        values <- values[!duplicated(names(values), fromLast = TRUE)]
        as.numeric(values[dates])
      }
      wd <- data.frame(DATE = dates, SRAD = get_values("SRAD"), TMAX = get_values("TMAX"),
                       TMIN = get_values("TMIN"), RAIN = get_values("RAIN"), TDEW = get_values("TDEW"),
                       RH2M = get_values("RH2M"), WIND = get_values("WIND"))
      forcing <- as.matrix(wd[, names(.agera5_timeseries_vars), drop = FALSE])
      if (any(!is.finite(forcing) | forcing == -99)) stop("Incomplete required AgERA5 forcing")
      parsed <- as.Date(wd$DATE, format = "%Y%j")
      wd$YEAR <- lubridate::year(parsed); wd$MM <- lubridate::month(parsed)

      staging_file <- basename(tempfile(paste0(pid, ".WTH.tmp-"), tmpdir = output_dir))
      staging_path <- file.path(output_dir, staging_file)
      final_path <- file.path(output_dir, sprintf("%s.WTH", pid))
      if (file.exists(staging_path)) unlink(staging_path)

      .agera5_write_wth(wd, pid, lats[i], lons[i], output_dir, filename = staging_file)

      if (!file.rename(staging_path, final_path)) {
        unlink(staging_path)
        stop("Failed to atomically publish weather file.")
      }
      written <- written + 1L
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\nAgERA5 time-series point %s (%0.3f,%0.3f): %s\n",
                     ids[i], lats[i], lons[i], conditionMessage(e))
      cat(msg)
      if (!is.null(log_file)) write(msg, file = log_file, append = TRUE)
      if (exists("staging_path") && file.exists(staging_path)) unlink(staging_path)
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
                                   agera5_timeseries_chunk_degrees = AGERA5_TIMESERIES_DEFAULT_CHUNK_DEG,
                                   cache_only = FALSE) {
  backend <- gsub("-", "_", tolower(if (is.null(agera5_backend)) "gridded" else agera5_backend), fixed = TRUE)
  if (!backend %in% c("gridded", "grid", "classic", "timeseries", "time_series", "ts")) {
    stop("agera5_backend must be 'gridded' or 'timeseries'.", call. = FALSE)
  }
  if (!requireNamespace("ecmwfr", quietly = TRUE))
    stop("AgERA5 needs the 'ecmwfr' package + a Copernicus CDS key. install.packages('ecmwfr')")
  if (!isTRUE(cache_only)) .agera5_ensure_ecmwfr_key()
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(agera5_cache_dir)) dir.create(agera5_cache_dir, recursive = TRUE)
  if (backend %in% c("timeseries", "time_series", "ts")) {
    return(.process_weather_agera5_timeseries(
      shapefile, start_year, end_year, output_dir, id_col, lat_col, lon_col,
      n_cores, log_file, agera5_cache_dir, agera5_data_format,
      agera5_timeseries_chunk_degrees,
      cache_only = cache_only
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
