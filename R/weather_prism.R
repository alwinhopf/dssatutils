PRISM_URL <- "https://services.nacse.org/prism/data/get/us/4km/%s/%s"
PRISM_VARS <- c(ppt = "RAIN", tmax = "TMAX", tmin = "TMIN", tdmean = "TDEW")
# Polite spacing between NACSE requests (seconds) to avoid throttle responses.
PRISM_REQUEST_DELAY <- 1.0

# ACIS Web Services GridData endpoint and PRISM grid code (grid 21)
ACIS_GRIDDATA_URL <- "https://data.rcc-acis.org/GridData"
ACIS_PRISM_GRID <- "21"

.prism_is_zip <- function(path) {
  # A valid zip starts with the magic bytes "PK" (0x50 0x4B). NACSE throttle
  # responses are small HTML/text bodies that fail this check.
  if (!file.exists(path) || file.info(path)$size < 4) return(FALSE)
  magic <- readBin(path, "raw", n = 2)
  identical(magic, as.raw(c(0x50, 0x4B)))
}

.prism_download_grid <- function(var, day, cache_dir) {
  ymd <- format(day, "%Y%m%d")
  out_dir <- file.path(cache_dir, var, ymd)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  existing <- list.files(out_dir, pattern = "\\.(bil|tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
  if (length(existing)) return(existing[1])
  zip_path <- file.path(out_dir, paste0(var, "_", ymd, ".zip"))
  ok <- tryCatch({
    Sys.sleep(PRISM_REQUEST_DELAY)
    utils::download.file(sprintf(PRISM_URL, var, ymd), zip_path, mode = "wb", quiet = TRUE)
    if (.prism_is_zip(zip_path)) {
      utils::unzip(zip_path, exdir = out_dir)
      TRUE
    } else {
      message(sprintf("  PRISM %s %s: response was not a valid zip (likely throttled); skipping.", var, ymd))
      FALSE
    }
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) return(NA_character_)
  existing <- list.files(out_dir, pattern = "\\.(bil|tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
  if (length(existing)) existing[1] else NA_character_
}

.prism_acis_request <- function(lat, lon, sdate, edate) {
  list(version = 2L, loc = sprintf("%.4f,%.4f", lon, lat), grid = ACIS_PRISM_GRID,
       sdate = sdate, edate = edate, elements = "maxt,mint,pcpn", units = "degreeC,degreeC,mm")
}

.prism_acis_cache_file <- function(cache_dir, lat, lon, sdate, edate) {
  file.path(cache_dir, sprintf("prism_acis_v2_21_%.4f_%.4f_%s_%s.json", lon, lat, sdate, edate))
}

.prism_complete_acis_matrix <- function(mat, sdate, edate) {
  tryCatch({
    if (as.Date(sdate) > as.Date(edate)) return(FALSE)
    dates <- format(seq(as.Date(sdate), as.Date(edate), by = "day"), "%Y-%m-%d")
    if (!is.matrix(mat) || ncol(mat) != 4L || nrow(mat) != length(dates) ||
        !identical(as.character(mat[, 1]), dates)) return(FALSE)
    tmax <- suppressWarnings(as.numeric(mat[, 2]))
    tmin <- suppressWarnings(as.numeric(mat[, 3]))
    rain <- suppressWarnings(as.numeric(mat[, 4]))
    trace <- !is.na(mat[, 4]) & tolower(trimws(as.character(mat[, 4]))) == "t"
    rain[trace] <- 0
    all(is.finite(tmax) & is.finite(tmin) & is.finite(rain) &
        tmin >= -90 & tmax <= 70 & tmax >= tmin & rain >= 0 & rain <= 2000)
  }, error = function(e) FALSE)
}

.prism_publish_acis_wth <- function(df, pid, lat, lon, output_dir, source_label,
                                   sdate, edate, srad_method) {
  required <- c("TMAX", "TMIN", "RAIN", if (srad_method != "none") "SRAD")
  stage <- tempfile(".prism-", tmpdir = output_dir)
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  weather_write_wth(df, pid, lat, lon, stage, source_label, "PRSM", wndht = -99)
  path <- file.path(stage, paste0(pid, ".WTH"))
  if (!is_wth_valid(path, start_date = sdate, end_date = edate, required_columns = required)) return(FALSE)
  if (!file.rename(path, file.path(output_dir, paste0(pid, ".WTH")))) {
    stop("Failed to atomically publish PRISM weather file")
  }
  TRUE
}

.prism_fetch_acis_point <- function(lat, lon, sdate, edate, cache_file = NULL, max_retries = 3) {
  request <- .prism_acis_request(lat, lon, sdate, edate)
  if (!is.null(cache_file) && file.exists(cache_file)) {
    content_text <- paste(readLines(cache_file, warn = FALSE), collapse = "\n")
    parsed <- tryCatch(jsonlite::fromJSON(content_text), error = function(e) NULL)
    if (is.list(parsed) && identical(parsed$request, request) &&
        .prism_complete_acis_matrix(parsed$data, sdate, edate)) {
      return(parsed$data)
    }
  }

  payload <- list(
    loc = sprintf("%.4f,%.4f", lon, lat),
    sdate = sdate,
    edate = edate,
    grid = ACIS_PRISM_GRID,
    elems = list(
      list(name = "maxt", units = "degreeC"),
      list(name = "mint", units = "degreeC"),
      list(name = "pcpn", units = "mm")
    )
  )

  attempt <- 1
  body_json <- jsonlite::toJSON(payload, auto_unbox = TRUE)
  while (attempt <= max_retries) {
    res <- tryCatch({
      httr::POST(
        url = ACIS_GRIDDATA_URL,
        body = body_json,
        httr::content_type_json(),
        httr::timeout(60)
      )
    }, error = function(e) NULL)

    if (!is.null(res) && httr::status_code(res) == 200) {
      txt <- httr::content(res, "text", encoding = "UTF-8")
      if (nzchar(trimws(txt))) {
        parsed <- tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
        if (is.list(parsed) && .prism_complete_acis_matrix(parsed$data, sdate, edate)) {
          if (!is.null(cache_file)) {
            dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
            stage <- tempfile(".prism-", tmpdir = dirname(cache_file), fileext = ".json")
            tryCatch({
              writeLines(jsonlite::toJSON(list(request = request, data = parsed$data),
                                         auto_unbox = TRUE, digits = NA), stage)
              if (!file.rename(stage, cache_file)) stop("Failed to atomically publish PRISM cache")
            }, finally = unlink(stage))
          }
          return(parsed$data)
        }
      }
    }
    if (attempt < max_retries) Sys.sleep(1.0 * attempt)
    attempt <- attempt + 1
  }
  return(NULL)
}

.prism_parse_acis_matrix <- function(mat, lat, srad_method = "bristow_campbell") {
  if (is.null(mat) || !is.matrix(mat) || nrow(mat) == 0 || ncol(mat) < 4) {
    return(NULL)
  }
  dates_raw <- as.character(mat[, 1])
  d <- as.Date(dates_raw)
  valid_dates <- !is.na(d)
  if (!any(valid_dates)) return(NULL)

  dates_raw <- dates_raw[valid_dates]
  d <- d[valid_dates]
  tmax_raw <- mat[valid_dates, 2]
  tmin_raw <- mat[valid_dates, 3]
  rain_raw <- mat[valid_dates, 4]

  tmax <- suppressWarnings(as.numeric(tmax_raw))
  tmin <- suppressWarnings(as.numeric(tmin_raw))
  rain <- suppressWarnings(as.numeric(rain_raw))

  # ACIS trace precipitation marker "T" or "t" -> 0 mm
  trace_mask <- !is.na(rain_raw) & tolower(trimws(as.character(rain_raw))) == "t"
  rain[trace_mask] <- 0.0

  # ACIS missing / sentinel handling
  tmax[!is.finite(tmax) | tmax == -999] <- NA
  tmin[!is.finite(tmin) | tmin == -999] <- NA
  rain[!is.finite(rain) | rain == -999] <- NA

  valid_rows <- !is.na(tmax) & !is.na(tmin)
  if (!any(valid_rows)) return(NULL)

  d <- d[valid_rows]
  tmax <- tmax[valid_rows]
  tmin <- tmin[valid_rows]
  rain <- rain[valid_rows]

  DATE <- sprintf("%d%03d", as.integer(format(d, "%Y")), as.integer(format(d, "%j")))
  YEAR <- as.integer(format(d, "%Y"))
  MM   <- as.integer(format(d, "%m"))

  if (srad_method == "bristow_campbell") {
    SRAD <- estimate_srad_bristow_campbell(DATE, tmax, tmin, lat)
  } else {
    SRAD <- rep(-99, length(DATE))
  }

  data.frame(
    DATE = DATE,
    YEAR = YEAR,
    MM   = MM,
    SRAD = SRAD,
    TMAX = tmax,
    TMIN = tmin,
    RAIN = ifelse(is.na(rain), -99, rain),
    TDEW = rep(-99, length(DATE)),
    RH2M = rep(-99, length(DATE)),
    WIND = rep(-99, length(DATE)),
    stringsAsFactors = FALSE
  )
}

process_weather_prism <- function(shapefile, start_year, end_year, output_dir,
                                  id_col, lat_col, lon_col, n_cores, log_file,
                                  prism_cache_dir = NULL, srad_method = "bristow_campbell",
                                  backend = "acis") {
  backend_norm <- tolower(trimws(as.character(if (is.null(backend) || !length(backend)) "acis" else backend)[[1]]))
  if (!backend_norm %in% c("acis", "point", "nacse", "grid")) {
    stop(sprintf("Unknown PRISM backend: '%s'. Expected 'acis' or 'nacse'.", backend))
  }
  if (backend_norm %in% c("nacse", "grid") && (is.null(prism_cache_dir) || !nzchar(prism_cache_dir))) {
    stop("prism_cache_dir is required for the NACSE grid backend")
  }
  if (!srad_method %in% c("bristow_campbell", "none")) {
    stop(sprintf("Unknown srad_method: '%s'. Expected 'bristow_campbell' or 'none'.", srad_method))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(prism_cache_dir) && nzchar(as.character(prism_cache_dir))) {
    dir.create(prism_cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  coords <- .extract_coords(shapefile, id_col, lat_col, lon_col)
  ids <- coords$ids
  lats <- coords$lats
  lons <- coords$lons

  source_label <- if (srad_method == "bristow_campbell") {
    "PRISM 4km (SRAD estimated: Bristow-Campbell 1984)"
  } else {
    "PRISM 4km"
  }

  if (backend_norm %in% c("acis", "point")) {
    message(sprintf("--- Starting PRISM ACIS Point Processing (Years: %d-%d) ---", start_year, end_year))
    if (srad_method == "bristow_campbell") {
      message("  PRISM: Estimating solar radiation using Bristow-Campbell (1984)")
    } else {
      message("  PRISM: Solar radiation estimation disabled (srad_method='none')")
    }

    sdate <- sprintf("%04d-01-01", as.integer(start_year))
    latest_safe <- Sys.Date() - 2
    edate <- format(min(as.Date(sprintf("%04d-12-31", as.integer(end_year))), latest_safe), "%Y-%m-%d")
    if (sdate > edate) stop("PRISM start date is after the available end date")

    cache_dir <- if (!is.null(prism_cache_dir) && nzchar(as.character(prism_cache_dir))) {
      file.path(prism_cache_dir, "acis")
    } else {
      NULL
    }
    if (!is.null(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    written <- 0
    for (i in seq_along(ids)) {
      pid <- ids[i]
      lat <- lats[i]
      lon <- lons[i]
      cache_file <- if (!is.null(cache_dir)) {
        .prism_acis_cache_file(cache_dir, lat, lon, sdate, edate)
      } else {
        NULL
      }
      mat <- .prism_fetch_acis_point(lat, lon, sdate, edate, cache_file)
      if (!.prism_complete_acis_matrix(mat, sdate, edate)) {
        write(sprintf("PRISM point %s: incomplete or invalid ACIS history for %s through %s; not published",
                      pid, sdate, edate), file = log_file, append = TRUE)
        next
      }
      df <- .prism_parse_acis_matrix(mat, lat, srad_method)
      if (is.null(df) || nrow(df) == 0) {
        write(sprintf("PRISM point %s: no valid TMAX/TMIN data extracted", pid), file = log_file, append = TRUE)
        next
      }
      if (.prism_publish_acis_wth(df, pid, lat, lon, output_dir, source_label, sdate, edate, srad_method)) {
        written <- written + 1
      } else {
        write(sprintf("PRISM point %s: invalid WTH forcing; not published", pid), file = log_file, append = TRUE)
      }
    }
    message(sprintf("\nPRISM processing complete: %d/%d point(s) written.\n", written, length(ids)))
    return(invisible(written))
  }

  # Legacy NACSE CONUS grid download / raster extraction backend
  if (!requireNamespace("terra", quietly = TRUE)) stop("package 'terra' required for PRISM grid backend")
  dates <- seq(as.Date(sprintf("%d-01-01", start_year)),
               min(as.Date(sprintf("%d-12-31", end_year)), Sys.Date() - 2),
               by = "day")
  pts <- sf::st_transform(shapefile, 4326)
  pts_vect <- terra::vect(pts)
  frames <- setNames(vector("list", length(ids)), ids)
  for (i in seq_along(frames)) frames[[i]] <- list()
  message(sprintf("--- Starting PRISM Processing (Years: %d-%d) ---", start_year, end_year))
  if (srad_method == "bristow_campbell") {
    message("  PRISM: Estimating solar radiation using Bristow-Campbell (1984)")
  } else {
    message("  PRISM: Solar radiation estimation disabled (srad_method='none')")
  }
  for (di in seq_along(dates)) {
    day <- dates[di]
    day_vals <- list()
    for (var in names(PRISM_VARS)) {
      p <- .prism_download_grid(var, day, prism_cache_dir)
      if (!is.na(p)) day_vals[[PRISM_VARS[[var]]]] <- as.numeric(terra::extract(terra::rast(p), pts_vect, ID = FALSE)[, 1])
    }
    for (i in seq_along(ids)) {
      frames[[i]][[length(frames[[i]]) + 1]] <- data.frame(
        DATE = sprintf("%d%03d", as.integer(format(day, "%Y")), as.integer(format(day, "%j"))),
        YEAR = as.integer(format(day, "%Y")), MM = as.integer(format(day, "%m")),
        SRAD = -99, TMAX = if (!is.null(day_vals$TMAX)) day_vals$TMAX[i] else NA,
        TMIN = if (!is.null(day_vals$TMIN)) day_vals$TMIN[i] else NA,
        RAIN = if (!is.null(day_vals$RAIN)) day_vals$RAIN[i] else NA,
        TDEW = if (!is.null(day_vals$TDEW)) day_vals$TDEW[i] else -99,
        RH2M = -99, WIND = -99)
    }
  }
  written <- 0
  for (i in seq_along(ids)) {
    df <- do.call(rbind, frames[[i]])
    df <- df[!is.na(df$TMAX) & !is.na(df$TMIN), ]
    if (!nrow(df)) {
      write(sprintf("PRISM point %s: no valid TMAX/TMIN data extracted", ids[i]), file = log_file, append = TRUE)
      next
    }
    if (srad_method == "bristow_campbell") {
      df$SRAD <- estimate_srad_bristow_campbell(df$DATE, df$TMAX, df$TMIN, lats[i])
    } else {
      df$SRAD <- -99
    }
    weather_write_wth(df, ids[i], lats[i], lons[i], output_dir, source_label, "PRSM", wndht = -99)
    written <- written + 1
  }
  message(sprintf("\nPRISM processing complete: %d/%d point(s) written.\n", written, length(ids)))
}
