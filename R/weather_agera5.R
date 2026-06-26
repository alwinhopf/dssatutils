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
#   2. Store your key (ecmwfr::wf_set_key) or via ~/.cdsapirc.
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

process_weather_agera5 <- function(shapefile, start_year, end_year, output_dir,
                                   id_col, lat_col, lon_col, n_cores, log_file,
                                   agera5_cache_dir) {
  if (!requireNamespace("ecmwfr", quietly = TRUE))
    stop("AgERA5 needs the 'ecmwfr' package + a Copernicus CDS key. install.packages('ecmwfr')")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(agera5_cache_dir)) dir.create(agera5_cache_dir, recursive = TRUE)

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
      c(".agera5_data_files", ".agera5_download_job"),
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
        r <- terra::rast(res$data_files)
        tvals <- terra::time(r)
        date_codes <- sprintf("%d%03d", lubridate::year(tvals), lubridate::yday(tvals))
        ex <- terra::extract(r, pts, ID = FALSE)
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
      wd[is.na(wd)] <- -99
      wd$YEAR <- as.integer(substr(wd$DATE, 1, 4))
      wd$MM <- lubridate::month(as.Date(wd$DATE, format = "%Y%j"))
      wd$TAVG <- (wd$TMAX + wd$TMIN) / 2
      tav <- mean(wd$TAVG, na.rm = TRUE)
      monthly <- wd %>% dplyr::group_by(YEAR, MM) %>% dplyr::summarise(m = mean(TAVG), .groups = "drop")
      amp <- mean((monthly %>% dplyr::group_by(YEAR) %>% dplyr::summarise(a = max(m) - min(m), .groups = "drop"))$a, na.rm = TRUE)

      hdr <- sprintf(
        "$WEATHER DATA: AgERA5 (Point ID: %s)\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  AGE5 %8.4f %8.4f   -99 %5.1f %5.1f   2.0  10.0\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
        pid, lats[i], lons[i], tav, amp)
      # Guard against values that would overflow a %6.1f field and shift every
      # downstream column (see weather_nasapower.R). Local so it is visible in
      # the parallel worker; corrupt readings become the DSSAT missing value.
      clamp_wth <- function(x) ifelse(!is.na(x) & (x >= 9999.95 | x <= -999.95), -99, x)
      lines <- with(wd, sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                                DATE, clamp_wth(SRAD), clamp_wth(TMAX), clamp_wth(TMIN),
                                clamp_wth(RAIN), clamp_wth(TDEW), clamp_wth(RH2M), clamp_wth(WIND)))
      lines <- gsub("-99.0", "  -99", lines, fixed = TRUE)
      writeLines(c(hdr, lines), con = file.path(output_dir, sprintf("%s.WTH", pid)))
      written <- written + 1
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\nAgERA5 point %s: %s\n", pid, conditionMessage(e))
      cat(msg); write(msg, file = log_file, append = TRUE)
    })
  }
  message(sprintf("\nAgERA5 processing complete: %d/%d points written to '%s'.\n",
                  written, length(ids), output_dir))
}
