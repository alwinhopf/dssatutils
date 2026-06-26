# Weather source: GHCN-Daily station observations (NOAA, live download).
# Snaps each grid point to the nearest station with Tmax/Tmin over the period and
# writes a DSSAT .WTH. SRAD/RH/wind are not core GHCN elements -> written -99.

GHCN_STATIONS_URL <- "https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"
# NCEI data service: server-side filter by station/period/element, metric units.
GHCN_DATA_SERVICE <- "https://www.ncei.noaa.gov/access/services/data/v1"
# Descriptive UA (good-citizen practice). NCEI also rate-limits by IP, so heavy
# runs should reuse the cache rather than re-query in tight loops.
GHCN_USER_AGENT <- "dssatutils/0.4 (DSSAT weather pipeline; research use)"

.ghcn_load_stations <- function(cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(cache_dir, "ghcnd-stations.txt")
  if (!file.exists(path))
    utils::download.file(GHCN_STATIONS_URL, path, mode = "wb", quiet = TRUE)
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nchar(ln) >= 31]
  data.frame(
    sid = trimws(substr(ln, 1, 11)),
    lat = suppressWarnings(as.numeric(substr(ln, 13, 20))),
    lon = suppressWarnings(as.numeric(substr(ln, 22, 30))),
    stringsAsFactors = FALSE)
}

.ghcn_haversine <- function(lat0, lon0, lats, lons) {
  R <- 6371
  dlat <- (lats - lat0) * pi / 180
  dlon <- (lons - lon0) * pi / 180
  a <- sin(dlat / 2)^2 + cos(lat0 * pi / 180) * cos(lats * pi / 180) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

.ghcn_fetch_station <- function(sid, start_year, end_year) {
  # NCEI data service with metric units: Tmax/Tmin in degC, precip in mm (no /10).
  url <- sprintf(paste0(GHCN_DATA_SERVICE,
                        "?dataset=daily-summaries&stations=%s&startDate=%d-01-01",
                        "&endDate=%d-12-31&dataTypes=TMAX,TMIN,PRCP&units=metric&format=csv"),
                 sid, start_year, end_year)
  df <- tryCatch(utils::read.csv(url, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !all(c("DATE", "TMAX", "TMIN") %in% names(df)) || !nrow(df)) return(NULL)
  d <- as.Date(df$DATE)
  out <- data.frame(
    DATE = sprintf("%d%03d", as.integer(format(d, "%Y")), as.integer(format(d, "%j"))),
    YEAR = as.integer(format(d, "%Y")), MM = as.integer(format(d, "%m")),
    SRAD = -99,
    TMAX = suppressWarnings(as.numeric(df$TMAX)),
    TMIN = suppressWarnings(as.numeric(df$TMIN)),
    RAIN = if ("PRCP" %in% names(df)) suppressWarnings(as.numeric(df$PRCP)) else -99,
    TDEW = -99, RH2M = -99, WIND = -99, stringsAsFactors = FALSE)
  out[!is.na(out$TMAX) & !is.na(out$TMIN), ]
}

process_weather_ghcn <- function(shapefile, start_year, end_year, output_dir,
                                 id_col, lat_col, lon_col, n_cores, log_file,
                                 ghcn_cache_dir, max_candidates = 8) {
  old_ua <- getOption("HTTPUserAgent")
  options(HTTPUserAgent = GHCN_USER_AGENT)
  on.exit(options(HTTPUserAgent = old_ua), add = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stations <- .ghcn_load_stations(ghcn_cache_dir)
  pts <- sf::st_transform(shapefile, 4326)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  xy <- sf::st_coordinates(pts); lats <- xy[, 2]; lons <- xy[, 1]
  message(sprintf("--- Starting GHCN-Daily Processing (Years: %d-%d) ---", start_year, end_year))
  written <- 0
  for (k in seq_along(ids)) {
    pid <- ids[k]; lat <- lats[k]; lon <- lons[k]
    tryCatch({
      ord <- order(.ghcn_haversine(lat, lon, stations$lat, stations$lon))
      frame <- NULL; chosen <- NA_character_
      for (j in ord[seq_len(min(max_candidates, length(ord)))]) {
        frame <- .ghcn_fetch_station(stations$sid[j], start_year, end_year)
        if (!is.null(frame) && nrow(frame)) { chosen <- stations$sid[j]; break }
      }
      if (is.null(frame) || !nrow(frame))
        stop(sprintf("no GHCN station with Tmax/Tmin near (%.3f,%.3f)", lat, lon))
      weather_write_wth(frame, pid, lat, lon, output_dir,
                        sprintf("GHCN-Daily %s", chosen), "GHCN")
      written <- written + 1
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\nGHCN point %s (%.3f,%.3f): %s\n", pid, lat, lon, conditionMessage(e))
      cat(msg); write(msg, file = log_file, append = TRUE)
    })
  }
  message(sprintf("\nGHCN-Daily processing complete: %d/%d point(s) written.\n", written, length(ids)))
}
