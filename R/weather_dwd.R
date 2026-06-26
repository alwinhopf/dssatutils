# File: weather_dwd.R   (R twin of python/dssatutils/weather_dwd.py)
#
# Weather source: DWD (Deutscher Wetterdienst) Open Data — daily climate station
# observations for Germany -> DSSAT .WTH.
#
# WHY: for German / Central-European sites the national service's quality-
# controlled gauge network beats a global reanalysis. The daily "kl" product
# gives max/min/mean temperature, precipitation, sunshine duration, wind,
# humidity and vapour pressure directly.
#
# ACCESS (fully open, no key): https://opendata.dwd.de/ (CDC). DSSAT's daily
# SOLAR RADIATION is estimated from sunshine duration (SDK) with the
# Angstrom-Prescott relation (FAO-56); where SDK is missing, SRAD is -99.
#
# Coverage: Germany.

DWD_KL_URL <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/"
DWD_STATION_DESC <- "KL_Tageswerte_Beschreibung_Stationen.txt"


# --- Solar radiation from sunshine duration (Angstrom-Prescott / FAO-56) ------
.dwd_extraterrestrial <- function(lat_deg, doy) {
  phi <- lat_deg * pi / 180
  dr <- 1 + 0.033 * cos(2 * pi / 365 * doy)
  decl <- 0.409 * sin(2 * pi / 365 * doy - 1.39)
  ws <- acos(pmin(pmax(-tan(phi) * tan(decl), -1), 1))
  Gsc <- 0.0820
  Ra <- (24 * 60 / pi) * Gsc * dr * (ws * sin(phi) * sin(decl) +
                                     cos(phi) * cos(decl) * sin(ws))
  list(Ra = pmax(Ra, 0), ws = ws)
}

.dwd_srad_from_sunshine <- function(lat_deg, doy, sunshine_h, a_s = 0.25, b_s = 0.50) {
  e <- .dwd_extraterrestrial(lat_deg, doy)
  N <- 24 / pi * e$ws
  frac <- ifelse(N > 0, pmin(pmax(sunshine_h / N, 0), 1), 0)
  Rs <- (a_s + b_s * frac) * e$Ra
  Rs <- pmin(pmax(Rs, 0), 0.8 * e$Ra)
  ifelse(is.na(sunshine_h), NA_real_, Rs)
}

.dwd_tdew_from_vp <- function(vpm_hpa) {
  ln <- ifelse(vpm_hpa > 0, log(vpm_hpa / 6.1094), NA_real_)
  (243.04 * ln) / (17.625 - ln)
}


# --- Station metadata + per-station daily data (cached on disk) ---------------
.dwd_stations <- function(cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  local <- file.path(cache_dir, DWD_STATION_DESC)
  if (!file.exists(local) || file.size(local) == 0) {
    utils::download.file(paste0(DWD_KL_URL, "historical/", DWD_STATION_DESC),
                         local, mode = "wb", quiet = TRUE)
  }
  lines <- readLines(local, encoding = "latin1", warn = FALSE)
  rows <- list()
  for (ln in lines[-(1:2)]) {                       # skip header + dashed rule
    if (!nzchar(trimws(ln))) next
    # First 6 fields are whitespace-delimited numbers; the station name follows.
    parts <- strsplit(trimws(ln), "\\s+")[[1]]
    if (length(parts) < 6) next
    sid <- suppressWarnings(as.integer(parts[1]))
    von <- suppressWarnings(as.integer(parts[2]))
    bis <- suppressWarnings(as.integer(parts[3]))
    elev <- suppressWarnings(as.numeric(parts[4]))
    lat <- suppressWarnings(as.numeric(parts[5]))
    lon <- suppressWarnings(as.numeric(parts[6]))
    if (is.na(sid) || is.na(lat) || is.na(lon)) next
    rows[[length(rows) + 1]] <- data.frame(
      station_id = sprintf("%05d", sid), von = von, bis = bis,
      elev = elev, lat = lat, lon = lon, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

.dwd_historical_index <- function(cache_dir) {
  local <- file.path(cache_dir, "_hist_index.txt")
  if (file.exists(local) && file.size(local) > 0) {
    names <- readLines(local, warn = FALSE)
  } else {
    html <- paste(readLines(paste0(DWD_KL_URL, "historical/"), warn = FALSE), collapse = "\n")
    names <- regmatches(html, gregexpr("tageswerte_KL_[^\"]+_hist\\.zip", html))[[1]]
    writeLines(names, local)
  }
  idx <- list()
  for (n in names) {
    p <- strsplit(n, "_")[[1]]
    if (length(p) >= 3) idx[[sprintf("%05d", as.integer(p[3]))]] <- n
  }
  idx
}

.dwd_parse_product <- function(zip_path) {
  files <- utils::unzip(zip_path, list = TRUE)$Name
  prod <- files[grepl("^produkt", files)]
  if (length(prod) == 0) return(data.frame())
  con <- unz(zip_path, prod[1])
  df <- tryCatch(utils::read.csv2(con, sep = ";", stringsAsFactors = FALSE,
                                  strip.white = TRUE, na.strings = c("-999", "-999.0")),
                 error = function(e) data.frame())
  names(df) <- trimws(names(df))
  if (!"MESS_DATUM" %in% names(df)) return(data.frame())
  df$DATE <- as.Date(as.character(df$MESS_DATUM), format = "%Y%m%d")
  df[!is.na(df$DATE), ]
}

.dwd_fetch_station <- function(station_id, hist_name, cache_dir) {
  frames <- list()
  if (!is.null(hist_name)) {
    local <- file.path(cache_dir, hist_name)
    if (!file.exists(local) || file.size(local) == 0) {
      try(utils::download.file(paste0(DWD_KL_URL, "historical/", hist_name),
                               local, mode = "wb", quiet = TRUE), silent = TRUE)
    }
    if (file.exists(local) && file.size(local) > 0)
      frames[[length(frames) + 1]] <- .dwd_parse_product(local)
  }
  akt <- sprintf("tageswerte_KL_%s_akt.zip", station_id)
  local_akt <- file.path(cache_dir, akt)
  if (!file.exists(local_akt) || file.size(local_akt) == 0) {
    try(utils::download.file(paste0(DWD_KL_URL, "recent/", akt),
                             local_akt, mode = "wb", quiet = TRUE), silent = TRUE)
  }
  if (file.exists(local_akt) && file.size(local_akt) > 0)
    frames[[length(frames) + 1]] <- .dwd_parse_product(local_akt)

  frames <- Filter(function(f) is.data.frame(f) && nrow(f) > 0, frames)
  if (length(frames) == 0) return(data.frame())
  out <- do.call(rbind, lapply(frames, function(f) f[, intersect(
    c("DATE", "TXK", "TNK", "RSK", "SDK", "FM", "UPM", "VPM"), names(f)), drop = FALSE]))
  out <- out[!duplicated(out$DATE), ]
  out[order(out$DATE), ]
}


# --- Build a DSSAT-ready frame for one station over [start,end] ---------------
.dwd_build_frame <- function(daily, lat, start_year, end_year) {
  d <- daily[!is.na(daily$DATE), ]
  yr <- as.integer(format(d$DATE, "%Y"))
  d <- d[yr >= start_year & yr <= end_year, ]
  if (nrow(d) == 0) return(d)
  doy <- as.integer(format(d$DATE, "%j"))
  num <- function(col) if (col %in% names(d)) suppressWarnings(as.numeric(d[[col]])) else rep(NA_real_, nrow(d))
  out <- data.frame(
    DATE = sprintf("%d%03d", as.integer(format(d$DATE, "%Y")), doy),
    YEAR = as.integer(format(d$DATE, "%Y")),
    MM   = as.integer(format(d$DATE, "%m")),
    TMAX = num("TXK"), TMIN = num("TNK"), RAIN = num("RSK"),
    RH2M = num("UPM"), WIND = num("FM"),
    stringsAsFactors = FALSE)
  out$SRAD <- .dwd_srad_from_sunshine(lat, doy, num("SDK"))
  out$TDEW <- .dwd_tdew_from_vp(num("VPM"))
  out[!is.na(out$TMAX) & !is.na(out$TMIN), ]
}


# --- .WTH writer (shared format with the other weather modules) ---------------
.dwd_write_wth <- function(df, pid, lat, lon, elev, output_dir) {
  tavg <- (df$TMAX + df$TMIN) / 2
  tav <- mean(tavg, na.rm = TRUE)
  mon <- tapply(tavg, list(df$YEAR, df$MM), mean, na.rm = TRUE)
  amps <- apply(mon, 1, function(r) { r <- r[is.finite(r)]; if (length(r)) max(r) - min(r) else NA })
  amp <- mean(amps, na.rm = TRUE)
  elev_str <- if (!is.null(elev) && is.finite(elev)) sprintf("%5.0f", elev) else "  -99"
  header <- sprintf(
    "$WEATHER DATA: DWD (Point ID: %s)\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  DWD  %8.4f %8.4f %s %5.1f %5.1f   2.0  10.0\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
    pid, lat, lon, elev_str, tav, amp)
  clamp <- function(x) ifelse(!is.na(x) & (x >= 9999.95 | x <= -999.95), -99, x)
  d <- df; for (c in c("SRAD","TMAX","TMIN","RAIN","TDEW","RH2M","WIND")) d[[c]][is.na(d[[c]])] <- -99
  lines <- sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                   d$DATE, clamp(d$SRAD), clamp(d$TMAX), clamp(d$TMIN),
                   clamp(d$RAIN), clamp(d$TDEW), clamp(d$RH2M), clamp(d$WIND))
  lines <- gsub("-99.0", "  -99", lines, fixed = TRUE)
  writeLines(c(header, lines), file.path(output_dir, sprintf("%s.WTH", pid)))
}


# --- Public entry point -------------------------------------------------------
process_weather_dwd <- function(shapefile, start_year, end_year, output_dir,
                                id_col, lat_col, lon_col, n_cores, log_file,
                                dwd_cache_dir, max_station_km = 70) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dwd_cache_dir, recursive = TRUE, showWarnings = FALSE)
  end_year <- min(end_year, as.integer(format(Sys.Date(), "%Y")))
  message(sprintf("--- Starting DWD Download (Years: %d-%d) ---", start_year, end_year))

  stations <- .dwd_stations(dwd_cache_dir)
  hist_idx <- .dwd_historical_index(dwd_cache_dir)
  stations <- stations[stations$von %/% 10000 <= end_year & stations$bis %/% 10000 >= start_year, ]
  if (nrow(stations) == 0) { message("DWD: no stations cover the requested period."); return(invisible()) }
  message(sprintf("DWD: %d candidate stations cover %d-%d.", nrow(stations), start_year, end_year))

  pts <- sf::st_transform(shapefile, 4326)
  xy <- sf::st_coordinates(pts)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  lats <- xy[, 2]; lons <- xy[, 1]

  station_cache <- new.env(parent = emptyenv())
  get_daily <- function(sid) {
    if (is.null(station_cache[[sid]]))
      station_cache[[sid]] <- .dwd_fetch_station(sid, hist_idx[[sid]], dwd_cache_dir)
    station_cache[[sid]]
  }

  written <- 0
  for (k in seq_along(ids)) {
    pid <- ids[k]; lat <- lats[k]; lon <- lons[k]
    out_path <- file.path(output_dir, sprintf("%s.WTH", pid))
    if (file.exists(out_path)) { written <- written + 1; next }
    tryCatch({
      dlat <- (stations$lat - lat) * pi / 180
      dlon <- (stations$lon - lon) * pi / 180
      a <- sin(dlat / 2)^2 + cos(lat * pi / 180) * cos(stations$lat * pi / 180) * sin(dlon / 2)^2
      dist <- 6371 * 2 * asin(pmin(1, sqrt(a)))
      ord <- order(dist)
      frame <- NULL; used <- NULL
      for (j in head(ord, 8)) {
        if (dist[j] > max_station_km) break
        sid <- stations$station_id[j]
        daily <- get_daily(sid)
        if (!is.data.frame(daily) || nrow(daily) == 0) next
        f <- .dwd_build_frame(daily, stations$lat[j], start_year, end_year)
        if (nrow(f) >= 30 && sum(is.na(f$RAIN)) < nrow(f)) {
          frame <- f; used <- stations[j, ]; break
        }
      }
      if (is.null(frame))
        stop(sprintf("no DWD station within %.0f km with data for %d-%d", max_station_km, start_year, end_year))
      .dwd_write_wth(frame, pid, used$lat, used$lon, used$elev, output_dir)
      written <- written + 1
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\nDWD point %s (%.3f,%.3f): %s\n", pid, lat, lon, conditionMessage(e))
      cat(msg); write(msg, file = log_file, append = TRUE)
    })
  }
  message(sprintf("\nDWD processing complete: %d/%d points written to '%s'.\n", written, length(ids), output_dir))
}
