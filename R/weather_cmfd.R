# File: weather_cmfd.R   (R twin of python/dssatutils/weather_cmfd.py)
#
# Weather source: CMFD (China Meteorological Forcing Dataset) -> DSSAT .WTH.
#
# WHY: CMFD is the standard 0.1° near-surface forcing dataset for CHINA
# (1979-2018). It is 3-HOURLY, so this module aggregates each variable to the
# daily statistics DSSAT needs.
#
# CMFD variables (one NetCDF per variable per month; var name == token):
#   temp (K)  prec (mm/hr)  srad (W/m^2)  shum (kg/kg)  pres (Pa)  wind (m/s)
#
# Daily reduction:  TMAX/TMIN = max/min of temp-273.15;  RAIN = mean(prec)*24;
#   SRAD = mean(srad)*0.0864 MJ/m^2/day;  WIND = mean(wind);
#   RH/TDEW from daily-mean shum + temp + pres.
#
# ACCESS: free, but requires a (free) TPDC account (data.tpdc.ac.cn) to download
# the NetCDFs. Point cmfd_nc_dir at the folder. Coverage: China.

CMFD_TOKENS <- c("temp", "prec", "srad", "shum", "pres", "wind")

.cmfd_rh_tdew <- function(shum, temp_c, pres_pa) {
  q <- shum; p_hpa <- pres_pa / 100; tc <- temp_c
  e <- q * p_hpa / (0.622 + 0.378 * q)
  es <- 6.112 * exp(17.67 * tc / (tc + 243.5))
  rh <- pmin(pmax(100 * e / es, 0), 100)
  ln <- ifelse(e > 0, log(e / 6.112), NA_real_)
  tdew <- (243.5 * ln) / (17.67 - ln)
  list(rh = rh, tdew = tdew)
}

.cmfd_find_files <- function(nc_dir, token) {
  if (is.null(nc_dir) || !dir.exists(nc_dir)) return(character(0))
  files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  base <- tolower(basename(files))
  files[startsWith(base, paste0(token, "_")) | grepl(paste0("_", token, "_"), base) |
        grepl(paste0(token, "-"), base)]
}

# Per-point daily reductions for one variable: returns {pid: data.frame(date,...)}.
.cmfd_extract_daily <- function(paths, token, ids, pts_vect) {
  acc <- setNames(lapply(ids, function(.) NULL), ids)   # per-pid list of (time,val)
  for (path in paths) {
    r <- terra::rast(path)
    tt <- terra::time(r)
    if (all(is.na(tt))) next
    ex <- terra::extract(r, pts_vect, ID = FALSE)        # points x timesteps
    for (j in seq_along(ids)) {
      df <- data.frame(t = as.POSIXct(tt), v = as.numeric(ex[j, ]))
      acc[[ids[j]]] <- rbind(acc[[ids[j]]], df)
    }
  }
  out <- setNames(vector("list", length(ids)), ids)
  for (pid in ids) {
    d <- acc[[pid]]
    if (is.null(d) || nrow(d) == 0) { out[[pid]] <- data.frame(); next }
    d <- d[!duplicated(d$t), ]
    day <- as.Date(d$t)
    agg <- switch(token,
      temp = data.frame(date = sort(unique(day)),
                        TMAX = tapply(d$v, day, max, na.rm = TRUE) - 273.15,
                        TMIN = tapply(d$v, day, min, na.rm = TRUE) - 273.15,
                        TMEAN = tapply(d$v, day, mean, na.rm = TRUE) - 273.15),
      prec = data.frame(date = sort(unique(day)), RAIN = tapply(d$v, day, mean, na.rm = TRUE) * 24),
      srad = data.frame(date = sort(unique(day)), SRAD = tapply(d$v, day, mean, na.rm = TRUE) * 0.0864),
      wind = data.frame(date = sort(unique(day)), WIND = tapply(d$v, day, mean, na.rm = TRUE)),
      shum = data.frame(date = sort(unique(day)), SHUM = tapply(d$v, day, mean, na.rm = TRUE)),
      pres = data.frame(date = sort(unique(day)), PRES = tapply(d$v, day, mean, na.rm = TRUE)),
      data.frame(date = sort(unique(day)), V = tapply(d$v, day, mean, na.rm = TRUE)))
    rownames(agg) <- NULL
    out[[pid]] <- agg
  }
  out
}

.cmfd_write_wth <- function(df, pid, lat, lon, output_dir) {
  tavg <- (df$TMAX + df$TMIN) / 2
  tav <- mean(tavg, na.rm = TRUE)
  mon <- tapply(tavg, list(df$YEAR, df$MM), mean, na.rm = TRUE)
  amps <- apply(mon, 1, function(r) { r <- r[is.finite(r)]; if (length(r)) max(r) - min(r) else NA })
  amp <- mean(amps, na.rm = TRUE)
  header <- sprintf(
    "$WEATHER DATA: CMFD (Point ID: %s)\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  CMFD %8.4f %8.4f   -99 %5.1f %5.1f   2.0  10.0\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
    pid, lat, lon, tav, amp)
  clamp <- function(x) ifelse(!is.na(x) & (x >= 9999.95 | x <= -999.95), -99, x)
  d <- df; for (c in c("SRAD","TMAX","TMIN","RAIN","TDEW","RH2M","WIND")) d[[c]][is.na(d[[c]])] <- -99
  lines <- sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                   d$DATE, clamp(d$SRAD), clamp(d$TMAX), clamp(d$TMIN),
                   clamp(d$RAIN), clamp(d$TDEW), clamp(d$RH2M), clamp(d$WIND))
  lines <- gsub("-99.0", "  -99", lines, fixed = TRUE)
  writeLines(c(header, lines), file.path(output_dir, sprintf("%s.WTH", pid)))
}

process_weather_cmfd <- function(shapefile, start_year, end_year, output_dir,
                                 id_col, lat_col, lon_col, n_cores, log_file,
                                 cmfd_nc_dir = "") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  end_year <- min(end_year, as.integer(format(Sys.Date(), "%Y")))
  if (!nzchar(cmfd_nc_dir) || !dir.exists(cmfd_nc_dir))
    stop("CMFD needs cmfd_nc_dir pointing at a folder of CMFD NetCDF files ",
         "(temp/prec/srad/shum/pres/wind). Download from the TPDC (free account).")

  pts <- sf::st_transform(shapefile, 4326)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  xy <- sf::st_coordinates(pts); lats <- xy[, 2]; lons <- xy[, 1]
  pts_vect <- terra::vect(pts)

  message(sprintf("--- Starting CMFD Processing (Years: %d-%d) ---", start_year, end_year))
  var_files <- lapply(CMFD_TOKENS, function(tok) .cmfd_find_files(cmfd_nc_dir, tok))
  names(var_files) <- CMFD_TOKENS
  if (length(var_files$temp) == 0) stop("CMFD requires at least the temp NetCDF files.")

  daily <- list()
  for (tok in CMFD_TOKENS) if (length(var_files[[tok]]) > 0)
    daily[[tok]] <- .cmfd_extract_daily(var_files[[tok]], tok, ids, pts_vect)

  written <- 0
  for (k in seq_along(ids)) {
    pid <- ids[k]; lat <- lats[k]; lon <- lons[k]
    tryCatch({
      tdf <- daily[["temp"]][[pid]]
      if (is.null(tdf) || nrow(tdf) == 0) stop("No CMFD temperature data extracted (point outside China grid).")
      frame <- tdf
      jn <- function(tok) { d <- daily[[tok]][[pid]]; if (!is.null(d) && nrow(d)) merge(frame, d, by = "date", all.x = TRUE) else frame }
      for (tok in c("prec", "srad", "wind")) frame <- jn(tok)
      sh <- daily[["shum"]][[pid]]; pr <- daily[["pres"]][[pid]]
      if (!is.null(sh) && nrow(sh) && !is.null(pr) && nrow(pr)) {
        m <- merge(merge(frame[, c("date", "TMEAN")], sh, by = "date", all.x = TRUE), pr, by = "date", all.x = TRUE)
        rt <- .cmfd_rh_tdew(m$SHUM, m$TMEAN, m$PRES)
        frame <- merge(frame, data.frame(date = m$date, RH2M = rt$rh, TDEW = rt$tdew), by = "date", all.x = TRUE)
      }
      yr <- as.integer(format(frame$date, "%Y"))
      frame <- frame[yr >= start_year & yr <= end_year, ]
      if (nrow(frame) == 0) stop("No CMFD data in the requested year range.")
      frame$DATE <- sprintf("%d%03d", as.integer(format(frame$date, "%Y")), as.integer(format(frame$date, "%j")))
      frame$YEAR <- as.integer(format(frame$date, "%Y")); frame$MM <- as.integer(format(frame$date, "%m"))
      for (c in c("SRAD","RAIN","RH2M","WIND","TDEW")) if (is.null(frame[[c]])) frame[[c]] <- -99
      frame <- frame[!is.na(frame$TMAX) & !is.na(frame$TMIN), ]
      .cmfd_write_wth(frame, pid, lat, lon, output_dir)
      written <- written + 1
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\nCMFD point %s (%.3f,%.3f): %s\n", pid, lat, lon, conditionMessage(e))
      cat(msg); write(msg, file = log_file, append = TRUE)
    })
  }
  message(sprintf("\nCMFD processing complete: %d/%d points written to '%s'.\n", written, length(ids), output_dir))
}
