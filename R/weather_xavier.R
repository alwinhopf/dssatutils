# File: weather_xavier.R   (R twin of python/dssatutils/weather_xavier.py)
#
# Weather source: BR-DWGD (Xavier et al.) Brazilian Daily Weather Gridded Data
# -> DSSAT .WTH.
#
# WHY: BR-DWGD is the standard 0.1° daily gauge-interpolated weather product for
# BRAZIL (1961-present). It already provides daily solar radiation in MJ/m^2/day
# (no estimation), plus Tmax/Tmin, precipitation, RH and 2 m wind.
#
# Variables (one NetCDF per variable): Tmax/Tmin (°C), pr (mm), Rs (MJ/m^2/day),
# RH (%), u2 (m/s @2m).
#
# ACCESS (open, NO key): download the NetCDFs from the BR-DWGD site and point
# xavier_nc_dir at the folder. Coverage: Brazil.

XAVIER_VARS <- c(TMAX = "Tmax", TMIN = "Tmin", RAIN = "pr",
                 SRAD = "Rs", RH2M = "RH", WIND = "u2")

.xavier_tdew_from_rh <- function(tmean_c, rh_pct) {
  rh <- pmin(pmax(rh_pct, 1), 100); a <- 17.625; b <- 243.04
  gamma <- log(rh / 100) + (a * tmean_c) / (b + tmean_c)
  (b * gamma) / (a - gamma)
}

.xavier_find_files <- function(nc_dir, token) {
  if (is.null(nc_dir) || !dir.exists(nc_dir)) return(character(0))
  files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  base <- basename(files)
  files[startsWith(base, paste0(token, "_")) | grepl(paste0("_", token, "_"), base) |
        base == paste0(token, ".nc")]
}

# Per-point daily series across one variable's (possibly several) files.
.xavier_extract <- function(paths, dssat_var, ids, pts_vect, year_start, year_end) {
  out <- setNames(lapply(ids, function(.) numeric(0)), ids)
  for (path in paths) {
    r <- terra::rast(path)
    tt <- terra::time(r)
    if (all(is.na(tt))) next
    tt <- as.Date(tt)
    yr <- as.integer(format(tt, "%Y"))
    keep <- which(yr >= year_start & yr <= year_end)
    if (length(keep) == 0) next
    r <- r[[keep]]; tt <- tt[keep]
    ex <- terra::extract(r, pts_vect, ID = FALSE)
    codes <- sprintf("%d%03d", as.integer(format(tt, "%Y")), as.integer(format(tt, "%j")))
    for (j in seq_along(ids)) {
      vals <- as.numeric(ex[j, ])
      good <- is.finite(vals)
      v <- vals[good]; names(v) <- codes[good]
      out[[ids[j]]] <- c(out[[ids[j]]], v)
    }
  }
  out
}

.xavier_write_wth <- function(df, pid, lat, lon, output_dir) {
  tavg <- (df$TMAX + df$TMIN) / 2
  tav <- mean(tavg, na.rm = TRUE)
  mon <- tapply(tavg, list(df$YEAR, df$MM), mean, na.rm = TRUE)
  amps <- apply(mon, 1, function(r) { r <- r[is.finite(r)]; if (length(r)) max(r) - min(r) else NA })
  amp <- mean(amps, na.rm = TRUE)
  header <- sprintf(
    "$WEATHER DATA: BR-DWGD/Xavier (Point ID: %s)\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  XAVR %8.4f %8.4f   -99 %5.1f %5.1f   2.0   2.0\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
    pid, lat, lon, tav, amp)
  clamp <- function(x) ifelse(!is.na(x) & (x >= 9999.95 | x <= -999.95), -99, x)
  d <- df; for (c in c("SRAD","TMAX","TMIN","RAIN","TDEW","RH2M","WIND")) d[[c]][is.na(d[[c]])] <- -99
  lines <- sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                   d$DATE, clamp(d$SRAD), clamp(d$TMAX), clamp(d$TMIN),
                   clamp(d$RAIN), clamp(d$TDEW), clamp(d$RH2M), clamp(d$WIND))
  lines <- gsub("-99.0", "  -99", lines, fixed = TRUE)
  writeLines(c(header, lines), file.path(output_dir, sprintf("%s.WTH", pid)))
}

process_weather_xavier <- function(shapefile, start_year, end_year, output_dir,
                                   id_col, lat_col, lon_col, n_cores, log_file,
                                   xavier_nc_dir = "") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  end_year <- min(end_year, as.integer(format(Sys.Date(), "%Y")))
  if (!nzchar(xavier_nc_dir) || !dir.exists(xavier_nc_dir))
    stop("Xavier needs xavier_nc_dir pointing at a folder of BR-DWGD NetCDF files ",
         "(Tmax/Tmin/pr/Rs/RH/u2). Download from the BR-DWGD site (no key).")

  pts <- sf::st_transform(shapefile, 4326)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  xy <- sf::st_coordinates(pts); lats <- xy[, 2]; lons <- xy[, 1]
  pts_vect <- terra::vect(pts)

  message(sprintf("--- Starting BR-DWGD/Xavier Processing (Years: %d-%d) ---", start_year, end_year))
  var_files <- lapply(XAVIER_VARS, function(tok) .xavier_find_files(xavier_nc_dir, tok))
  names(var_files) <- names(XAVIER_VARS)
  if (length(var_files$TMAX) == 0 || length(var_files$TMIN) == 0)
    stop("Xavier requires at least Tmax and Tmin NetCDF files.")

  per_var <- list()
  for (v in names(XAVIER_VARS)) {
    if (length(var_files[[v]]) == 0) {
      message(sprintf("  Xavier: no file for %s (%s); it will be written as -99.", v, XAVIER_VARS[[v]]))
      next
    }
    per_var[[v]] <- .xavier_extract(var_files[[v]], v, ids, pts_vect, start_year, end_year)
  }

  written <- 0
  for (k in seq_along(ids)) {
    pid <- ids[k]; lat <- lats[k]; lon <- lons[k]
    tryCatch({
      tmax <- per_var[["TMAX"]][[pid]]; tmin <- per_var[["TMIN"]][[pid]]
      if (is.null(tmax) || length(tmax) == 0) stop("No Xavier data extracted (point outside Brazil grid).")
      dates <- names(tmax)
      grab <- function(v) { s <- per_var[[v]][[pid]]; if (is.null(s)) rep(NA_real_, length(dates)) else as.numeric(s[dates]) }
      df <- data.frame(DATE = dates, YEAR = as.integer(substr(dates, 1, 4)),
                       MM = as.integer(format(as.Date(dates, "%Y%j"), "%m")),
                       TMAX = as.numeric(tmax), TMIN = grab("TMIN"), RAIN = grab("RAIN"),
                       SRAD = grab("SRAD"), WIND = grab("WIND"), RH2M = grab("RH2M"),
                       stringsAsFactors = FALSE)
      tmean <- (df$TMAX + df$TMIN) / 2
      df$TDEW <- if (all(is.na(df$RH2M))) -99 else .xavier_tdew_from_rh(tmean, df$RH2M)
      df <- df[!is.na(df$TMAX) & !is.na(df$TMIN), ]
      .xavier_write_wth(df, pid, lat, lon, output_dir)
      written <- written + 1
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\nXavier point %s (%.3f,%.3f): %s\n", pid, lat, lon, conditionMessage(e))
      cat(msg); write(msg, file = log_file, append = TRUE)
    })
  }
  message(sprintf("\nXavier processing complete: %d/%d points written to '%s'.\n", written, length(ids), output_dir))
}
