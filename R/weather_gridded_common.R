# Internal helpers for local/cache-backed gridded weather sources.

weather_calc_tav <- function(df) mean((df$TMAX + df$TMIN) / 2, na.rm = TRUE)

weather_calc_amp <- function(df) {
  tavg <- (df$TMAX + df$TMIN) / 2
  mon <- tapply(tavg, list(df$YEAR, df$MM), mean, na.rm = TRUE)
  amps <- apply(mon, 1, function(r) { r <- r[is.finite(r)]; if (length(r)) max(r) - min(r) else NA })
  mean(amps, na.rm = TRUE)
}

weather_tdew_from_rh <- function(tmean_c, rh_pct) {
  rh <- pmin(pmax(rh_pct, 1), 100)
  a <- 17.625; b <- 243.04
  gamma <- log(rh / 100) + (a * tmean_c) / (b + tmean_c)
  (b * gamma) / (a - gamma)
}

weather_convert_units <- function(vals, units, kind) {
  u <- tolower(ifelse(is.null(units), "", units))
  if (kind == "temp") {
    if (grepl("\\bk\\b|kelvin", u) || stats::median(vals, na.rm = TRUE) > 100) vals <- vals - 273.15
  } else if (kind == "rain") {
    if (grepl("s-1|/s", u)) vals <- vals * 86400
  } else if (kind == "srad") {
    if (grepl("w", u) && grepl("m", u)) vals <- vals * 0.0864
    else if (grepl("j", u) && grepl("m", u)) vals <- vals / 1000000
  } else if (kind == "wind") {
    # 10 m reanalysis wind -> 2 m (FAO-56 log profile factor ~0.748).
    vals <- vals * 0.748
  } else if (kind == "vp") {
    # Vapour pressure (hPa) -> dewpoint (degC), inverse Magnus formula.
    e <- pmax(vals, 1e-3)
    ln <- log(e / 6.1094)
    vals <- (243.04 * ln) / (17.625 - ln)
  }
  vals
}

weather_write_wth <- function(df, pid, lat, lon, output_dir, source_label,
                              insi, refht = 2.0, wndht = 2.0) {
  tav <- weather_calc_tav(df); amp <- weather_calc_amp(df)
  header <- sprintf(
    "$WEATHER DATA: %s (Point ID: %s)\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  %-4s %8.4f %8.4f   -99 %5.1f %5.1f %5.1f %5.1f\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
    source_label, pid, insi, lat, lon, tav, amp, refht, wndht)
  d <- df
  for (nm in c("SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND")) d[[nm]][is.na(d[[nm]])] <- -99
  lines <- sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                   d$DATE, d$SRAD, d$TMAX, d$TMIN, d$RAIN, d$TDEW, d$RH2M, d$WIND)
  lines <- gsub("-99.0", "  -99", lines, fixed = TRUE)
  writeLines(c(header, lines), file.path(output_dir, sprintf("%s.WTH", pid)))
}

weather_find_nc_files <- function(nc_dir, tokens) {
  if (!nzchar(nc_dir) || !dir.exists(nc_dir)) return(character())
  files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  stems <- tolower(tools::file_path_sans_ext(basename(files)))
  components <- strsplit(stems, "[^a-z0-9]+")
  hit <- vapply(components, function(parts) any(tolower(tokens) %in% parts), logical(1))
  files[hit]
}

weather_find_nc_file <- function(nc_dir, tokens) {
  files <- weather_find_nc_files(nc_dir, tokens)
  if (length(files)) files[1] else NA_character_
}

weather_extract_netcdf_series <- function(path, ids, pts_vect, start_year, end_year, kind) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("package 'terra' required for gridded NetCDF weather")
  r <- terra::rast(path)
  tt <- as.Date(terra::time(r))
  if (all(is.na(tt))) tt <- as.Date(terra::time(r), origin = "1970-01-01")
  yr <- as.integer(format(tt, "%Y"))
  keep <- which(yr >= start_year & yr <= end_year)
  if (!length(keep)) return(setNames(vector("list", length(ids)), ids))
  r <- r[[keep]]; tt <- tt[keep]
  units <- tryCatch(terra::units(r)[1], error = function(e) "")
  ex <- terra::extract(r, pts_vect, ID = FALSE)
  codes <- sprintf("%d%03d", as.integer(format(tt, "%Y")), as.integer(format(tt, "%j")))
  out <- setNames(vector("list", length(ids)), ids)
  for (i in seq_along(ids)) {
    vals <- weather_convert_units(as.numeric(ex[i, ]), units, kind)
    good <- is.finite(vals)
    v <- vals[good]; names(v) <- codes[good]
    out[[ids[i]]] <- v
  }
  out
}

process_local_netcdf_weather <- function(shapefile, start_year, end_year, output_dir,
                                         id_col, lat_col, lon_col, log_file,
                                         nc_dir, var_specs, source_label, insi,
                                         refht = 2.0, wndht = 2.0) {
  if (!nzchar(nc_dir) || !dir.exists(nc_dir)) stop(sprintf("%s needs local NetCDF directory: %s", source_label, nc_dir))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pts <- sf::st_transform(shapefile, 4326)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  xy <- sf::st_coordinates(pts)
  lats <- xy[, 2]; lons <- xy[, 1]
  pts_vect <- terra::vect(pts)
  per_var <- list()
  for (v in names(var_specs)) {
    spec <- var_specs[[v]]
    paths <- weather_find_nc_files(nc_dir, spec$tokens)
    if (!length(paths)) {
      if (isTRUE(spec$required)) stop(sprintf("%s required variable %s not found in %s", source_label, v, nc_dir))
      message(sprintf("  %s: no NetCDF for %s; writing -99 where needed.", source_label, v))
      next
    }
    per_var[[v]] <- weather_extract_netcdf_series(paths, ids, pts_vect, start_year, end_year, spec$kind)
  }
  missing_forcing <- setdiff(c("TMAX", "TMIN", "RAIN", "SRAD"), names(per_var))
  if (length(missing_forcing)) {
    stop(sprintf("%s requires %s; refusing to write WTH with missing forcing",
                 source_label, paste(missing_forcing, collapse = ", ")), call. = FALSE)
  }
  written <- 0
  for (k in seq_along(ids)) {
    pid <- ids[k]; lat <- lats[k]; lon <- lons[k]
    tryCatch({
      tmax <- per_var[["TMAX"]][[pid]]
      tmin <- per_var[["TMIN"]][[pid]]
      if (is.null(tmax) || !length(tmax) || is.null(tmin)) stop("No required TMAX/TMIN series extracted.")
      dates <- names(tmax)
      expected_end <- as.Date(sprintf("%d-12-31", end_year))
      if (end_year == as.integer(format(Sys.Date(), "%Y"))) expected_end <- Sys.Date()
      expected <- format(seq(as.Date(sprintf("%d-01-01", start_year)), expected_end, by = "day"), "%Y%j")
      for (required in c("TMAX", "TMIN", "RAIN", "SRAD")) {
        actual <- names(per_var[[required]][[pid]])
        missing <- setdiff(expected, actual)
        if (length(missing)) stop(sprintf("%s incomplete (%d missing day(s))", required, length(missing)))
      }
      grab <- function(v) { s <- per_var[[v]][[pid]]; if (is.null(s)) rep(NA_real_, length(dates)) else as.numeric(s[dates]) }
      df <- data.frame(DATE = dates, YEAR = as.integer(substr(dates, 1, 4)),
                       MM = as.integer(format(as.Date(dates, "%Y%j"), "%m")),
                       TMAX = as.numeric(tmax), TMIN = grab("TMIN"),
                       RAIN = grab("RAIN"), SRAD = grab("SRAD"),
                       TDEW = grab("TDEW"), RH2M = grab("RH2M"), WIND = grab("WIND"))
      if (all(is.na(df$TDEW)) && !is.null(per_var[["TMEAN"]]) && !is.null(per_var[["RH2M"]]))
        df$TDEW <- weather_tdew_from_rh(grab("TMEAN"), df$RH2M)
      df <- df[!is.na(df$TMAX) & !is.na(df$TMIN), ]
      weather_write_wth(df, pid, lat, lon, output_dir, source_label, insi, refht, wndht)
      written <- written + 1
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\n%s point %s (%.3f,%.3f): %s\n", source_label, pid, lat, lon, conditionMessage(e))
      cat(msg); write(msg, file = log_file, append = TRUE)
    })
  }
  written
}
