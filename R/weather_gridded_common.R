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

weather_convert_units <- function(vals, units, kind, wind_height_m = NULL) {
  u <- gsub("per", "/", gsub("[[:space:]_^*]", "", tolower(if (is.null(units)) "" else units)), fixed = TRUE)
  unsupported <- function() stop(sprintf("Unsupported or missing %s units: %s", kind, units))
  if (kind == "temp") {
    if (u %in% c("k", "kelvin", "degk", "degreek", "degreeskelvin")) return(vals - 273.15)
    if (u %in% c("c", "degc", "celsius", "degreecelsius", "degreescelsius", "°c")) return(vals)
    unsupported()
  }
  if (kind == "rain") {
    if (u %in% c("mm", "mm/day", "mmday-1", "mmd-1", "kgm-2", "kg/m2")) return(vals)
    if (u %in% c("kgm-2s-1", "kg/m2/s", "mms-1", "mm/s")) return(vals * 86400)
    if (u %in% c("m", "m/day", "mday-1")) return(vals * 1000)
    unsupported()
  }
  if (kind == "srad") {
    if (u %in% c("wm-2", "w/m2")) return(vals * 0.0864)
    for (prefix in c("mj", "kj", "j")) {
      if (u %in% paste0(prefix, c("m-2", "/m2", "m-2day-1", "m-2d-1", "/m2/day", "m-2/day"))) {
        return(vals * c(mj = 1, kj = .001, j = 1e-6)[[prefix]])
      }
    }
    unsupported()
  }
  if (kind == "wind") {
    if (u %in% c("km/h", "kmh-1")) vals <- vals / 3.6
    else if (!u %in% c("m/s", "ms-1")) unsupported()
    if (is.null(wind_height_m) || !is.finite(wind_height_m) || wind_height_m < .1) stop("Wind requires an explicit height_m in metres")
    if (wind_height_m == 2) return(vals)
    return(vals * log(67.8 * 2 - 5.42) / log(67.8 * wind_height_m - 5.42))
  }
  if (kind == "vp") {
    if (u == "pa") vals <- vals / 100
    else if (u == "kpa") vals <- vals * 10
    else if (!u %in% c("hpa", "mbar")) unsupported()
    vals[vals <= 0] <- NA_real_
    ln <- log(vals / 6.1094)
    return(243.04 * ln / (17.625 - ln))
  }
  if (kind %in% c("rh", "rh2m")) {
    if (u %in% c("1", "fraction")) return(vals * 100)
    if (u %in% c("%", "percent", "percentage")) return(vals)
    unsupported()
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
  lines <- sprintf("%7s%s%s%s%s%s%s%s",
                   d$DATE, .format_wth_value(d$SRAD), .format_wth_value(d$TMAX), .format_wth_value(d$TMIN), .format_wth_value(d$RAIN), .format_wth_value(d$TDEW), .format_wth_value(d$RH2M), .format_wth_value(.weather_wind_run(d$WIND)))
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

weather_extract_netcdf_series <- function(path, ids, pts_vect, start_year, end_year, kind, wind_height_m = NULL) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("package 'terra' required for gridded NetCDF weather")
  r <- terra::rast(path)
  if (!terra::is.lonlat(r) || !terra::is.lonlat(pts_vect)) stop("Weather coordinates must be geographic longitude/latitude")
  timestamps <- terra::time(r)
  tt <- as.Date(timestamps)
  if (all(is.na(tt))) tt <- as.Date(terra::time(r), origin = "1970-01-01")
  yr <- as.integer(format(tt, "%Y"))
  keep <- which(yr >= start_year & yr <= end_year)
  if (!length(keep)) return(setNames(vector("list", length(ids)), ids))
  if (anyDuplicated(timestamps[keep])) {
    stop("Duplicate weather timestamps; remove overlapping records before extraction")
  }
  if (anyDuplicated(tt[keep])) {
    stop("Subdaily weather data is unsupported; provide one daily record per date")
  }
  r <- r[[keep]]; tt <- tt[keep]
  units <- tryCatch(terra::units(r)[1], error = function(e) "")
  if (kind == "wind" && is.null(wind_height_m)) {
    nc <- ncdf4::nc_open(path[1]); on.exit(ncdf4::nc_close(nc), add = TRUE)
    vars <- names(nc$var)
    heights <- lapply(vars, function(v) ncdf4::ncatt_get(nc, v, "height_m"))
    valid <- vapply(heights, function(x) isTRUE(x$hasatt), logical(1))
    if (sum(valid) == 1) wind_height_m <- as.numeric(heights[[which(valid)]]$value)
  }
  e <- terra::ext(r)
  res <- terra::res(r)
  crds <- terra::crds(pts_vect)
  lat_res <- 0
  lon_res <- 0
  if (e$xmin >= 0) {
    crds[, 1] <- ifelse(crds[, 1] < 0, crds[, 1] + 360, crds[, 1])
  } else if (e$xmax <= 180) {
    crds[, 1] <- ifelse(crds[, 1] > 180, crds[, 1] - 360, crds[, 1])
  }
  bounds <- c(e$xmin + res[1]/2, e$xmax - res[1]/2, e$ymin + res[2]/2, e$ymax - res[2]/2)
  for (i in seq_along(ids)) {
    lon_val <- crds[i, 1]; lat_val <- crds[i, 2]
    if (lat_val < (bounds[3] - lat_res) || lat_val > (bounds[4] + lat_res) ||
        lon_val < (bounds[1] - lon_res) || lon_val > (bounds[2] + lon_res)) {
      stop(sprintf("Point %s (%.4f, %.4f) is outside grid domain: lat [%.4f, %.4f], lon [%.4f, %.4f]",
                   ids[i], lat_val, lon_val, e$ymin, e$ymax, e$xmin, e$xmax))
    }
  }
  pts_vect <- terra::vect(crds, crs = terra::crs(pts_vect))
  ex <- terra::extract(r, pts_vect, ID = FALSE)
  codes <- sprintf("%d%03d", as.integer(format(tt, "%Y")), as.integer(format(tt, "%j")))
  out <- setNames(vector("list", length(ids)), ids)
  for (i in seq_along(ids)) {
    vals <- weather_convert_units(as.numeric(ex[i, ]), units, kind, wind_height_m)
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
    per_var[[v]] <- weather_extract_netcdf_series(paths, ids, pts_vect, start_year, end_year, spec$kind, wind_height_m = spec$height_m)
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
