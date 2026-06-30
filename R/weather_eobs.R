# File: weather_eobs.R   (R twin of python/dssatutils/weather_eobs.py)
#
# Weather source: E-OBS (ECA&D European gridded daily observations) -> DSSAT .WTH.
#
# WHY: E-OBS is the standard high-resolution (0.1°, ~11 km) daily gauge-based
# gridded dataset for Europe, 1950-present, and (unlike most gridded products)
# includes daily GLOBAL RADIATION (qq) so DSSAT's SRAD comes straight from data.
#
# E-OBS variables (one NetCDF per variable): tx/tn/tg (°C), rr (mm),
# qq (global radiation, W/m²), fg (wind m/s), hu (relative humidity %).
#
# ACCESS — two modes:
#   (A) LOCAL (default, no key): point eobs_nc_dir at pre-downloaded E-OBS NetCDF
#       files (www.ecad.eu/download/ensembles/download.php).
#   (B) CDS (optional, eobs_use_cds=TRUE): area/time subset via the Copernicus
#       CDS dataset "insitu-gridded-observations-europe" (needs ~/.cdsapirc, like
#       AgERA5). Uses ecmwfr and the shared setup_cds_credentials() helper.
#
# The NetCDF extraction (.eobs_extract_points) and .WTH writer (.eobs_write_wth)
# are isolated from any network/credential dependency and unit-testable.

# DSSAT var -> filename token used to locate the E-OBS NetCDF for that variable.
EOBS_VARS <- c(TMAX = "tx", TMIN = "tn", TMEAN = "tg", RAIN = "rr",
               SRAD = "qq", WIND = "fg", RH2M = "hu")


.eobs_tdew_from_rh <- function(tmean_c, rh_pct) {
  rh <- pmin(pmax(rh_pct, 1), 100)
  a <- 17.625; b <- 243.04
  gamma <- log(rh / 100) + (a * tmean_c) / (b + tmean_c)
  (b * gamma) / (a - gamma)
}

.eobs_find_var_file <- function(eobs_nc_dir, token) {
  if (is.null(eobs_nc_dir) || !dir.exists(eobs_nc_dir)) return(NA_character_)
  files <- list.files(eobs_nc_dir, pattern = "\\.nc$", full.names = TRUE)
  base <- tolower(basename(files))
  hit <- files[startsWith(base, paste0(token, "_")) |
               grepl(paste0("_", token, "_"), base) |
               grepl(paste0(token, "_ens"), base)]
  if (length(hit)) hit[1] else NA_character_
}

.eobs_download_cds <- function(token, year_start, year_end, area, cache_dir,
                               cds_user = "ecmwfr") {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(cache_dir, sprintf("eobs_%s_%d_%d.nc", token, year_start, year_end))
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)

  .dssatutils_ensure_cds_credentials(user = cds_user, require_ecmwfr = TRUE)
  name_map <- c(
    tx = "maximum_temperature",
    tn = "minimum_temperature",
    tg = "mean_temperature",
    rr = "precipitation_amount",
    qq = "surface_shortwave_downwelling_radiation",
    fg = "wind_speed",
    hu = "relative_humidity"
  )
  req <- list(
    dataset_short_name = "insitu-gridded-observations-europe",
    product_type = "ensemble_mean",
    variable = unname(if (token %in% names(name_map)) name_map[[token]] else token),
    grid_resolution = "0.1deg",
    period = "full_period",
    version = "30.0e",
    format = "netcdf",
    area = as.numeric(area),
    target = basename(dest)
  )

  err <- NULL
  tryCatch(
    ecmwfr::wf_request(request = req, user = cds_user, transfer = TRUE,
                       path = cache_dir, verbose = FALSE),
    error = function(e) err <<- conditionMessage(e)
  )
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)
  message(sprintf("  E-OBS CDS download failed (%s): %s",
                  token, if (is.null(err)) "no data file returned" else err))
  NA_character_
}


# --- Per-point daily series ({pid -> named vector by YYYYDOY}) for one var ----
.eobs_extract_points <- function(path, dssat_var, ids, pts_vect, year_start, year_end) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("package 'terra' required for E-OBS")
  r <- terra::rast(path)
  tt <- terra::time(r)
  if (all(is.na(tt))) tt <- as.Date(terra::time(r), origin = "1970-01-01")
  tt <- as.Date(tt)
  yr <- as.integer(format(tt, "%Y"))
  keep <- which(yr >= year_start & yr <= year_end)
  if (length(keep) == 0) return(setNames(vector("list", length(ids)), ids))
  r <- r[[keep]]; tt <- tt[keep]
  ex <- terra::extract(r, pts_vect, ID = FALSE)          # rows = points, cols = times
  date_codes <- sprintf("%d%03d", as.integer(format(tt, "%Y")), as.integer(format(tt, "%j")))
  out <- setNames(vector("list", length(ids)), ids)
  for (j in seq_along(ids)) {
    vals <- as.numeric(ex[j, ])
    if (dssat_var == "SRAD") vals <- vals * 0.0864        # W/m² -> MJ/m²/day
    good <- is.finite(vals)
    v <- vals[good]; names(v) <- date_codes[good]
    out[[ids[j]]] <- v
  }
  out
}


# --- .WTH writer (shared format with the other weather modules) ---------------
.eobs_write_wth <- function(df, pid, lat, lon, output_dir) {
  tavg <- (df$TMAX + df$TMIN) / 2
  tav <- mean(tavg, na.rm = TRUE)
  mon <- tapply(tavg, list(df$YEAR, df$MM), mean, na.rm = TRUE)
  amps <- apply(mon, 1, function(r) { r <- r[is.finite(r)]; if (length(r)) max(r) - min(r) else NA })
  amp <- mean(amps, na.rm = TRUE)
  header <- sprintf(
    "$WEATHER DATA: E-OBS (Point ID: %s)\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n  EOBS %8.4f %8.4f   -99 %5.1f %5.1f   2.0  10.0\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
    pid, lat, lon, tav, amp)
  clamp <- function(x) ifelse(!is.na(x) & (x >= 9999.95 | x <= -999.95), -99, x)
  d <- df; for (c in c("SRAD","TMAX","TMIN","RAIN","TDEW","RH2M","WIND")) d[[c]][is.na(d[[c]])] <- -99
  lines <- sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                   d$DATE, clamp(d$SRAD), clamp(d$TMAX), clamp(d$TMIN),
                   clamp(d$RAIN), clamp(d$TDEW), clamp(d$RH2M), clamp(d$WIND))
  lines <- gsub("-99.0", "  -99", lines, fixed = TRUE)
  writeLines(c(header, lines), file.path(output_dir, sprintf("%s.WTH", pid)))
}


# --- Public entry point -------------------------------------------------------
process_weather_eobs <- function(shapefile, start_year, end_year, output_dir,
                                 id_col, lat_col, lon_col, n_cores, log_file,
                                 eobs_nc_dir = "", eobs_cache_dir = "",
                                 eobs_use_cds = FALSE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  end_year <- min(end_year, as.integer(format(Sys.Date(), "%Y")))

  coords <- .extract_coords(shapefile, id_col, lat_col, lon_col)
  ids <- coords$ids
  lats <- coords$lats
  lons <- coords$lons
  pts_vect <- terra::vect(data.frame(lon = lons, lat = lats),
                          geom = c("lon", "lat"), crs = "EPSG:4326")

  message(sprintf("--- Starting E-OBS Processing (Years: %d-%d) ---", start_year, end_year))

  if (isTRUE(eobs_use_cds)) {
    if (!nzchar(eobs_cache_dir)) eobs_cache_dir <- file.path(output_dir, "eobs_cache")
    pad <- 0.2
    area <- c(max(lats) + pad, min(lons) - pad, min(lats) - pad, max(lons) + pad)
    message("  E-OBS via Copernicus CDS (requires a configured CDS token).")
    var_paths <- vapply(EOBS_VARS, function(tok) {
      .eobs_download_cds(tok, start_year, end_year, area, eobs_cache_dir)
    }, character(1))
    names(var_paths) <- names(EOBS_VARS)
  } else {
    if (!nzchar(eobs_nc_dir) || !dir.exists(eobs_nc_dir))
      stop("E-OBS local mode needs eobs_nc_dir pointing at a folder of E-OBS NetCDF ",
           "files (tx/tn/rr/qq...). Download from www.ecad.eu.")
    var_paths <- vapply(EOBS_VARS, function(tok) .eobs_find_var_file(eobs_nc_dir, tok), character(1))
    names(var_paths) <- names(EOBS_VARS)
  }
  if (is.na(var_paths[["TMAX"]]) || is.na(var_paths[["TMIN"]]))
    stop("E-OBS requires at least tx (TMAX) and tn (TMIN) NetCDF files.")

  per_var <- list()
  for (v in names(EOBS_VARS)) {
    p <- var_paths[[v]]
    if (is.na(p)) {
      message(sprintf("  E-OBS: no file for %s (%s); it will be written as -99.", v, EOBS_VARS[[v]]))
      next
    }
    per_var[[v]] <- .eobs_extract_points(p, v, ids, pts_vect, start_year, end_year)
  }

  written <- 0
  for (k in seq_along(ids)) {
    pid <- ids[k]; lat <- lats[k]; lon <- lons[k]
    tryCatch({
      tmax <- per_var[["TMAX"]][[pid]]; tmin <- per_var[["TMIN"]][[pid]]
      if (is.null(tmax) || length(tmax) == 0)
        stop("No E-OBS data extracted for this point (outside grid / land mask).")
      dates <- names(tmax)
      grab <- function(v) { s <- per_var[[v]][[pid]]; if (is.null(s)) rep(NA_real_, length(dates)) else as.numeric(s[dates]) }
      df <- data.frame(
        DATE = dates,
        YEAR = as.integer(substr(dates, 1, 4)),
        MM   = as.integer(format(as.Date(dates, "%Y%j"), "%m")),
        TMAX = as.numeric(tmax),
        TMIN = grab("TMIN"), RAIN = grab("RAIN"),
        SRAD = grab("SRAD"), WIND = grab("WIND"), RH2M = grab("RH2M"),
        stringsAsFactors = FALSE)
      tmean <- grab("TMEAN")
      df$TDEW <- if (all(is.na(tmean)) || all(is.na(df$RH2M))) -99 else .eobs_tdew_from_rh(tmean, df$RH2M)
      df <- df[!is.na(df$TMAX) & !is.na(df$TMIN), ]
      .eobs_write_wth(df, pid, lat, lon, output_dir)
      written <- written + 1
    }, error = function(e) {
      msg <- sprintf("\n--- ERROR ---\nE-OBS point %s (%.3f,%.3f): %s\n", pid, lat, lon, conditionMessage(e))
      cat(msg); write(msg, file = log_file, append = TRUE)
    })
  }
  message(sprintf("\nE-OBS processing complete: %d/%d points written to '%s'.\n", written, length(ids), output_dir))
}
