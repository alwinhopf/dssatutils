# ==============================================================================
#  SOIL HELPER: LUCAS Topsoil (Europe, via ESDAC) — DSSAT .SOL
#  Filename: soil_lucas.R   (R twin of python/dssatutils/soil_lucas.py)
# ==============================================================================
#
#  LUCAS is the harmonised, ground-MEASURED topsoil survey for the EU (texture,
#  organic carbon, pH, ... at tens of thousands of georeferenced points).
#
#  LIMITATION: LUCAS samples the TOPSOIL only (0-20 cm). DSSAT needs a profile,
#  so the measured 0-20 cm layer is carried down to LUCAS_ROOTING_MAX_CM as an
#  explicit extrapolation (flagged in the .SOL header). Bulk density is not in
#  LUCAS; it is estimated from Saxton & Rawls porosity: BD = (1 - SSAT) * 2.65.
#
#  ACCESS: free, but ESDAC distributes LUCAS behind a one-off request form
#  (esdac.jrc.ec.europa.eu -> LUCAS Topsoil). Download the point table and point
#  lucas_csv at it. Coverage: EU.
# ==============================================================================

LUCAS_ROOTING_MAX_CM <- 150

.lucas_aliases <- list(
  id   = c("POINTID", "POINT_ID", "Point_ID", "id", "ID"),
  lat  = c("TH_LAT", "GPS_LAT", "lat", "Latitude", "latitude", "Y", "POINT_Y"),
  lon  = c("TH_LONG", "GPS_LONG", "lon", "Longitude", "longitude", "X", "POINT_X"),
  clay = c("clay", "Clay", "Clay_content", "clay_content"),
  sand = c("sand", "Sand", "Sand_content", "sand_content"),
  silt = c("silt", "Silt", "Silt_content", "silt_content"),
  oc   = c("OC", "oc", "OC_gkg", "organic_carbon"),
  ph   = c("pH_H2O", "pH_in_H2O", "pH_CaCl2", "pH", "ph")
)

.lucas_saxton_rawls <- function(sand_pct, clay_pct, om_pct) {
  S <- sand_pct / 100; C <- clay_pct / 100; OM <- om_pct / 100
  t1500 <- -0.024*S + 0.487*C + 0.006*OM + 0.005*(S*OM) - 0.013*(C*OM) + 0.068*(S*C) + 0.031
  # Floor wilting point at 0.02 (see soil_ssurgo.R): Saxton-Rawls can yield
  # SLLL <= 0 on very sandy soils, which SIGFPEs DSSAT's water balance.
  SLLL  <- pmax(t1500 + (0.14*t1500 - 0.02), 0.02)
  t33   <- -0.251*S + 0.195*C + 0.011*OM + 0.006*(S*OM) - 0.027*(C*OM) + 0.452*(S*C) + 0.299
  # Keep DUL-LL >= 0.04 (see soil_ssurgo.R): a near-zero gap SIGFPEs DSSAT.
  SDUL  <- pmax(t33 + (1.283*t33^2 - 0.374*t33 - 0.015), SLLL + 0.04)
  ts33t <- 0.278*S + 0.034*C + 0.022*OM - 0.018*(S*OM) - 0.027*(C*OM) - 0.584*(S*C) + 0.078
  ts33  <- ts33t + (0.636*ts33t - 0.107)
  SSAT  <- SDUL + ts33 - 0.097*S + 0.043
  c(SLLL = SLLL, SDUL = SDUL, SSAT = SSAT)
}

.lucas_resolve <- function(df, key, col_map) {
  if (!is.null(col_map[[key]]) && col_map[[key]] %in% names(df)) return(col_map[[key]])
  hit <- intersect(.lucas_aliases[[key]], names(df))
  if (length(hit)) hit[1] else NA_character_
}

.lucas_load <- function(lucas_csv, col_map) {
  raw <- if (grepl("\\.(xlsx|xls)$", lucas_csv, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("reading .xlsx LUCAS needs the 'readxl' package")
    as.data.frame(readxl::read_excel(lucas_csv))
  } else {
    utils::read.csv(lucas_csv, stringsAsFactors = FALSE, check.names = TRUE)
  }
  cols <- setNames(lapply(names(.lucas_aliases), function(k) .lucas_resolve(raw, k, col_map)), names(.lucas_aliases))
  for (need in c("lat", "lon", "clay", "sand")) {
    if (is.na(cols[[need]]))
      stop(sprintf("LUCAS table missing a '%s' column (looked for %s). Pass col_map.",
                   need, paste(.lucas_aliases[[need]], collapse = ", ")))
  }
  num <- function(col) suppressWarnings(as.numeric(raw[[col]]))
  out <- data.frame(
    src_id = if (!is.na(cols$id)) as.character(raw[[cols$id]]) else paste0("L", seq_len(nrow(raw))),
    lat = num(cols$lat), lon = num(cols$lon), clay = num(cols$clay), sand = num(cols$sand),
    stringsAsFactors = FALSE)
  out$silt <- if (!is.na(cols$silt)) num(cols$silt) else (100 - out$clay - out$sand)
  out$oc   <- if (!is.na(cols$oc)) num(cols$oc) else NA_real_
  out$ph   <- if (!is.na(cols$ph)) num(cols$ph) else NA_real_
  out <- out[is.finite(out$lat) & is.finite(out$lon) & is.finite(out$clay) & is.finite(out$sand), ]
  if (nrow(out) == 0) stop("No usable LUCAS rows after parsing (check column mapping / units).")
  out
}

format_dssat_soil_lucas <- function(profile_data, output_dir) {
  soil_id <- as.character(profile_data$ID[1])
  filename <- file.path(output_dir, paste0(soil_id, ".SOL"))
  if (file.exists(filename)) return()
  cat("*SOILS: Europe LUCAS Topsoil Profiles\n", file = filename)
  cat("! Generated from LUCAS topsoil (0-20 cm MEASURED; subsoil EXTRAPOLATED), Saxton & Rawls\n\n",
      file = filename, append = TRUE)
  cat(sprintf("*%-6s  LUCAS         %9.3f %9.3f\n",
              soil_id, profile_data$latitude[1], profile_data$longitude[1]),
      file = filename, append = TRUE)
  cat("@SITE        COUNTRY          LAT     LONG SCS FAMILY\n", file = filename, append = TRUE)
  cat(sprintf(" %-11s EU          %9.3f %9.3f \n",
              soil_id, profile_data$latitude[1], profile_data$longitude[1]),
      file = filename, append = TRUE)
  cat("@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE\n", file = filename, append = TRUE)
  cat("    BN   .13     6    .6    73     1     1 IB001 IB001 IB001\n", file = filename, append = TRUE)
  cat("@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC\n",
      file = filename, append = TRUE)
  prof <- profile_data[order(profile_data$depth_bottom), ]
  for (k in seq_len(nrow(prof))) {
    layer <- prof[k, ]
    slll <- sub("^0", " ", sprintf("%5.3f", layer$SLLL))
    sdul <- sub("^0", " ", sprintf("%5.3f", layer$SDUL))
    ssat <- sub("^0", " ", sprintf("%5.3f", layer$SSAT))
    cat(sprintf("%6d   -99 %s %s %s  1.00   -99 %5.2f %5.2f %5.1f %5.1f   -99   -99   -99   -99   -99   -99\n",
                as.integer(layer$depth_bottom), slll, sdul, ssat,
                layer$bulk_density, layer$om_pct / 1.724, layer$clay_pct, layer$silt_pct),
        file = filename, append = TRUE)
  }
  cat("\n", file = filename, append = TRUE)
}

process_soils_lucas <- function(grid_points, output_dir_csv, output_dir_individual, n_cores,
                                id_col, lat_col, long_col, format_sql_func = NULL,
                                lucas_csv = "", max_dist_km = 50, col_map = list()) {
  message("Starting LUCAS Topsoil Processing (Smart Resume Mode)...")
  if (!nzchar(lucas_csv) || !file.exists(lucas_csv))
    stop("LUCAS needs lucas_csv pointing at the downloaded ESDAC LUCAS topsoil table ",
         "(CSV/XLSX). Request it free at esdac.jrc.ec.europa.eu.")
  dir.create(output_dir_individual, recursive = TRUE, showWarnings = FALSE)
  lucas <- .lucas_load(lucas_csv, col_map)
  message(sprintf("LUCAS: loaded %d topsoil samples from %s", nrow(lucas), basename(lucas_csv)))

  pts <- sf::st_transform(grid_points, 4326)
  grid_df <- sf::st_drop_geometry(pts)
  xy <- sf::st_coordinates(pts)
  grid_df[[lat_col]] <- xy[, 2]; grid_df[[long_col]] <- xy[, 1]

  existing <- tools::file_path_sans_ext(list.files(output_dir_individual, pattern = "\\.SOL$"))
  todo <- grid_df[!(as.character(grid_df[[id_col]]) %in% existing), , drop = FALSE]
  message(sprintf("Resume Check: Found %d existing profiles. Processing %d remaining points.",
                  nrow(grid_df) - nrow(todo), nrow(todo)))
  if (nrow(todo) == 0) { message("All soil profiles already exist. Skipping LUCAS processing."); return(TRUE) }

  results <- list(); fails <- list()
  for (i in seq_len(nrow(todo))) {
    ID <- as.character(todo[[id_col]][i]); lat <- as.numeric(todo[[lat_col]][i]); lon <- as.numeric(todo[[long_col]][i])
    if (file.exists(file.path(output_dir_individual, paste0(ID, ".SOL")))) next
    dlat <- (lucas$lat - lat) * pi / 180; dlon <- (lucas$lon - lon) * pi / 180
    a <- sin(dlat / 2)^2 + cos(lat * pi / 180) * cos(lucas$lat * pi / 180) * sin(dlon / 2)^2
    dist <- 6371 * 2 * asin(pmin(1, sqrt(a)))
    k <- which.min(dist)
    if (dist[k] > max_dist_km) {
      fails[[length(fails) + 1]] <- list(ID = ID, latitude = lat, longitude = lon,
        reason = sprintf("no-coverage: nearest LUCAS sample is %.0f km away (> %.0f km; outside EU survey)", dist[k], max_dist_km))
      next
    }
    rec <- lucas[k, ]
    clay <- rec$clay; sand <- rec$sand
    silt <- if (is.finite(rec$silt)) rec$silt else max(0, 100 - clay - sand)
    om <- if (is.finite(rec$oc)) rec$oc / 10 * 1.724 else 1
    sr <- .lucas_saxton_rawls(sand, clay, om)
    bd <- max(0.9, min(1.8, (1 - sr[["SSAT"]]) * 2.65))
    rows <- lapply(list(c(0, 20), c(20, LUCAS_ROOTING_MAX_CM)), function(d) {
      data.frame(ID = ID, latitude = lat, longitude = lon, depth_top = d[1], depth_bottom = d[2],
                 clay_pct = clay, sand_pct = sand, silt_pct = silt, om_pct = om, bulk_density = bd,
                 SLLL = sr[["SLLL"]], SDUL = sr[["SDUL"]], SSAT = sr[["SSAT"]], stringsAsFactors = FALSE)
    })
    profile_df <- do.call(rbind, rows)
    format_dssat_soil_lucas(profile_df, output_dir_individual)
    results[[length(results) + 1]] <- profile_df
  }

  if (length(results) > 0)
    readr::write_csv(do.call(rbind, results), output_dir_csv, append = file.exists(output_dir_csv))
  if (length(fails) > 0) {
    fail_df <- data.frame(
      ID = vapply(fails, function(f) f$ID, character(1)),
      latitude = vapply(fails, function(f) f$latitude, numeric(1)),
      longitude = vapply(fails, function(f) f$longitude, numeric(1)),
      reason = vapply(fails, function(f) f$reason, character(1)), stringsAsFactors = FALSE)
    failure_log <- file.path(dirname(output_dir_csv),
        paste0(tools::file_path_sans_ext(basename(output_dir_csv)), "_download_failures.csv"))
    readr::write_csv(fail_df, failure_log)
    message(sprintf("[LUCAS] %d of %d point(s) had no LUCAS sample within %.0f km. Details -> %s",
                    nrow(fail_df), nrow(todo), max_dist_km, failure_log))
  }
  message("LUCAS Topsoil Processing Complete.")
  return(TRUE)
}
