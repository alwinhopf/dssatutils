# ==============================================================================
#  SOIL HELPER: POLARIS (30 m probabilistic disaggregation of SSURGO, CONUS)
#  Filename: soil_polaris.R
#  Description: R twin of python/dssatutils/soil_polaris.py. Streams POLARIS
#    GeoTIFF tiles via GDAL /vsicurl (terra), derives DSSAT soil physics from
#    POLARIS's van Genuchten retention curve, and writes per-point .SOL files.
# ==============================================================================

POLARIS_BASE <- "http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0"

.polaris_vars   <- c("clay", "sand", "silt", "bd", "om", "ph",
                     "theta_r", "theta_s", "alpha", "n", "ksat")
.polaris_log10  <- c("om", "alpha", "ksat", "hb")
.polaris_dlabel <- c("0_5", "5_15", "15_30", "30_60", "60_100", "100_200")
.polaris_dbot   <- c(5, 15, 30, 60, 100, 200)
.polaris_dctr   <- c(2.5, 10, 22.5, 45, 80, 150)

PSI_DUL_KPA <- 33
PSI_LL_KPA  <- 1500
.LL_FLOOR <- 0.02
.PAW_MIN  <- 0.04

polaris_tile <- function(lat, lon) {
  s <- floor(lat); w <- floor(lon)
  sprintf("lat%d%d_lon%d%d", s, s + 1L, w, w + 1L)
}

polaris_backtransform <- function(var, value) {
  if (is.na(value)) return(NA_real_)
  if (var %in% .polaris_log10) 10^value else value
}

vg_theta <- function(psi_kpa, theta_r, theta_s, alpha, n) {
  if (is.na(n) || n <= 1 || any(is.na(c(theta_r, theta_s, alpha)))) return(NA_real_)
  m <- 1 - 1 / n
  theta_r + (theta_s - theta_r) / (1 + (alpha * abs(psi_kpa))^n)^m
}

saxton_rawls <- function(sand_pct, clay_pct, om_pct) {
  S <- sand_pct / 100; C <- clay_pct / 100; OM <- om_pct / 100
  t1500 <- -0.024*S + 0.487*C + 0.006*OM + 0.005*S*OM - 0.013*C*OM + 0.068*S*C + 0.031
  lll <- t1500 + (0.14 * t1500 - 0.02)
  t33 <- -0.251*S + 0.195*C + 0.011*OM + 0.006*S*OM - 0.027*C*OM + 0.452*S*C + 0.299
  dul <- t33 + (1.283 * t33^2 - 0.374 * t33 - 0.015)
  ts33t <- 0.278*S + 0.034*C + 0.022*OM - 0.018*S*OM - 0.027*C*OM - 0.584*S*C + 0.078
  ts33 <- ts33t + (0.636 * ts33t - 0.107)
  sat <- dul + ts33 - 0.097 * S + 0.043
  list(SLLL = lll, SDUL = dul, SSAT = sat)
}

water_limits <- function(theta_r, theta_s, alpha, n,
                         sand = NA, clay = NA, om_pct = NA) {
  lll <- vg_theta(PSI_LL_KPA, theta_r, theta_s, alpha, n)
  dul <- vg_theta(PSI_DUL_KPA, theta_r, theta_s, alpha, n)
  sat <- if (is.na(theta_s)) NA_real_ else theta_s
  if (any(is.na(c(lll, dul, sat)))) {
    if (!any(is.na(c(sand, clay, om_pct)))) {
      sr <- saxton_rawls(sand, clay, om_pct)
      if (is.na(lll)) lll <- sr$SLLL
      if (is.na(dul)) dul <- sr$SDUL
      if (is.na(sat)) sat <- sr$SSAT
    } else stop("No usable van Genuchten or texture data for layer.")
  }
  lll <- max(lll, .LL_FLOOR)
  dul <- max(dul, lll + .PAW_MIN)
  sat <- max(sat, dul + .PAW_MIN)
  list(SLLL = lll, SDUL = dul, SSAT = sat)
}

ssks_cmhr <- function(ksat) {
  if (is.na(ksat)) return(-99)
  min(999, max(0, ksat))
}

format_dssat_sol_file_polaris <- function(site_data, output_dir,
                                          source_name = "POLARIS v1.0",
                                          source_tag = "p50") {
  if (nrow(site_data) == 0) stop("No soil layers found for this ID.")
  if (all(is.na(site_data$clay)) || all(is.na(site_data$silt)) || all(is.na(site_data$bd)))
    stop("Critical soil data (clay/silt/bulk density) all NA.")
  soil_id <- as.character(site_data$ID[1])
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  filename <- file.path(output_dir, paste0(soil_id, ".SOL"))
  cat(sprintf("*SOILS: %s\n", source_name), file = filename)
  cat(sprintf("! Source: %s (statistic=%s)\n\n", source_name, source_tag), file = filename, append = TRUE)
  cat(sprintf("*%-10s  POLARIS       %9.3f %9.3f\n", substr(soil_id, 1, 10), site_data$latitude[1], site_data$longitude[1]), file = filename, append = TRUE)
  cat("@SITE        COUNTRY          LAT     LONG SCS FAMILY\n", file = filename, append = TRUE)
  cat(sprintf(" %-11s USA           %9.3f %9.3f \n", substr(soil_id, 1, 11), site_data$latitude[1], site_data$longitude[1]), file = filename, append = TRUE)
  cat("@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE\n", file = filename, append = TRUE)
  cat("    BN   .13     6    .6    73     1     1 IB001 IB001 IB001\n", file = filename, append = TRUE)
  cat("@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC\n", file = filename, append = TRUE)
  site_data <- site_data[order(site_data$depth_bottom), ]
  for (i in seq_len(nrow(site_data))) {
    layer <- site_data[i, ]
    srgf <- exp(-0.02 * layer$depth_center); if (srgf < 0.02) srgf <- 0
    slhw <- if (is.na(layer$ph)) -99 else layer$ph
    cat(sprintf("%6d   -99 %5.3f %5.3f %5.3f %5.2f %5.1f %5.2f %5.2f %5.1f %5.1f   -99   -99 %5.1f   -99   -99   -99\n",
                as.integer(layer$depth_bottom), layer$SLLL, layer$SDUL, layer$SSAT,
                srgf, layer$SSKS, layer$bd, layer$oc_pct, layer$clay, layer$silt, slhw),
        file = filename, append = TRUE)
  }
  cat("\n", file = filename, append = TRUE)
}

polaris_timeout_sec <- function() {
  x <- suppressWarnings(as.numeric(Sys.getenv("POLARIS_TIMEOUT_SEC", "45")))
  if (!is.finite(x) || x <= 0) 45 else x
}
polaris_retries <- function() {
  x <- suppressWarnings(as.integer(Sys.getenv("POLARIS_RETRIES", "3")))
  if (is.na(x) || x < 1L) 3L else x
}
polaris_progress_every <- function() {
  x <- suppressWarnings(as.integer(Sys.getenv("POLARIS_PROGRESS_EVERY", "10")))
  if (is.na(x) || x < 1L) 10L else x
}

polaris_tile_source <- function(var, stat, depth_label, tile, cache_dir = NULL) {
  url <- sprintf("%s/%s/%s/%s/%s.tif", POLARIS_BASE, var, stat, depth_label, tile)
  if (is.null(cache_dir)) return(paste0("/vsicurl/", url))
  local <- file.path(cache_dir, var, stat, depth_label, paste0(tile, ".tif"))
  if (!file.exists(local)) {
    dir.create(dirname(local), recursive = TRUE, showWarnings = FALSE)
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = max(old_timeout, polaris_timeout_sec()))
    tmp <- paste0(local, ".tmp")
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
    if (!file.rename(tmp, local)) stop("Could not move downloaded POLARIS tile into cache: ", local)
  }
  local
}

# Normalize matrix coordinates to SpatVector so terra's ID=FALSE behavior is
# explicit and consistent across terra versions/platforms. The matrix extract
# method has a different formal signature and should not receive ID directly.
polaris_extract_values <- function(src, pts, n_expected, retries = polaris_retries()) {
  last_error <- NULL
  for (attempt in seq_len(retries)) {
    ans <- tryCatch({
      r <- terra::rast(src)
      p <- terra::vect(data.frame(x = pts[, 1], y = pts[, 2]),
                       geom = c("x", "y"), crs = terra::crs(r))
      z <- terra::extract(r, p, ID = FALSE)
      if (is.null(z) || ncol(z) < 1L) stop("terra::extract returned no value column")
      vals <- as.numeric(z[[1]])
      if (length(vals) != n_expected)
        stop(sprintf("terra::extract returned %d values for %d points", length(vals), n_expected))
      vals
    }, error = function(e) { last_error <<- e; NULL })
    if (!is.null(ans)) return(ans)
    if (attempt < retries) Sys.sleep(min(2^(attempt - 1L), 4))
  }
  warning(sprintf("POLARIS raster read failed after %d attempt(s): %s",
                  retries, conditionMessage(last_error)), call. = FALSE)
  rep(NA_real_, n_expected)
}

fetch_polaris <- function(gridfile, id_col, stat, cache_dir = NULL) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("terra is required for POLARIS extraction.")
  gdf <- sf::st_transform(gridfile, 4326)
  coords <- sf::st_coordinates(gdf)
  lons <- coords[, 1]; lats <- coords[, 2]
  ids <- gridfile[[id_col]]
  tile_key <- mapply(polaris_tile, lats, lons)
  tiles <- unique(tile_key)
  total <- length(tiles) * length(.polaris_vars) * length(.polaris_dlabel)
  done <- 0L

  old_gdal_timeout <- Sys.getenv("GDAL_HTTP_TIMEOUT", unset = NA_character_)
  old_gdal_retry <- Sys.getenv("GDAL_HTTP_MAX_RETRY", unset = NA_character_)
  on.exit({
    if (is.na(old_gdal_timeout)) Sys.unsetenv("GDAL_HTTP_TIMEOUT") else Sys.setenv(GDAL_HTTP_TIMEOUT = old_gdal_timeout)
    if (is.na(old_gdal_retry)) Sys.unsetenv("GDAL_HTTP_MAX_RETRY") else Sys.setenv(GDAL_HTTP_MAX_RETRY = old_gdal_retry)
  }, add = TRUE)
  Sys.setenv(GDAL_HTTP_TIMEOUT = as.character(polaris_timeout_sec()), GDAL_HTTP_MAX_RETRY = "1")

  message(sprintf("POLARIS: %d point(s), %d tile(s), %d remote raster read(s); timeout=%ss, retries=%d",
                  length(ids), length(tiles), total, polaris_timeout_sec(), polaris_retries()))
  out <- vector("list", total)
  out_i <- 0L
  for (tile_i in seq_along(tiles)) {
    tile <- tiles[[tile_i]]
    sel <- which(tile_key == tile)
    pts <- cbind(lons[sel], lats[sel])
    message(sprintf("POLARIS tile %d/%d: %s (%d point%s)", tile_i, length(tiles), tile,
                    length(sel), if (length(sel) == 1L) "" else "s"))
    for (var in .polaris_vars) {
      for (d in seq_along(.polaris_dlabel)) {
        src <- polaris_tile_source(var, stat, .polaris_dlabel[d], tile, cache_dir)
        vals <- polaris_extract_values(src, pts, length(sel))
        vals <- vapply(vals, function(v) polaris_backtransform(var, v), numeric(1))
        out_i <- out_i + 1L
        out[[out_i]] <- data.frame(ID = ids[sel], prop = var,
                                   depth_bottom = .polaris_dbot[d], depth_center = .polaris_dctr[d],
                                   value = vals, stringsAsFactors = FALSE)
        done <- done + 1L
        if (done %% polaris_progress_every() == 0L || done == total)
          message(sprintf("POLARIS progress: %d/%d raster reads (%.1f%%)", done, total, 100 * done / total))
      }
    }
  }
  do.call(rbind, out[seq_len(out_i)])
}

process_soils_polaris <- function(gridfile, soilfile_csv_path, output_sol_dir,
                                  id_col, stat = "p50", cache_dir = NULL) {
  if (!stat %in% c("p50", "mean", "mode", "p5", "p95")) stop(sprintf("Unknown POLARIS statistic: %s", stat))
  message(sprintf("--- POLARIS extraction (statistic=%s, CONUS 30 m) ---", stat))
  grid_wgs84 <- sf::st_transform(gridfile, 4326)
  cc <- sf::st_coordinates(grid_wgs84)
  long_df <- fetch_polaris(gridfile, id_col, stat, cache_dir)
  if (is.null(long_df) || nrow(long_df) == 0) stop("No POLARIS data extracted (coords outside CONUS, or server unreachable).")

  wide <- as.data.frame(tidyr::pivot_wider(long_df,
    id_cols = c("ID", "depth_bottom", "depth_center"), names_from = "prop", values_from = "value",
    values_fn = function(x) x[1]))
  for (v in .polaris_vars) if (!v %in% names(wide)) wide[[v]] <- NA_real_
  usable_ids <- unique(wide$ID[!is.na(wide$sand) | !is.na(wide$clay)])
  wide <- wide[wide$ID %in% usable_ids, , drop = FALSE]
  if (nrow(wide) == 0) stop("No usable POLARIS data for any point.")

  wide$oc_pct <- wide$om / 1.724
  wide$SSKS <- vapply(wide$ksat, ssks_cmhr, numeric(1))
  lim <- mapply(function(tr, ts, al, nn, sa, cl, om)
    water_limits(tr, ts, al, nn, sand = sa, clay = cl, om_pct = om),
    wide$theta_r, wide$theta_s, wide$alpha, wide$n, wide$sand, wide$clay, wide$om,
    SIMPLIFY = FALSE)
  wide$SLLL <- vapply(lim, `[[`, numeric(1), "SLLL")
  wide$SDUL <- vapply(lim, `[[`, numeric(1), "SDUL")
  wide$SSAT <- vapply(lim, `[[`, numeric(1), "SSAT")

  coords_df <- data.frame(ID = grid_wgs84[[id_col]], longitude = cc[, 1], latitude = cc[, 2], stringsAsFactors = FALSE)
  final_df <- merge(wide, coords_df, by = "ID")
  utils::write.csv(data.frame(ID = grid_wgs84[[id_col]], SOIL_ID = grid_wgs84[[id_col]]), soilfile_csv_path, row.names = FALSE)
  if (!dir.exists(output_sol_dir)) dir.create(output_sol_dir, recursive = TRUE)
  log_path <- file.path(output_sol_dir, "soil_processing_errors.log")
  cat(sprintf("Log started: %s\n", Sys.time()), file = log_path)

  success <- 0; errors <- 0
  for (uid in unique(final_df$ID)) {
    subset_df <- final_df[final_df$ID == uid, , drop = FALSE]
    res <- tryCatch({ format_dssat_sol_file_polaris(subset_df, output_sol_dir, source_tag = stat); TRUE },
                    error = function(e) { cat(sprintf("ID: %s | Error: %s\n", uid, conditionMessage(e)), file = log_path, append = TRUE); FALSE })
    if (isTRUE(res)) success <- success + 1 else errors <- errors + 1
  }
  message(sprintf("POLARIS processing complete. Success: %d, Errors: %d", success, errors))
  invisible(NULL)
}
