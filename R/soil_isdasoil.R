# ==============================================================================
#  SOIL HELPER: iSDAsoil (Africa, 30 m) — DSSAT .SOL
#  Filename: soil_isdasoil.R   (R twin of python/dssatutils/soil_isdasoil.py)
# ==============================================================================
#
#  iSDAsoil: 30 m predicted soil properties for AFRICA at two depths (0-20 cm,
#  20-50 cm) — the highest-resolution open soil product for the continent.
#
#  ACCESS (fully open, NO key): cloud-optimised GeoTIFFs in a public S3 bucket,
#  streamed per-point via GDAL /vsicurl (terra). Each property file has 4 bands:
#    band 1 = mean 0-20 cm, band 2 = mean 20-50 cm, bands 3-4 = std-dev.
#  Native CRS EPSG:3857, nodata = 255.
#    https://isdasoil.s3.amazonaws.com/soil_data/<property>/<property>.tif
#
#  Stored uint8 -> real units (verified against the iSDA / GEE catalogue):
#    clay_content / sand_content : value as-is (%)
#    carbon_organic              : exp(x/10) - 1  (g/kg) -> OC% -> OM%
#    bulk_density                : x / 100         (g/cm^3)
#
#  DSSAT physics (Saxton & Rawls) + .SOL layout match soil_ssurgo.R, so an
#  iSDAsoil profile is comparable to a SSURGO/gNATSGO one. Coverage: Africa.
# ==============================================================================

ISDA_BASE <- "https://isdasoil.s3.amazonaws.com/soil_data"
ISDA_NODATA <- 255
# DSSAT layers from the two iSDA depths; the 20-50 cm value is carried to
# ISDA_ROOTING_MAX_CM so DSSAT has a usable rooting profile.
ISDA_ROOTING_MAX_CM <- 150


.isda_back_transform <- function(prop, raw) {
  v <- as.numeric(raw)
  v[raw == ISDA_NODATA] <- NA_real_
  if (prop == "carbon_organic") return(exp(v / 10) - 1)   # g/kg
  if (prop == "bulk_density")   return(v / 100)            # g/cm^3
  if (prop == "ph")             return(v / 10)
  v                                                        # clay/sand: %
}

# Sample the two mean bands (0-20, 20-50 cm) of one property at points (EPSG:3857
# SpatVector). Returns an (n x 2) matrix of real-unit values, or NULL on error.
.isda_sample_property <- function(prop, pts_vect) {
  url <- sprintf("/vsicurl/%s/%s/%s.tif", ISDA_BASE, prop, prop)
  r <- tryCatch(terra::rast(url), error = function(e) e)
  if (inherits(r, "error")) {
    warning(sprintf("iSDAsoil: could not read %s COG: %s", prop, conditionMessage(r)))
    return(NULL)
  }
  ex <- tryCatch(terra::extract(r[[1:2]], pts_vect, ID = FALSE),
                 error = function(e) e)
  if (inherits(ex, "error")) {
    warning(sprintf("iSDAsoil: extract failed for %s: %s", prop, conditionMessage(ex)))
    return(NULL)
  }
  cbind(.isda_back_transform(prop, ex[[1]]), .isda_back_transform(prop, ex[[2]]))
}


# --- DSSAT .SOL writer (iSDAsoil-labelled; same column layout as SSURGO) -------
format_dssat_soil_isdasoil <- function(profile_data, output_dir) {
  soil_id <- as.character(profile_data$ID[1])
  filename <- file.path(output_dir, paste0(soil_id, ".SOL"))
  if (file.exists(filename)) return()
  cat("*SOILS: Africa iSDAsoil Soil Profiles\n", file = filename)
  cat("! Generated from iSDAsoil 30 m (0-20 & 20-50 cm), Saxton & Rawls physics\n\n",
      file = filename, append = TRUE)
  cat(sprintf("*%-6s  ISDA          %9.3f %9.3f\n",
              soil_id, profile_data$latitude[1], profile_data$longitude[1]),
      file = filename, append = TRUE)
  cat("@SITE        COUNTRY          LAT     LONG SCS FAMILY\n", file = filename, append = TRUE)
  cat(sprintf(" %-11s -99         %9.3f %9.3f \n",
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
    ssks_str <- if ("SSKS" %in% names(layer) && !is.na(layer$SSKS) && layer$SSKS > 0) {
      sprintf("%6.2f", min(999.0, layer$SSKS))
    } else {
      "   -99"
    }
    cat(sprintf("%6d   -99 %s %s %s  1.00%s %5.2f %5.2f %5.1f %5.1f   -99   -99   -99   -99   -99   -99\n",
                as.integer(layer$depth_bottom), slll, sdul, ssat, ssks_str,
                layer$bulk_density, layer$om_pct / 1.724, layer$clay_pct, layer$silt_pct),
        file = filename, append = TRUE)
  }
  cat("\n", file = filename, append = TRUE)
}


# --- Saxton & Rawls (2006) — identical to soil_ssurgo.R's inline formulas ------
.isda_saxton_rawls <- function(sand_pct, clay_pct, om_pct, bd = 1.4) {
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
  SSKS  <- .saxton_rawls_ssks(SSAT, SDUL, SLLL, coarse_fraction = 0, bulk_density = bd)
  c(SLLL = SLLL, SDUL = SDUL, SSAT = SSAT, SSKS = SSKS)
}


# --- Public entry point — signature mirrors process_soils_ssurgo --------------
process_soils_isdasoil <- function(grid_points, output_dir_csv, output_dir_individual, n_cores,
                                   id_col, lat_col, long_col, format_sql_func = NULL) {
  message("Starting iSDAsoil Processing (Smart Resume Mode)...")
  dir.create(output_dir_individual, recursive = TRUE, showWarnings = FALSE)

  pts <- sf::st_transform(grid_points, 4326)
  grid_df <- sf::st_drop_geometry(pts)
  xy <- sf::st_coordinates(pts)
  grid_df[[lat_col]] <- xy[, 2]; grid_df[[long_col]] <- xy[, 1]

  existing <- tools::file_path_sans_ext(list.files(output_dir_individual, pattern = "\\.SOL$"))
  todo_mask <- !(as.character(grid_df[[id_col]]) %in% existing)
  todo <- grid_df[todo_mask, , drop = FALSE]
  message(sprintf("Resume Check: Found %d existing profiles. Processing %d remaining points.",
                  nrow(grid_df) - nrow(todo), nrow(todo)))
  if (nrow(todo) == 0) {
    message("All soil profiles already exist. Skipping iSDAsoil processing.")
    return(TRUE)
  }

  ids  <- as.character(todo[[id_col]])
  lats <- as.numeric(todo[[lat_col]]); lons <- as.numeric(todo[[long_col]])
  pts3857 <- terra::project(terra::vect(cbind(lons, lats), type = "points", crs = "EPSG:4326"),
                            "EPSG:3857")

  props <- list()
  for (prop in c("clay_content", "sand_content", "carbon_organic", "bulk_density")) {
    props[[prop]] <- .isda_sample_property(prop, pts3857)
  }
  if (is.null(props$clay_content) || is.null(props$sand_content)) {
    message("iSDAsoil: clay/sand COGs unreadable; aborting (check network / S3 access).")
    return(FALSE)
  }

  layer_spec <- list(c(0, 20, 1), c(20, 50, 2), c(50, ISDA_ROOTING_MAX_CM, 2))  # top, bottom, band-col
  results <- list(); fails <- list()

  for (i in seq_along(ids)) {
    ID <- ids[i]
    if (file.exists(file.path(output_dir_individual, paste0(ID, ".SOL")))) next
    if (!is.finite(props$clay_content[i, 1])) {
      fails[[length(fails) + 1]] <- list(.fail = TRUE, ID = ID, latitude = lats[i], longitude = lons[i],
        reason = "no-coverage: no iSDAsoil value at this location (outside Africa / water / nodata)")
      next
    }
    rows <- lapply(layer_spec, function(sp) {
      top <- sp[1]; bottom <- sp[2]; b <- sp[3]
      clay <- props$clay_content[i, b]; sand <- props$sand_content[i, b]
      oc <- if (!is.null(props$carbon_organic)) props$carbon_organic[i, b] else NA
      bd <- if (!is.null(props$bulk_density)) props$bulk_density[i, b] else NA
      clay <- if (is.finite(clay)) clay else 20
      sand <- if (is.finite(sand)) sand else 40
      om <- if (is.finite(oc)) oc / 10 * 1.724 else 1     # g/kg -> OC% -> OM%
      bd <- if (is.finite(bd)) bd else 1.4
      silt <- max(0, 100 - clay - sand)
      sr <- .isda_saxton_rawls(sand, clay, om, bd = bd)
      data.frame(ID = ID, latitude = lats[i], longitude = lons[i],
                 depth_top = top, depth_bottom = bottom,
                 clay_pct = clay, sand_pct = sand, silt_pct = silt,
                 om_pct = om, bulk_density = bd,
                 SLLL = sr[["SLLL"]], SDUL = sr[["SDUL"]], SSAT = sr[["SSAT"]],
                 SSKS = sr[["SSKS"]],
                 stringsAsFactors = FALSE)
    })
    profile_df <- do.call(rbind, rows)
    format_dssat_soil_isdasoil(profile_df, output_dir_individual)
    results[[length(results) + 1]] <- profile_df
  }

  if (length(results) > 0) {
    out <- do.call(rbind, results)
    readr::write_csv(out, output_dir_csv, append = file.exists(output_dir_csv))
  }
  if (length(fails) > 0) {
    fail_df <- data.frame(
      ID = vapply(fails, function(f) f$ID, character(1)),
      latitude = vapply(fails, function(f) f$latitude, numeric(1)),
      longitude = vapply(fails, function(f) f$longitude, numeric(1)),
      reason = vapply(fails, function(f) f$reason, character(1)), stringsAsFactors = FALSE)
    failure_log <- file.path(dirname(output_dir_csv),
        paste0(tools::file_path_sans_ext(basename(output_dir_csv)), "_download_failures.csv"))
    readr::write_csv(fail_df, failure_log)
    message(sprintf("[iSDAsoil] %d of %d point(s) produced NO soil profile (no-coverage). Details -> %s",
                    nrow(fail_df), length(ids), failure_log))
  }
  message("iSDAsoil Processing Complete.")
  return(TRUE)
}
