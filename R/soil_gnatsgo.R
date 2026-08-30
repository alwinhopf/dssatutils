# ==============================================================================
#  SOIL HELPER: gNATSGO (USDA gridded National Soil Survey) -- SMART RESUME
#  Filename: soil_gnatsgo.R   (R twin of python/dssatutils/soil_gnatsgo.py)
# ==============================================================================
#
#  WHY gNATSGO over SSURGO: gNATSGO is a complete, gap-free 30 m soil grid for the
#  conterminous US (SSURGO + STATSGO2 + Raster Soil Surveys). It has a map unit
#  EVERYWHERE there is land, filling the "no-coverage" holes plain SSURGO leaves.
#  Where the underlying map unit is the SSURGO one (most cropland) the profile is
#  IDENTICAL to SSURGO -- this module reuses soil_ssurgo.R's tabular query, soil
#  physics (Saxton & Rawls) and .SOL layout so the two are byte-comparable there.
#
#  ACCESS (two public services, no key):
#   1. Map-unit key at a point -- SoilWeb Web Coverage Service via soilDB::mukey.wcs
#      (same 30 m gNATSGO mukey grid the Python module fetches over raw WCS).
#   2. Tabular soil properties for that mukey -- USDA Soil Data Access (SDA), the
#      same endpoint soil_ssurgo.R uses (SDA hosts SSURGO + STATSGO2 tabular).
#  A few gNATSGO Raster-Soil-Survey mukeys are not in SDA; those are logged as
#  "no-tabular".
#
#  Depends on soil_ssurgo.R being sourced first (robust_SDA_query,
#  calculate_soil_properties) -- guaranteed inside the package namespace.
# ==============================================================================


# --- gNATSGO mukey at a WGS84 point via the SoilWeb WCS -----------------------
# Returns list(ok, mukey, error). mukey = NA means "no map unit here" (water /
# outside CONUS). Uses soilDB::mukey.wcs (terra SpatRaster) + a tiny AOI box.
gnatsgo_mukey_at_point <- function(lat, lon, buffer_m = 45, max_retries = 3,
                                   retry_delay_seconds = 5) {
  if (!requireNamespace("terra", quietly = TRUE))
    return(list(ok = FALSE, mukey = NA_character_, error = "package 'terra' not installed"))

  pt5070 <- tryCatch(
    sf::st_transform(sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326), 5070),
    error = function(e) e)
  if (inherits(pt5070, "error"))
    return(list(ok = FALSE, mukey = NA_character_, error = sprintf("reprojection failed: %s", conditionMessage(pt5070))))
  xy <- sf::st_coordinates(pt5070)
  aoi <- sf::st_as_sfc(sf::st_bbox(c(xmin = xy[1] - buffer_m, ymin = xy[2] - buffer_m,
                                     xmax = xy[1] + buffer_m, ymax = xy[2] + buffer_m),
                                   crs = sf::st_crs(5070)))
  last_err <- NA_character_
  for (attempt in seq_len(max_retries)) {
    r <- tryCatch(soilDB::mukey.wcs(aoi = aoi, db = "gnatsgo", res = 30),
                  error = function(e) e)
    if (!inherits(r, "error")) {
      # mukey.wcs returns a CATEGORICAL raster: the cell stores a factor whose
      # LABEL is the mukey string. as.character() yields the label; as.integer()
      # would wrongly return the factor level code. Extract at the exact point.
      val <- tryCatch(terra::extract(r, terra::vect(pt5070))[1, 2],
                      error = function(e) NA)
      mukey_chr <- if (is.na(val) || is.null(val)) NA_character_ else as.character(val)
      mukey <- suppressWarnings(as.integer(mukey_chr))
      if (is.na(mukey) || mukey <= 0)
        return(list(ok = TRUE, mukey = NA_character_, error = NA_character_))
      return(list(ok = TRUE, mukey = as.character(mukey), error = NA_character_))
    }
    last_err <- trimws(conditionMessage(r))
    if (attempt < max_retries) Sys.sleep(retry_delay_seconds)
  }
  list(ok = FALSE, mukey = NA_character_, error = last_err)
}


# --- DSSAT .SOL writer (gNATSGO-labelled; identical column layout to SSURGO) ---
format_dssat_soil_gnatsgo <- function(profile_data, output_dir) {
  soil_id <- as.character(profile_data$ID[1])
  filename <- file.path(output_dir, paste0(soil_id, ".SOL"))
  if (file.exists(filename)) return()

  cat("*SOILS: USA gNATSGO Soil Profiles\n", file = filename)
  cat("! Generated from gNATSGO (SoilWeb WCS mukey grid + USDA SDA tabular)\n\n",
      file = filename, append = TRUE)
  cat(sprintf("*%-6s  gNATSGO       %9.3f %9.3f\n",
              soil_id, profile_data$latitude[1], profile_data$longitude[1]),
      file = filename, append = TRUE)
  cat("@SITE        COUNTRY          LAT     LONG SCS FAMILY\n", file = filename, append = TRUE)
  cat(sprintf(" %-11s USA         %9.3f %9.3f \n",
              soil_id, profile_data$latitude[1], profile_data$longitude[1]),
      file = filename, append = TRUE)
  cat("@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE\n", file = filename, append = TRUE)
  cat("    BN   .13     6    .6    73     1     1 IB001 IB001 IB001\n", file = filename, append = TRUE)
  cat("@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC\n",
      file = filename, append = TRUE)

  profile_data %>%
    dplyr::mutate(depth_num = as.numeric(sub(".*-", "", sub("cm", "", depth_range)))) %>%
    dplyr::arrange(depth_num) %>%
    dplyr::group_walk(function(layer, key) {
      depth_val <- layer$depth_num
      slll <- sub("^0", " ", sprintf("%5.3f", layer$SLLL))
      sdul <- sub("^0", " ", sprintf("%5.3f", layer$SDUL))
      ssat <- sub("^0", " ", sprintf("%5.3f", layer$SSAT))
      depth_format <- sprintf("%6d", depth_val)
      ssks_val <- if ("SSKS" %in% names(layer)) layer$SSKS else rep(NA_real_, nrow(layer))
      ssks_str <- ifelse(!is.na(ssks_val) & ssks_val > 0,
                         sprintf("%6.2f", pmin(999.0, ssks_val)),
                         "   -99")
      cat(paste0(sprintf("%s   -99 %s %s %s  1.00%s %5.2f %5.2f %5.1f %5.1f   -99   -99   -99   -99   -99   -99\n",
                         depth_format, slll, sdul, ssat, ssks_str,
                         layer$bulk_density, layer$om_pct / 1.724, layer$clay_pct, layer$silt_pct),
                 collapse = ""),
          file = filename, append = TRUE)
    })
  cat("\n", file = filename, append = TRUE)
}


# --- Main processing function -------------------------------------------------
# Signature mirrors process_soils_ssurgo exactly.
process_soils_gnatsgo <- function(grid_points, output_dir_csv, output_dir_individual, n_cores,
                                  id_col, lat_col, long_col, format_sql_func = NULL) {

  message("Starting gNATSGO Processing (Smart Resume Mode)...")
  dir.create(output_dir_individual, recursive = TRUE, showWarnings = FALSE)

  if (missing(format_sql_func) || is.null(format_sql_func) || !is.function(format_sql_func)) {
    format_sql_func <- function(vec) {
      vals <- unique(stats::na.omit(as.character(vec)))
      if (length(vals) == 0) return("('')")
      paste0("('", paste(gsub("'", "''", vals), collapse = "','"), "')")
    }
  }

  grid_df <- grid_points %>% sf::st_drop_geometry() %>% as.data.frame()
  if (!lat_col %in% names(grid_df))  grid_df[[lat_col]]  <- sf::st_coordinates(grid_points)[, 2]
  if (!long_col %in% names(grid_df)) grid_df[[long_col]] <- sf::st_coordinates(grid_points)[, 1]

  if (dir.exists(output_dir_individual)) {
    existing_files <- tools::file_path_sans_ext(list.files(output_dir_individual, pattern = "\\.SOL$"))
    missing_mask <- !(as.character(grid_df[[id_col]]) %in% existing_files)
    points_to_process <- grid_df[missing_mask, ]
    n_total <- nrow(grid_df); n_process <- nrow(points_to_process)
    message(sprintf("Resume Check: Found %d existing profiles. Processing %d remaining points.",
                    n_total - n_process, n_process))
    if (n_process == 0) {
      message("All soil profiles already exist. Skipping gNATSGO processing.")
      return(TRUE)
    }
  } else {
    points_to_process <- grid_df
    message(sprintf("First Run: Processing all %d points.", nrow(points_to_process)))
  }

  # --- Worker: WCS mukey lookup, then SSURGO-identical tabular path ----------
  process_point_wrapper <- function(point_data_row) {
    ID   <- as.character(point_data_row[[id_col]])
    LATv <- point_data_row[[lat_col]]
    LONv <- point_data_row[[long_col]]
    fail <- function(reason) list(.fail = TRUE, ID = ID,
                                  latitude = LATv, longitude = LONv, reason = reason)

    # 1. gNATSGO mukey at the point (WCS grid)
    mk <- gnatsgo_mukey_at_point(LATv, LONv)
    if (!isTRUE(mk$ok))
      return(fail(sprintf("network: gNATSGO WCS lookup failed after retries (%s)",
                          if (is.na(mk$error)) "no detail" else mk$error)))
    if (is.na(mk$mukey))
      return(fail("no-coverage: no gNATSGO map unit at this location (water / outside CONUS)"))
    mukey <- mk$mukey

    # 2. Bedrock depth
    q_bedrock <- sprintf("SELECT mukey, brockdepmin FROM muaggatt WHERE mukey IN %s",
                         format_sql_func(mukey))
    bq <- robust_SDA_query(q_bedrock)
    bedrock_depth <- 200
    if (isTRUE(bq$ok) && !is.null(bq$data)) {
      bd <- as.data.frame(bq$data)
      if (nrow(bd) > 0 && !all(is.na(bd$brockdepmin))) bedrock_depth <- min(bd$brockdepmin, na.rm = TRUE)
    }
    if (is.infinite(bedrock_depth)) bedrock_depth <- 200

    all_layers <- list("0-5cm"=c(0,5), "5-20cm"=c(5,20), "20-35cm"=c(20,35),
                       "35-50cm"=c(35,50), "50-65cm"=c(50,65), "65-80cm"=c(65,80),
                       "80-95cm"=c(80,95), "95-110cm"=c(95,110), "110-125cm"=c(110,125),
                       "125-140cm"=c(125,140), "140-155cm"=c(140,155), "155-170cm"=c(155,170),
                       "170-185cm"=c(170,185), "185-200cm"=c(185,200))
    valid_layers <- all_layers[sapply(all_layers, function(x) x[1] < bedrock_depth)]
    if (length(valid_layers) > 0) {
      last <- length(valid_layers)
      if (valid_layers[[last]][2] > bedrock_depth) valid_layers[[last]][2] <- bedrock_depth
      names(valid_layers)[last] <- paste0(valid_layers[[last]][1], "-", valid_layers[[last]][2], "cm")
    } else {
      valid_layers <- list(); valid_layers[[paste0("0-", bedrock_depth, "cm")]] <- c(0, bedrock_depth)
    }

    # 3. Horizon properties (SDA)
    q_soil <- sprintf("SELECT component.mukey, component.cokey, component.comppct_r,
                       chorizon.hzdept_r, chorizon.hzdepb_r, chorizon.claytotal_r,
                       chorizon.sandtotal_r, chorizon.om_r, chorizon.dbthirdbar_r
                       FROM component INNER JOIN chorizon ON component.cokey = chorizon.cokey
                       WHERE component.mukey IN %s", format_sql_func(mukey))
    pq <- robust_SDA_query(q_soil)
    if (!isTRUE(pq$ok))
      return(fail(sprintf("network: SDA horizon query failed after retries (%s)",
                          if (is.na(pq$error)) "no detail" else pq$error)))
    props <- as.data.frame(pq$data)
    if (nrow(props) == 0) {
      # Distinguish a non-soil unit from a Raster-Soil-Survey mukey absent in SDA.
      nm <- robust_SDA_query(sprintf("SELECT mukey, muname FROM mapunit WHERE mukey IN %s",
                                     format_sql_func(mukey)))
      has_mu <- isTRUE(nm$ok) && !is.null(nm$data) && nrow(as.data.frame(nm$data)) > 0
      if (!has_mu)
        return(fail(sprintf("no-tabular: gNATSGO mukey %s not found in SDA (likely a Raster Soil Survey unit)", mukey)))
      muname <- paste(unique(stats::na.omit(as.data.frame(nm$data)$muname)), collapse = "; ")
      return(fail(sprintf("no-soil: map unit has no soil horizons%s -- typically Water / Urban / Pits / Rock outcrop (mukey %s)",
                          if (nzchar(muname)) sprintf(" [%s]", muname) else "", mukey)))
    }

    results_list <- lapply(names(valid_layers), function(layer_name) {
      d <- valid_layers[[layer_name]]
      cp <- calculate_soil_properties(props, d[1], d[2])
      if (nrow(cp) > 0) cp$depth_range <- layer_name
      cp
    })
    results_list <- results_list[sapply(results_list, function(x) !is.null(x) && nrow(x) > 0)]
    if (length(results_list) == 0)
      return(fail(sprintf("no-layers: horizon data present but no usable layers after depth filtering (bedrock %s cm)", bedrock_depth)))

    results_df <- do.call(rbind, results_list)
    results_df <- results_df %>% dplyr::mutate(
      ID = ID, longitude = LONv, latitude = LATv, bedrock_depth_cm = bedrock_depth,
      sand_dec = sand_pct/100, clay_dec = clay_pct/100, om_dec = om_pct/100,
      theta_1500t = -0.024*sand_dec + 0.487*clay_dec + 0.006*om_dec + 0.005*(sand_dec*om_dec) - 0.013*(clay_dec*om_dec) + 0.068*(sand_dec*clay_dec) + 0.031,
      SLLL_raw = theta_1500t + (0.14*theta_1500t - 0.02),
      # Floor the wilting point (see soil_ssurgo.R): Saxton & Rawls yields SLLL<=0
      # on very sandy layers, which SIGFPEs DSSAT. Clamp to 0.02 cm3/cm3; SLLL_raw
      # is kept so the caller can tally and log the adjustment.
      SLLL = pmax(SLLL_raw, 0.02),
      theta_33t = -0.251*sand_dec + 0.195*clay_dec + 0.011*om_dec + 0.006*(sand_dec*om_dec) - 0.027*(clay_dec*om_dec) + 0.452*(sand_dec*clay_dec) + 0.299,
      SDUL_raw = theta_33t + (1.283*(theta_33t)^2 - 0.374*theta_33t - 0.015),
      # Floor field capacity so plant-available water (DUL-LL) stays usable (see
      # soil_ssurgo.R): near-pure-sand layers get SDUL barely above the floored
      # SLLL, which SIGFPEs DSSAT's water balance. Guarantee DUL-LL >= 0.04;
      # SDUL_raw is kept so the caller can tally and log the adjustment.
      SDUL = pmax(SDUL_raw, SLLL + 0.04),
      theta_s33t = 0.278*sand_dec + 0.034*clay_dec + 0.022*om_dec - 0.018*(sand_dec*om_dec) - 0.027*(clay_dec*om_dec) - 0.584*(sand_dec*clay_dec) + 0.078,
      theta_s33 = theta_s33t + (0.636*theta_s33t - 0.107),
      SSAT = SDUL + theta_s33 - 0.097*sand_dec + 0.043,
      SSKS = .saxton_rawls_ssks(SSAT, SDUL, SLLL, coarse_fraction = 0, bulk_density = bulk_density)
    )
    format_dssat_soil_gnatsgo(results_df, output_dir_individual)
    results_df
  }

  total_points <- nrow(points_to_process)
  CHUNK_SIZE <- 10000
  num_chunks <- ceiling(total_points / CHUNK_SIZE)
  message(sprintf("Processing batch of %d points in %d chunks...", total_points, num_chunks))

  # WCS + SDA are network-bound; cap concurrency below the SSURGO ceiling.
  n_cores <- max(1L, min(n_cores, 8L))
  if (.Platform$OS.type == "windows" && n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, { library(soilDB); library(sf); library(terra); library(dplyr) })
    parallel::clusterExport(cl, varlist = c("process_point_wrapper", "gnatsgo_mukey_at_point",
                                            ".saxton_rawls_ssks",
                                            "robust_SDA_query", "calculate_soil_properties",
                                            "format_dssat_soil_gnatsgo", "output_dir_individual",
                                            "id_col", "lat_col", "long_col", "format_sql_func"),
                            envir = environment())
  } else {
    cl <- if (n_cores > 1L) n_cores else NULL
  }

  all_fails <- list()
  for (i in 1:num_chunks) {
    start_idx <- (i - 1) * CHUNK_SIZE + 1
    end_idx <- min(i * CHUNK_SIZE, total_points)
    message(sprintf("  > Chunk %d/%d (Points %d - %d)", i, num_chunks, start_idx, end_idx))
    chunk_data <- points_to_process[start_idx:end_idx, ]
    chunk_list <- split(chunk_data, seq(nrow(chunk_data)))
    chunk_results <- pbapply::pblapply(chunk_list, process_point_wrapper, cl = cl)

    is_fail <- function(x) is.list(x) && isTRUE(x[[".fail"]])
    all_fails     <- c(all_fails, Filter(is_fail, chunk_results))
    valid_results <- Filter(function(x) !is.null(x) && !is_fail(x), chunk_results)
    if (length(valid_results) > 0) {
      chunk_df <- dplyr::bind_rows(valid_results)
      # Surface wilting-point clamps in the run log, then drop the helper column.
      if ("SLLL_raw" %in% names(chunk_df)) {
        clamped <- chunk_df[!is.na(chunk_df$SLLL_raw) & chunk_df$SLLL_raw < 0.02, ]
        if (nrow(clamped) > 0) {
          ids <- unique(as.character(clamped$ID))
          message(sprintf("[gNATSGO] SLLL floored to 0.020 on %d layer(s) across %d point(s) (negative/low Saxton-Rawls wilting point; sandy soils): %s",
                          nrow(clamped), length(ids), paste(ids, collapse = ", ")))
        }
        chunk_df$SLLL_raw <- NULL
      }
      # Surface field-capacity raises (minimum available-water enforcement) too.
      if ("SDUL_raw" %in% names(chunk_df)) {
        raised <- chunk_df[!is.na(chunk_df$SDUL_raw) & (chunk_df$SDUL - chunk_df$SDUL_raw) > 1e-9, ]
        if (nrow(raised) > 0) {
          ids <- unique(as.character(raised$ID))
          message(sprintf("[gNATSGO] SDUL raised to keep DUL-LL >= 0.040 on %d layer(s) across %d point(s) (near-zero plant-available water; sandy soils): %s",
                          nrow(raised), length(ids), paste(ids, collapse = ", ")))
        }
        chunk_df$SDUL_raw <- NULL
      }
      readr::write_csv(chunk_df, output_dir_csv, append = file.exists(output_dir_csv))
    }
    rm(chunk_list, chunk_results, valid_results); if (exists("chunk_df")) rm(chunk_df); gc()
  }
  if (.Platform$OS.type == "windows" && inherits(cl, "cluster")) parallel::stopCluster(cl)

  if (length(all_fails) > 0) {
    fail_df <- data.frame(
      ID        = vapply(all_fails, function(f) as.character(f$ID), character(1)),
      latitude  = vapply(all_fails, function(f) suppressWarnings(as.numeric(f$latitude)),  numeric(1)),
      longitude = vapply(all_fails, function(f) suppressWarnings(as.numeric(f$longitude)), numeric(1)),
      reason    = vapply(all_fails, function(f) as.character(f$reason), character(1)),
      stringsAsFactors = FALSE)
    failure_log <- file.path(dirname(output_dir_csv),
        paste0(tools::file_path_sans_ext(basename(output_dir_csv)), "_download_failures.csv"))
    readr::write_csv(fail_df, failure_log)
    message(sprintf("[gNATSGO] %d of %d processed point(s) produced NO soil profile:", nrow(fail_df), total_points))
    cats <- sub(":.*$", "", fail_df$reason)
    tb <- sort(table(cats), decreasing = TRUE)
    for (nm in names(tb)) message(sprintf("   - %-13s %d   (%s)", paste0(nm, ":"), tb[[nm]],
        switch(nm,
          "network"     = "transient WCS/SDA/timeout -- re-run to retry these",
          "no-coverage" = "no gNATSGO map unit here (water / outside CONUS)",
          "no-soil"     = "non-soil map unit (Water, Urban, Pits, Rock) -- no horizons exist",
          "no-tabular"  = "gNATSGO mukey not in SDA (Raster Soil Survey unit)",
          "no-layers"   = "horizons present but unusable after depth filtering",
          "other")))
    message(sprintf("   Per-point details (ID, lat, long, reason) -> %s", failure_log))
  } else {
    message("[gNATSGO] All processed points produced a soil profile.")
  }
  message("gNATSGO Processing Complete.")
  return(TRUE)
}
