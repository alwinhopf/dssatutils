# ==============================================================================
#  SOIL HELPER: SSURGO (USDA WEB SERVICE) - SMART RESUME & MEMORY OPTIMIZED
#  Filename: soil_ssurgo.R
# ==============================================================================



# --- Wrappers for Retries ---
# Return a STRUCTURED result: list(ok, data, error). Previously these returned
# NULL on failure, discarding the actual error (network vs. server vs. bad SQL),
# which made per-point failures impossible to diagnose. Callers can now report
# the real reason a point produced no soil profile.
robust_SDA_query <- function(query, max_retries = 3, retry_delay_seconds = 5) {
  last_err <- NA_character_
  for (attempt in 1:max_retries) {
    result <- try(soilDB::SDA_query(query), silent = TRUE)
    if (!inherits(result, "try-error")) return(list(ok = TRUE, data = result, error = NA_character_))
    last_err <- trimws(as.character(result))
    Sys.sleep(retry_delay_seconds)
  }
  list(ok = FALSE, data = NULL, error = last_err)
}

robust_SDA_spatialQuery <- function(point_sf, what, max_retries = 3, retry_delay_seconds = 5) {
  last_err <- NA_character_
  for (attempt in 1:max_retries) {
    result <- try(soilDB::SDA_spatialQuery(point_sf, what = what), silent = TRUE)
    if (!inherits(result, "try-error")) return(list(ok = TRUE, data = result, error = NA_character_))
    last_err <- trimws(as.character(result))
    Sys.sleep(retry_delay_seconds)
  }
  list(ok = FALSE, data = NULL, error = last_err)
}

# --- Calculation Logic ---
calculate_soil_properties <- function(soil_properties, top_depth, bottom_depth) {
  ungrouped <- soil_properties %>%
    dplyr::mutate(
      adj_top = pmax(hzdept_r, top_depth),
      adj_bottom = pmin(hzdepb_r, bottom_depth),
      thickness = adj_bottom - adj_top,
      weighted_clay = claytotal_r * thickness * comppct_r,
      weighted_sand = sandtotal_r * thickness * comppct_r,
      weighted_om = om_r * thickness * comppct_r,
      weighted_bd = dbthirdbar_r * thickness * comppct_r,
      depth_range = paste(top_depth, "-", bottom_depth, "cm")
    ) %>%
    dplyr::filter(thickness > 0)
  
  grouped <- ungrouped %>%
    dplyr::group_by(mukey, depth_range) %>%
    dplyr::summarize(
      clay_pct = sum(weighted_clay, na.rm = TRUE) / sum(thickness * comppct_r, na.rm = TRUE),
      sand_pct = sum(weighted_sand, na.rm = TRUE) / sum(thickness * comppct_r, na.rm = TRUE),
      silt_pct = 100 - clay_pct - sand_pct,
      om_pct = sum(weighted_om, na.rm = TRUE) / sum(thickness * comppct_r, na.rm = TRUE),
      bulk_density = sum(weighted_bd, na.rm = TRUE) / sum(thickness * comppct_r, na.rm = TRUE),
      .groups = 'drop'
    )
  grouped
}

# # --- Saxton & Rawls (2006) SSKS (cm/h) ---------------------------------------
.saxton_rawls_ssks <- function(theta_s, theta_33, theta_1500, coarse_fraction = 0, bulk_density = 1.4) {
  theta_s <- pmax(theta_s, 0.02)
  theta_33 <- pmax(theta_33, 0.01)
  theta_1500 <- pmax(theta_1500, 0.005)
  denom_lambda <- log(1500) - log(33)
  lambda <- ifelse(theta_33 > theta_1500 & theta_33 > 0 & theta_1500 > 0,
                   (log(theta_33) - log(theta_1500)) / denom_lambda,
                   0.5)
  Ks <- 1930 * (pmax(theta_s - theta_33, 1e-4))^(3 - lambda)
  denom <- 1 - coarse_fraction * (1 - 3 * (bulk_density / 2.65) / 2)
  denom <- ifelse(denom == 0, 1e-4, denom)
  Kb <- Ks * (1 - coarse_fraction) / denom / 10
  pmin(999.0, pmax(Kb, 0.001))
}

# --- DSSAT Formatting ---
format_dssat_soil_single <- function(profile_data, output_dir) {
  soil_id <- as.character(profile_data$ID[1])
  filename <- file.path(output_dir, paste0(soil_id, ".SOL"))
  
  # Redundant check for safety, though main loop handles this now
  if (file.exists(filename)) { return() }
  
  cat("*SOILS: USA SSURGO Soil Profiles\n", file = filename)
  cat("! Generated from SSURGO database\n\n", file = filename, append = TRUE)
  
  profile_depth <- max(as.numeric(sub(".*-", "", sub("cm", "", profile_data$depth_range))))
  cat(sprintf("*%-10s  %-11s %-5s %5.0f %s\n",
              soil_id, "SSURGO", "-99", profile_depth, "SSURGO profile"),
      file = filename, append = TRUE)
  
  cat("@SITE        COUNTRY          LAT     LONG SCS FAMILY\n", file = filename, append = TRUE)
  cat(sprintf(" %-11s %-11s %8.3f %8.3f \n",
              soil_id, "USA", profile_data$latitude[1], profile_data$longitude[1]),
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
                         ifelse(ssks_val >= 100, sprintf("%6.1f", pmin(999.0, ssks_val)),
                                sprintf("%6.2f", ssks_val)),
                         "   -99")
      
      cat(paste0(sprintf("%s   -99 %s %s %s  1.00%s %5.2f %5.2f %5.1f %5.1f   -99   -99   -99   -99   -99   -99\n",
                         depth_format, slll, sdul, ssat, ssks_str,
                         layer$bulk_density, layer$om_pct/1.724, layer$clay_pct, layer$silt_pct),
                 collapse = ""),
          file = filename, append = TRUE)
    })
  cat("\n", file = filename, append = TRUE)
}

# --- Main Processing Function for SSURGO (Optimized) ---
process_soils_ssurgo <- function(grid_points, output_dir_csv, output_dir_individual, n_cores,
                                 id_col, lat_col, long_col, format_sql_func) {
  
  message("Starting SSURGO Processing (Smart Resume Mode)...")

  if (missing(format_sql_func) || is.null(format_sql_func) || !is.function(format_sql_func)) {
    format_sql_func <- function(vec) {
      vals <- unique(stats::na.omit(as.character(vec)))
      if (length(vals) == 0) return("('')")
      paste0("('", paste(gsub("'", "''", vals), collapse = "','"), "')")
    }
  }

  dir.create(output_dir_individual, recursive = TRUE, showWarnings = FALSE)

  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores < 1) {
    stop("n_cores must be a positive integer.")
  }
  
  # --- 1. FILTER: Identify Missing Points BEFORE Loop ---
  # Convert sf to standard dataframe to avoid geometry overhead during filtering
  grid_df <- grid_points %>% 
    sf::st_drop_geometry() %>% 
    as.data.frame()
  
  # Ensure coords are present in dataframe
  if(!lat_col %in% names(grid_df)) grid_df[[lat_col]] <- sf::st_coordinates(grid_points)[,2]
  if(!long_col %in% names(grid_df)) grid_df[[long_col]] <- sf::st_coordinates(grid_points)[,1]
  
  # Get list of existing SOL files (Fast Directory Scan)
  if(dir.exists(output_dir_individual)) {
    existing_files <- tools::file_path_sans_ext(list.files(output_dir_individual, pattern = "\\.SOL$"))
    
    # Identify which points are NOT in existing files
    # Force character conversion to ensure accurate matching
    missing_mask <- ! (as.character(grid_df[[id_col]]) %in% existing_files)
    
    points_to_process <- grid_df[missing_mask, ]
    
    n_total <- nrow(grid_df)
    n_process <- nrow(points_to_process)
    n_skip <- n_total - n_process
    
    message(sprintf("Resume Check: Found %d existing profiles. Processing %d remaining points.", n_skip, n_process))
    
    if (n_process == 0) {
      message("All soil profiles already exist. Skipping SSURGO processing.")
      return(TRUE)
    }
    
  } else {
    # First run, process everything
    points_to_process <- grid_df
    message(sprintf("First Run: Processing all %d points.", nrow(points_to_process)))
  }
  
  # --- Worker Function ---
  # Returns the per-layer results data.frame on SUCCESS, or a tagged failure
  # record list(.fail=TRUE, ID, latitude, longitude, reason) on FAILURE. The
  # reason is prefixed with a category ("network:", "no-coverage:", "no-soil:",
  # "no-layers:") so the main process can tally and log exactly why each point
  # produced no .SOL -- instead of the old silent return(NULL).
  process_point_wrapper <- function(point_data_row) {
    ID   <- as.character(point_data_row[[id_col]])
    LATv <- point_data_row[[lat_col]]
    LONv <- point_data_row[[long_col]]
    fail <- function(reason) list(.fail = TRUE, ID = ID,
                                  latitude = LATv, longitude = LONv, reason = reason)

    # Reconstruct sf object for spatial query (WGS84; earlier steps reprojected)
    point_sf <- sf::st_as_sf(point_data_row, coords = c(long_col, lat_col), crs = 4326)

    # 1. Spatial query -> map unit key(s)
    sp <- robust_SDA_spatialQuery(point_sf, what = 'mukey')
    if (!isTRUE(sp$ok))
      return(fail(sprintf("network: SDA spatial query failed after retries (%s)",
                          if (is.na(sp$error)) "no detail" else sp$error)))
    soil_data_query <- sp$data
    if (is.null(soil_data_query) || !("mukey" %in% names(soil_data_query)) ||
        length(soil_data_query$mukey) == 0 || all(is.na(soil_data_query$mukey)))
      return(fail("no-coverage: no SSURGO map unit at this location (outside surveyed area / offshore)"))

    muname <- if ("muname" %in% names(soil_data_query))
                paste(unique(stats::na.omit(soil_data_query$muname)), collapse = "; ") else NA_character_

    # 2. Bedrock depth (optional; default 200 cm when unavailable)
    q_bedrock <- sprintf("SELECT mukey, brockdepmin FROM muaggatt WHERE mukey IN %s",
                         format_sql_func(soil_data_query$mukey))
    bq <- robust_SDA_query(q_bedrock)
    bedrock_depth <- 200
    if (isTRUE(bq$ok) && !is.null(bq$data)) {
      bd <- as.data.frame(bq$data)
      if (nrow(bd) > 0 && !all(is.na(bd$brockdepmin))) bedrock_depth <- min(bd$brockdepmin, na.rm = TRUE)
    }
    if (is.infinite(bedrock_depth)) bedrock_depth <- 200

    # Define Layers
    all_layers <- list("0-5cm"=c(0,5), "5-20cm"=c(5,20), "20-35cm"=c(20,35),
                       "35-50cm"=c(35,50), "50-65cm"=c(50,65), "65-80cm"=c(65,80),
                       "80-95cm"=c(80,95), "95-110cm"=c(95,110), "110-125cm"=c(110,125),
                       "125-140cm"=c(125,140), "140-155cm"=c(140,155), "155-170cm"=c(155,170),
                       "170-185cm"=c(170,185), "185-200cm"=c(185,200))

    valid_layers <- all_layers[sapply(all_layers, function(x) x[1] < bedrock_depth)]
    if(length(valid_layers) > 0) {
      last <- length(valid_layers)
      if(valid_layers[[last]][2] > bedrock_depth) valid_layers[[last]][2] <- bedrock_depth
      names(valid_layers)[last] <- paste0(valid_layers[[last]][1], "-", valid_layers[[last]][2], "cm")
    } else {
      valid_layers <- list()
      valid_layers[[paste0("0-", bedrock_depth, "cm")]] <- c(0, bedrock_depth)
    }

    # 3. Horizon properties
    q_soil <- sprintf("SELECT component.mukey, component.cokey, component.comppct_r,
                       chorizon.hzdept_r, chorizon.hzdepb_r, chorizon.claytotal_r,
                       chorizon.sandtotal_r, chorizon.om_r, chorizon.dbthirdbar_r
                       FROM component INNER JOIN chorizon ON component.cokey = chorizon.cokey
                       WHERE component.mukey IN %s", format_sql_func(soil_data_query$mukey))

    pq <- robust_SDA_query(q_soil)
    if (!isTRUE(pq$ok))
      return(fail(sprintf("network: soil-properties query failed after retries (%s)",
                          if (is.na(pq$error)) "no detail" else pq$error)))
    props <- as.data.frame(pq$data)
    if (nrow(props) == 0)
      return(fail(sprintf("no-soil: map unit has no soil horizons%s -- typically Water / Urban / Pits / Rock outcrop (mukey %s)",
                          if (!is.na(muname) && nzchar(muname)) sprintf(" [%s]", muname) else "",
                          paste(soil_data_query$mukey, collapse = ","))))

    # Calc Props per Layer
    results_list <- lapply(names(valid_layers), function(layer_name) {
      d <- valid_layers[[layer_name]]
      cp <- calculate_soil_properties(props, d[1], d[2])
      if(nrow(cp) > 0) cp$depth_range <- layer_name
      return(cp)
    })
    results_list <- results_list[sapply(results_list, function(x) !is.null(x) && nrow(x) > 0)]

    if (length(results_list) == 0)
      return(fail(sprintf("no-layers: horizon data present but no usable layers after depth filtering (bedrock %s cm; muname %s)",
                          bedrock_depth, if (!is.na(muname)) muname else "NA")))

    results_df <- do.call(rbind, results_list)
    # Add derived DSSAT physics
    results_df <- results_df %>% dplyr::mutate(
      ID=ID, longitude=point_data_row[[long_col]], latitude=point_data_row[[lat_col]],
      bedrock_depth_cm=bedrock_depth,
      sand_dec=sand_pct/100, clay_dec=clay_pct/100, om_dec=om_pct/100,
      theta_1500t = -0.024*sand_dec + 0.487*clay_dec + 0.006*om_dec + 0.005*(sand_dec*om_dec) - 0.013*(clay_dec*om_dec) + 0.068*(sand_dec*clay_dec) + 0.031,
      SLLL_raw = theta_1500t + (0.14*theta_1500t - 0.02),
      # Floor the wilting point: Saxton & Rawls can return SLLL <= 0 for very
      # sandy, low-clay layers, which is unphysical and makes DSSAT's soil-water
      # balance divide-by-zero (SIGFPE). Clamp to 0.02 cm3/cm3; the raw value is
      # kept (SLLL_raw) so the caller can tally and log the adjustment.
      SLLL = pmax(SLLL_raw, 0.02),
      theta_33t = -0.251*sand_dec + 0.195*clay_dec + 0.011*om_dec + 0.006*(sand_dec*om_dec) - 0.027*(clay_dec*om_dec) + 0.452*(sand_dec*clay_dec) + 0.299,
      SDUL_raw = theta_33t + (1.283*(theta_33t)^2 - 0.374*theta_33t - 0.015),
      # Floor field capacity so plant-available water (DUL-LL) stays usable: on
      # near-pure-sand layers Saxton & Rawls returns SDUL barely above the floored
      # SLLL (DUL-LL ~ 0.005-0.014 cm3/cm3), which drives DSSAT's soil-water
      # balance to a divide-by-(DUL-LL) singularity and SIGFPEs mid-season.
      # Guarantee DUL-LL >= 0.04; SDUL_raw is kept so the caller can tally/log it.
      SDUL = pmax(SDUL_raw, SLLL + 0.04),
      theta_s33t = 0.278*sand_dec + 0.034*clay_dec + 0.022*om_dec - 0.018*(sand_dec*om_dec) - 0.027*(clay_dec*om_dec) - 0.584*(sand_dec*clay_dec) + 0.078,
      theta_s33 = theta_s33t + (0.636*theta_s33t - 0.107),
      SSAT = SDUL + theta_s33 - 0.097*sand_dec + 0.043,
      SSKS = .saxton_rawls_ssks(SSAT, SDUL, SLLL, coarse_fraction = 0, bulk_density = bulk_density)
    )

    # WRITE .SOL FILE IMMEDIATELY
    format_dssat_soil_single(results_df, output_dir_individual)
    return(results_df)
  }
  
  # --- CHUNK SETUP ---
  total_points <- nrow(points_to_process)
  CHUNK_SIZE <- 10000 # Smaller chunk size for better resume granularity
  num_chunks <- ceiling(total_points / CHUNK_SIZE)
  
  message(sprintf("Processing batch of %d points in %d chunks...", total_points, num_chunks))
  
  # --- 1. SETUP CLUSTER (Done ONCE, outside the loop) ---
  if (.Platform$OS.type == "windows" && n_cores > 1) {
    message(sprintf("Initializing cluster with %d cores...", n_cores))
    cl <- parallel::makeCluster(n_cores)
    # Safety net: release the cluster even if processing errors out before the
    # explicit stopCluster() later. try() keeps it harmless on the normal path.
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, { library(soilDB); library(sf); library(dplyr); library(tidyr) })
    parallel::clusterExport(cl, varlist=c("process_point_wrapper", ".saxton_rawls_ssks", "robust_SDA_query", 
                                "robust_SDA_spatialQuery", "calculate_soil_properties",
                                "format_dssat_soil_single", "output_dir_individual", 
                                "id_col", "lat_col", "long_col", "format_sql_func"), 
                  envir = environment())
  } else if (n_cores > 1) {
    cl <- n_cores 
  } else {
    # Keep single-core runs in-process so local mocks and lightweight tests do
    # not fall through to live SDA requests in a worker process.
    cl <- NULL
  }
  
  # --- 2. CHUNKED LOOP ---
  all_fails <- list()  # accumulates per-point failure records across chunks
  for (i in 1:num_chunks) {
    start_idx <- (i - 1) * CHUNK_SIZE + 1
    end_idx <- min(i * CHUNK_SIZE, total_points)

    message(sprintf("  > Chunk %d/%d (Points %d - %d)", i, num_chunks, start_idx, end_idx))

    chunk_data <- points_to_process[start_idx:end_idx, ]
    chunk_list <- split(chunk_data, seq(nrow(chunk_data)))

    # Run Parallel on Chunk
    chunk_results <- pbapply::pblapply(chunk_list, process_point_wrapper, cl = cl)

    # Separate tagged failure records from successful per-layer data frames.
    is_fail <- function(x) is.list(x) && isTRUE(x[[".fail"]])
    all_fails     <- c(all_fails, Filter(is_fail, chunk_results))
    valid_results <- Filter(function(x) !is.null(x) && !is_fail(x), chunk_results)

    if(length(valid_results) > 0) {
      chunk_df <- dplyr::bind_rows(valid_results)

      # Surface wilting-point clamps (Saxton-Rawls artifact on sandy soils) in the
      # run log instead of applying them silently, then drop the helper column.
      if ("SLLL_raw" %in% names(chunk_df)) {
        clamped <- chunk_df[!is.na(chunk_df$SLLL_raw) & chunk_df$SLLL_raw < 0.02, ]
        if (nrow(clamped) > 0) {
          ids <- unique(as.character(clamped$ID))
          message(sprintf("[SSURGO] SLLL floored to 0.020 on %d layer(s) across %d point(s) (negative/low Saxton-Rawls wilting point; sandy soils): %s",
                          nrow(clamped), length(ids), paste(ids, collapse = ", ")))
        }
        chunk_df$SLLL_raw <- NULL
      }

      # Surface field-capacity raises (minimum available-water enforcement) too.
      if ("SDUL_raw" %in% names(chunk_df)) {
        raised <- chunk_df[!is.na(chunk_df$SDUL_raw) & (chunk_df$SDUL - chunk_df$SDUL_raw) > 1e-9, ]
        if (nrow(raised) > 0) {
          ids <- unique(as.character(raised$ID))
          message(sprintf("[SSURGO] SDUL raised to keep DUL-LL >= 0.040 on %d layer(s) across %d point(s) (near-zero plant-available water; sandy soils): %s",
                          nrow(raised), length(ids), paste(ids, collapse = ", ")))
        }
        chunk_df$SDUL_raw <- NULL
      }

      # Write to CSV (Append Mode)
      readr::write_csv(chunk_df, output_dir_csv, append = file.exists(output_dir_csv))
    }

    # --- CRITICAL: CLEAR MEMORY ---
    rm(chunk_list, chunk_results, valid_results)
    if(exists("chunk_df")) rm(chunk_df)
    gc()
  }

  # --- 3. CLEANUP ---
  if (.Platform$OS.type == "windows" && !is.null(cl)) {
    parallel::stopCluster(cl)
  }

  # --- 4. FAILURE REPORT ---
  # Explain, per point, IF and WHY no soil profile was produced -- written next to
  # the soil CSV so missing .SOL files are auditable instead of silent.
  if (length(all_fails) > 0) {
    fail_df <- data.frame(
      ID        = vapply(all_fails, function(f) as.character(f$ID), character(1)),
      latitude  = vapply(all_fails, function(f) suppressWarnings(as.numeric(f$latitude)),  numeric(1)),
      longitude = vapply(all_fails, function(f) suppressWarnings(as.numeric(f$longitude)), numeric(1)),
      reason    = vapply(all_fails, function(f) as.character(f$reason), character(1)),
      stringsAsFactors = FALSE
    )
    failure_log <- file.path(dirname(output_dir_csv),
        paste0(tools::file_path_sans_ext(basename(output_dir_csv)), "_download_failures.csv"))
    readr::write_csv(fail_df, failure_log)

    message(sprintf("[SSURGO] %d of %d processed point(s) produced NO soil profile:",
                    nrow(fail_df), total_points))
    cats <- sub(":.*$", "", fail_df$reason)
    tb <- sort(table(cats), decreasing = TRUE)
    for (nm in names(tb)) message(sprintf("   - %-13s %d   (%s)", paste0(nm, ":"), tb[[nm]],
        switch(nm,
          "network"     = "transient SDA/server/timeout -- re-run to retry these",
          "no-coverage" = "no SSURGO map unit here (outside survey area / offshore)",
          "no-soil"     = "non-soil map unit (Water, Urban, Pits, Rock) -- no horizons exist",
          "no-layers"   = "horizons present but unusable after depth filtering",
          "other")))
    message(sprintf("   Per-point details (ID, lat, long, reason) -> %s", failure_log))
  } else {
    message("[SSURGO] All processed points produced a soil profile.")
  }

  message("SSURGO Processing Complete.")
  return(TRUE)
}
