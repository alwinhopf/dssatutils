# ==============================================================================
#  SOIL HELPER: SSURGO (USDA WEB SERVICE) - DSSAT profile generation (Alderman method)
#  Filename: soil_ssurgo_alderman.R
# ==============================================================================

# ---- Retry wrappers ----------------------------------------------------------
# NOTE: these are deliberately suffixed "_alderman" and return the RAW SDA result
# (data.frame, or NULL on failure). They must NOT be named robust_SDA_query /
# robust_SDA_spatialQuery: soil_ssurgo.R defines those names with a DIFFERENT
# contract (a list(ok=, data=, error=)), and because the package has no Collate:
# field, R sources files alphabetically — so an identically named definition here
# would silently override soil_ssurgo.R's and break SSURGO/gNATSGO (which expect
# the list form). Keep these names unique.
robust_SDA_query_alderman <- function(query, max_retries = 3, retry_delay_seconds = 5) {
  for (attempt in seq_len(max_retries)) {
    result <- try(soilDB::SDA_query(query), silent = TRUE)
    if (!inherits(result, "try-error")) return(result)
    Sys.sleep(retry_delay_seconds)
  }
  NULL
}

robust_SDA_spatialQuery_alderman <- function(point_sf, what, max_retries = 3, retry_delay_seconds = 5) {
  for (attempt in seq_len(max_retries)) {
    result <- try(soilDB::SDA_spatialQuery(point_sf, what = what), silent = TRUE)
    if (!inherits(result, "try-error")) return(result)
    Sys.sleep(retry_delay_seconds)
  }
  NULL
}

# ---- Utility helpers ---------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

safe_first_non_na <- function(x, default = NA) {
  x <- x[!is.na(x)]
  if (length(x) == 0) default else x[1]
}

coalesce_num <- function(...) {
  vals <- list(...)
  out <- vals[[1]]
  if (length(vals) == 1) return(out)
  for (i in 2:length(vals)) {
    replace_idx <- is.na(out)
    out[replace_idx] <- vals[[i]][replace_idx]
  }
  out
}

clip01 <- function(x) pmin(pmax(x, 0.001), 0.95)

sanitize_char <- function(x, default = "-99") {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- default
  x
}

append_log_line <- function(log_file, level = "INFO", context = "SSURGO", msg = "", id = NULL,
                            lock_timeout_seconds = 30, lock_retry_seconds = 0.05) {
  if (is.null(log_file) || length(log_file) == 0 || is.na(log_file)[1] || !nzchar(as.character(log_file)[1])) {
    return(invisible(NULL))
  }

  log_file <- as.character(log_file)[1]
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  lock_dir <- paste0(log_file, ".lockdir")
  start_time <- Sys.time()
  lock_acquired <- FALSE

  while (!lock_acquired) {
    lock_acquired <- isTRUE(suppressWarnings(dir.create(lock_dir, showWarnings = FALSE, recursive = FALSE)))
    if (lock_acquired) break

    waited <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (is.finite(waited) && waited >= lock_timeout_seconds) {
      stop(sprintf("Timed out waiting for log lock: %s", lock_dir))
    }
    Sys.sleep(lock_retry_seconds)
  }

  on.exit({
    if (dir.exists(lock_dir)) unlink(lock_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  clean_msg <- gsub("[\r\n\t]+", " ", paste(as.character(msg), collapse = " "))
  id_part <- if (!is.null(id) && length(id) > 0 && !is.na(id)[1] && nzchar(as.character(id)[1])) {
    sprintf(" [ID=%s]", as.character(id)[1])
  } else {
    ""
  }
  pid_part <- sprintf(" [PID=%s]", Sys.getpid())
  line <- sprintf("[%s] [%s] [%s]%s%s %s",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                  level,
                  context,
                  id_part,
                  pid_part,
                  clean_msg)
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(NULL)
}

soil_helper_log <- function(log_file = NULL, level = "INFO", context = "SSURGO", msg = "", point_id = NULL) {
  try(append_log_line(log_file, level, context, msg, id = point_id), silent = TRUE)
  invisible(NULL)
}

# ---- Pedotransfer functions (adapted from csmsoil repo) ----------------------
soil_ptf_saxton_slll <- function(silt, clay, soc, bulk_density, coarse_fraction) {
  sand <- 1 - silt / 100 - clay / 100
  clay <- clay / 100
  som <- soc * 1.72
  theta_1500t <- -0.024 * sand + 0.487 * clay + 0.006 * som + 0.005 * sand * som -
    0.013 * clay * som + 0.068 * sand * clay + 0.031
  theta_1500 <- theta_1500t + 0.14 * theta_1500t - 0.02
  a <- bulk_density / 2.65
  Rv <- (a * coarse_fraction) / (1 - coarse_fraction * (1 - a))
  theta_1500 * (1 - Rv)
}

soil_ptf_saxton_sdul <- function(silt, clay, soc, bulk_density, coarse_fraction) {
  sand <- 1 - silt / 100 - clay / 100
  clay <- clay / 100
  som <- soc * 1.72
  theta_33t <- -0.251 * sand + 0.195 * clay + 0.011 * som + 0.006 * sand * som -
    0.013 * clay * som + 0.452 * sand * clay + 0.299
  theta_33 <- theta_33t + 1.283 * (theta_33t^2) - 0.374 * theta_33t - 0.015
  a <- bulk_density / 2.65
  Rv <- (a * coarse_fraction) / (1 - coarse_fraction * (1 - a))
  theta_33 * (1 - Rv)
}

soil_ptf_saxton_ssat <- function(silt, clay, soc, bulk_density, coarse_fraction) {
  sand <- 1 - silt / 100 - clay / 100
  clay <- clay / 100
  som <- soc * 1.72
  theta_33t <- -0.251 * sand + 0.195 * clay + 0.011 * som + 0.006 * sand * som -
    0.013 * clay * som + 0.452 * sand * clay + 0.299
  theta_33 <- theta_33t + 1.283 * (theta_33t^2) - 0.374 * theta_33t - 0.015
  theta_S33t <- 0.278 * sand + 0.034 * clay + 0.022 * som - 0.018 * sand * som -
    0.027 * clay * som - 0.584 * sand * clay + 0.078
  theta_S33 <- theta_S33t + 0.636 * theta_S33t - 0.107
  theta_S <- theta_33 + theta_S33 - 0.097 * sand + 0.043
  a <- bulk_density / 2.65
  Rv <- (a * coarse_fraction) / (1 - coarse_fraction * (1 - a))
  theta_S * (1 - Rv)
}

soil_ptf_saxton_ssks <- function(theta_s, theta_33, theta_1500, coarse_fraction, bulk_density) {
  theta_s <- pmax(theta_s, 0.02)
  theta_33 <- pmax(theta_33, 0.01)
  theta_1500 <- pmax(theta_1500, 0.005)
  lambda <- (log(theta_33) - log(theta_1500)) / (log(1500) - log(33))
  Ks <- 1930 * (pmax(theta_s - theta_33, 1e-4))^(3 - lambda)
  Kb <- Ks * (1 - coarse_fraction) / (1 - coarse_fraction * (1 - 3 * (bulk_density / 2.65) / 2)) / 10
  pmax(Kb, 0.001)
}

soil_ptf_slu1 <- function(sat, pwp, depth) {
  if (any(depth > 15)) {
    lyr_gt15 <- which(depth > 15)[1]
    sat <- sat[1:lyr_gt15]
    pwp <- pwp[1:lyr_gt15]
    depth <- depth[1:lyr_gt15]
    depth[lyr_gt15] <- 15
  }
  sum((sat - pwp) / 2 * diff(c(0, depth)) * 10, na.rm = TRUE)
}

soil_ptf_nrcs_hsg <- function(ksat, depth) {
  if (any(depth > 50)) {
    gt50 <- which(depth > 50)[1]
    ksat_50 <- min(ksat[1:gt50], na.rm = TRUE)
  } else {
    ksat_50 <- min(ksat, na.rm = TRUE)
  }
  if (any(depth > 100)) {
    gt100 <- which(depth > 100)[1]
    ksat_100 <- min(ksat[1:gt100], na.rm = TRUE)
  } else {
    ksat_100 <- min(ksat, na.rm = TRUE)
  }
  sl_depth <- max(depth, na.rm = TRUE)
  hsg <- NA_real_
  if (!is.finite(ksat_50) || !is.finite(ksat_100) || !is.finite(sl_depth)) return(hsg)
  if ((ksat_50 > 40 * 0.36 && sl_depth >= 50 && sl_depth <= 100) || (ksat_100 > 10 * 0.36 && sl_depth > 100)) hsg <- 1
  else if ((ksat_50 <= 40 * 0.36 && ksat_50 > 10 * 0.36 && sl_depth >= 50 && sl_depth <= 100) || (ksat_100 <= 10 * 0.36 && ksat_100 > 4 * 0.36 && sl_depth > 100)) hsg <- 2
  else if ((ksat_50 <= 10 * 0.36 && ksat_50 > 1 * 0.36 && sl_depth >= 50 && sl_depth <= 100) || (ksat_100 <= 4 * 0.36 && ksat_100 > 0.4 * 0.36 && sl_depth > 100)) hsg <- 3
  else if ((ksat_50 <= 1 * 0.36 && sl_depth >= 50 && sl_depth <= 100) || sl_depth < 50 || (ksat_100 <= 0.4 * 0.36 && sl_depth > 100)) hsg <- 4
  hsg
}

soil_ptf_curve_number <- function(slope, hsg, ksat = NULL, depth = NULL) {
  if (missing(hsg) || all(is.na(hsg))) {
    hsg <- soil_ptf_nrcs_hsg(ksat, depth)
  }
  curve_number <- hsg
  curve_number[slope >= 0 & slope <= 2 & hsg == 1] <- 61
  curve_number[slope > 2 & slope <= 5 & hsg == 1] <- 64
  curve_number[slope > 5 & slope <= 10 & hsg == 1] <- 68
  curve_number[slope > 10 & hsg == 1] <- 71
  curve_number[slope > 0 & slope <= 2 & hsg == 2] <- 73
  curve_number[slope > 2 & slope <= 5 & hsg == 2] <- 76
  curve_number[slope > 5 & slope <= 10 & hsg == 2] <- 80
  curve_number[slope > 10 & hsg == 2] <- 83
  curve_number[slope >= 0 & slope <= 2 & hsg == 3] <- 81
  curve_number[slope > 2 & slope <= 5 & hsg == 3] <- 84
  curve_number[slope > 5 & slope <= 10 & hsg == 3] <- 88
  curve_number[slope > 10 & hsg == 3] <- 91
  curve_number[slope >= 0 & slope <= 2 & hsg == 4] <- 84
  curve_number[slope > 2 & slope <= 5 & hsg == 4] <- 87
  curve_number[slope > 5 & slope <= 10 & hsg == 4] <- 91
  curve_number[slope > 10 & hsg == 4] <- 94
  curve_number[is.na(curve_number) | !is.finite(curve_number)] <- 73
  curve_number
}

sldr_from_drainage <- function(drainage) {
  switch(as.character(drainage),
         "Excessively drained" = 0.85,
         "Somewhat excessively drained" = 0.75,
         "Well drained" = 0.60,
         "Moderately well drained" = 0.40,
         "Somewhat poorly drained" = 0.25,
         "Poorly drained" = 0.05,
         "Very poorly drained" = 0.01,
         0.60)
}

# ---- SSURGO query helpers ----------------------------------------------------
sql_in_from_values <- function(x) {
  vals <- unique(na.omit(as.character(x)))
  if (length(vals) == 0) return("('')")
  escaped <- gsub("'", "''", vals, fixed = TRUE)
  paste0("(", paste(sprintf("'%s'", escaped), collapse = ","), ")")
}

get_point_component_table <- function(point_sf, point_id = NULL, log_file = NULL) {
  spatial_hit <- robust_SDA_spatialQuery_alderman(point_sf, what = 'mukey')
  mukeys <- character(0)
  if (!is.null(spatial_hit) && 'mukey' %in% names(spatial_hit)) {
    mukeys <- unique(na.omit(as.character(spatial_hit$mukey)))
  }

  if (length(mukeys) > 0) {
    soil_helper_log(log_file, "INFO", "SSURGO_COMPONENTS",
                    sprintf("Spatial query returned %d mukey(s): %s",
                            length(mukeys), paste(head(mukeys, 5), collapse = ", ")),
                    point_id = point_id)
  }

  if (length(mukeys) == 0) {
    soil_helper_log(log_file, "WARN", "SSURGO_COMPONENTS",
                    "Spatial query returned no mukey; trying WKT intersection fallback",
                    point_id = point_id)
    wkt <- sf::st_as_text(sf::st_geometry(point_sf))
    mukey_query <- paste0(
      "SELECT DISTINCT mu.mukey AS mukey ",
      "FROM mapunit mu ",
      "WHERE mu.mukey IN (SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('", wkt, "'))"
    )
    mukey_out <- robust_SDA_query_alderman(mukey_query)
    if (!is.null(mukey_out) && nrow(mukey_out) > 0 && 'mukey' %in% names(mukey_out)) {
      mukeys <- unique(na.omit(as.character(mukey_out$mukey)))
      soil_helper_log(log_file, "INFO", "SSURGO_COMPONENTS",
                      sprintf("WKT fallback returned %d mukey(s): %s",
                              length(mukeys), paste(head(mukeys, 5), collapse = ", ")),
                      point_id = point_id)
    }
  }

  if (length(mukeys) == 0) {
    soil_helper_log(log_file, "ERROR", "SSURGO_COMPONENTS",
                    "No mukey could be resolved for this point",
                    point_id = point_id)
    return(NULL)
  }

  comp_query <- paste0(
    "SELECT compname, cokey, mukey, COALESCE(comppct_r,'') AS comppct_r, ",
    "COALESCE(hydgrp,'') AS hydgrp, COALESCE(slope_r,'') AS slope_r, ",
    "COALESCE(drainagecl,'') AS drainage, COALESCE(albedodry_r,'') AS albedodry_r ",
    "FROM component WHERE mukey IN ", sql_in_from_values(mukeys)
  )
  comp_tbl <- robust_SDA_query_alderman(comp_query)
  if (is.null(comp_tbl)) {
    soil_helper_log(log_file, "ERROR", "SSURGO_COMPONENTS",
                    "Component query returned NULL",
                    point_id = point_id)
    return(NULL)
  }
  if (nrow(comp_tbl) == 0) {
    soil_helper_log(log_file, "ERROR", "SSURGO_COMPONENTS",
                    sprintf("Component query returned 0 rows for mukey(s): %s", paste(mukeys, collapse = ";")),
                    point_id = point_id)
    return(NULL)
  }
  comp_tbl <- as.data.frame(comp_tbl)
  soil_helper_log(log_file, "INFO", "SSURGO_COMPONENTS",
                  sprintf("Component query returned %d row(s) across %d cokey(s)",
                          nrow(comp_tbl), dplyr::n_distinct(comp_tbl$cokey)),
                  point_id = point_id)
  comp_tbl
}

get_component_horizons <- function(cokey, point_id = NULL, log_file = NULL) {
  q <- paste0(
    "SELECT chorizon.hzdept_r, chorizon.hzdepb_r, chorizon.dbovendry_r, ",
    "chorizon.dbtenthbar_r, chorizon.dbthirdbar_r, chorizon.dbfifteenbar_r, ",
    "chorizon.wsatiated_r, chorizon.wtenthbar_r, chorizon.wthirdbar_r, ",
    "chorizon.partdensity, chorizon.ksat_r, chorizon.wfifteenbar_r, ",
    "chorizon.sandtotal_r, chorizon.claytotal_r, chorizon.silttotal_r, ",
    "chorizon.om_r, chorizon.hzname, chfrags.fragvol_r AS fragvol_r, ",
    "chorizon.cokey FROM chorizon LEFT JOIN chfrags ON chfrags.chkey = chorizon.chkey ",
    "WHERE chorizon.cokey IN ('", cokey, "') ORDER BY chorizon.hzdepb_r"
  )
  hz <- robust_SDA_query_alderman(q)
  if (is.null(hz)) {
    soil_helper_log(log_file, "ERROR", "SSURGO_HORIZONS",
                    sprintf("Horizon query returned NULL for cokey=%s", as.character(cokey)[1]),
                    point_id = point_id)
    return(NULL)
  }
  if (nrow(hz) == 0) {
    soil_helper_log(log_file, "ERROR", "SSURGO_HORIZONS",
                    sprintf("No horizons returned for cokey=%s", as.character(cokey)[1]),
                    point_id = point_id)
    return(NULL)
  }
  hz <- as.data.frame(hz)
  soil_helper_log(log_file, "INFO", "SSURGO_HORIZONS",
                  sprintf("Horizon query returned %d row(s) for cokey=%s", nrow(hz), as.character(cokey)[1]),
                  point_id = point_id)
  hz
}

choose_dominant_component <- function(comp_tbl) {
  comp_tbl <- comp_tbl %>%
    dplyr::mutate(comppct_r = suppressWarnings(as.numeric(comppct_r)),
                  slope_r = suppressWarnings(as.numeric(slope_r)),
                  albedodry_r = suppressWarnings(as.numeric(albedodry_r))) %>%
    dplyr::arrange(dplyr::desc(comppct_r), cokey)
  comp_tbl[1, , drop = FALSE]
}

# ---- Horizon-to-DSSAT conversion ---------------------------------------------
calc_sbdm <- function(h_tbl) {
  sbdm <- coalesce_num(h_tbl$dbtenthbar_r, h_tbl$dbthirdbar_r, h_tbl$dbovendry_r)
  as.numeric(sbdm)
}

is_coarse_texture <- function(sand, silt, clay) {
  (sand >= 85 & (silt + 1.5 * clay) <= 15) |
    ((sand >= 85 & sand < 90 & (silt + 1.5 * clay) >= 15) |
       (sand >= 70 & sand < 85 & (silt + 2 * clay) <= 30)) |
    ((clay <= 20 & sand >= 52 & (silt + 2 * clay) > 30) |
       (clay < 7 & silt < 50 & sand > 43 & sand < 52))
}

calc_sdul_measured <- function(h_tbl) {
  coarse <- is_coarse_texture(h_tbl$sandtotal_r, h_tbl$silttotal_r, h_tbl$claytotal_r)
  wthird <- h_tbl$wthirdbar_r / 100
  wtenth <- coalesce_num(h_tbl$wtenthbar_r / 100, wthird)
  ifelse(coarse, wtenth, wthird)
}

calc_ssat_measured <- function(h_tbl) {
  partdensity <- coalesce_num(h_tbl$partdensity, rep(2.65, nrow(h_tbl)))
  dbtenthbar_r <- ifelse(is.na(h_tbl$wtenthbar_r), NA, h_tbl$dbtenthbar_r)
  dbthirdbar_r <- ifelse(is.na(h_tbl$wthirdbar_r), NA, h_tbl$dbthirdbar_r)
  ssat_tenth <- 0.95 * (1 - dbtenthbar_r / partdensity)
  ssat_third <- 0.95 * (1 - dbthirdbar_r / partdensity)
  ssat_wsat <- h_tbl$wsatiated_r / 100
  coalesce_num(ssat_wsat, ssat_tenth, ssat_third, 0.95 * (1 - h_tbl$dbovendry_r / partdensity))
}

build_dssat_profile_from_component <- function(component_row, horizon_tbl, point_id, lat, lon, log_file = NULL) {
  soil_helper_log(log_file, "INFO", "SSURGO_DOMINANT",
                  sprintf("Building dominant profile from cokey=%s with %d raw horizon row(s)",
                          as.character(component_row$cokey[1]), nrow(horizon_tbl)),
                  point_id = point_id)
  hz <- horizon_tbl %>%
    dplyr::mutate(dplyr::across(c(hzdept_r, hzdepb_r, dbovendry_r, dbtenthbar_r, dbthirdbar_r,
                    dbfifteenbar_r, wsatiated_r, wtenthbar_r, wthirdbar_r,
                    partdensity, ksat_r, wfifteenbar_r, sandtotal_r, claytotal_r,
                    silttotal_r, om_r, fragvol_r), ~ suppressWarnings(as.numeric(.)))) %>%
    dplyr::arrange(hzdepb_r) %>%
    dplyr::group_by(hzdepb_r, hzname, cokey) %>%
    dplyr::summarize(
      hzdept_r = suppressWarnings(min(hzdept_r, na.rm = TRUE)),
      dbovendry_r = mean(dbovendry_r, na.rm = TRUE),
      dbtenthbar_r = mean(dbtenthbar_r, na.rm = TRUE),
      dbthirdbar_r = mean(dbthirdbar_r, na.rm = TRUE),
      dbfifteenbar_r = mean(dbfifteenbar_r, na.rm = TRUE),
      wsatiated_r = mean(wsatiated_r, na.rm = TRUE),
      wtenthbar_r = mean(wtenthbar_r, na.rm = TRUE),
      wthirdbar_r = mean(wthirdbar_r, na.rm = TRUE),
      partdensity = mean(partdensity, na.rm = TRUE),
      ksat_r = mean(ksat_r, na.rm = TRUE),
      wfifteenbar_r = mean(wfifteenbar_r, na.rm = TRUE),
      sandtotal_r = mean(sandtotal_r, na.rm = TRUE),
      claytotal_r = mean(claytotal_r, na.rm = TRUE),
      silttotal_r = mean(silttotal_r, na.rm = TRUE),
      om_r = mean(om_r, na.rm = TRUE),
      fragvol_r = if (all(is.na(fragvol_r))) NA_real_ else sum(fragvol_r, na.rm = TRUE),
      .groups = "drop"
    )

  # replace non-finite summaries from all-NA groups
  for (nm in names(hz)) {
    if (is.numeric(hz[[nm]])) hz[[nm]][!is.finite(hz[[nm]])] <- NA_real_
  }

  hz <- hz %>%
    dplyr::mutate(
      bedrock = stringr::str_detect(dplyr::coalesce(hzname, ""), "[rR]"),
      fragvol_raw = ifelse(is.na(fragvol_r), NA_real_, pmin(fragvol_r, 99)),
      fragvol_r = coalesce_num(fragvol_raw, 0),
      coarse_fraction = fragvol_r / 100,
      partdensity = coalesce_num(partdensity, 2.65),
      SBDM = calc_sbdm(.),
      SSAT_raw = calc_ssat_measured(.),
      SDUL_raw = calc_sdul_measured(.),
      SLLL_raw = wfifteenbar_r / 100,
      soc = om_r / 1.724,
      SLLL_ptf = soil_ptf_saxton_slll(
        silttotal_r, claytotal_r, soc,
        coalesce_num(SBDM, dbthirdbar_r, dbovendry_r, 1.4),
        coarse_fraction
      ),
      SDUL_ptf = soil_ptf_saxton_sdul(
        silttotal_r, claytotal_r, soc,
        coalesce_num(SBDM, dbthirdbar_r, dbovendry_r, 1.4),
        coarse_fraction
      ),
      SSAT_ptf = soil_ptf_saxton_ssat(
        silttotal_r, claytotal_r, soc,
        coalesce_num(SBDM, dbthirdbar_r, dbovendry_r, 1.4),
        coarse_fraction
      ),
      SLLL = clip01(coalesce_num(SLLL_raw, ifelse(bedrock, NA_real_, SLLL_ptf))),
      SDUL = clip01(coalesce_num(SDUL_raw, ifelse(bedrock, NA_real_, SDUL_ptf))),
      SSAT = clip01(coalesce_num(SSAT_raw, ifelse(bedrock, NA_real_, SSAT_ptf))),
      SDUL = pmax(SDUL, SLLL + 0.005),
      SSAT = pmax(SSAT, SDUL + 0.01),
      SLCF = ifelse(bedrock & is.na(fragvol_raw), 99, coalesce_num(fragvol_raw, 0)),
      SRGF = ifelse(bedrock, pmax(0.01, 1 - SLCF / 100), 1.0),
      SSKS = ifelse(
        (bedrock & is.na(ksat_r)) | (!is.na(ksat_r) & ksat_r < 0.001 / 60 / 60 * 10000),
        0.001,
        coalesce_num(
          ksat_r * 60 * 60 / 10000,
          soil_ptf_saxton_ssks(SSAT, SDUL, SLLL, coarse_fraction, coalesce_num(SBDM, 1.4))
        )
      ),
      SLOC = coalesce_num(soc, 0),
      SLCL = coalesce_num(claytotal_r, pmax(0, 100 - coalesce_num(sandtotal_r, 0) - coalesce_num(silttotal_r, 0))),
      SLSI = coalesce_num(silttotal_r, pmax(0, 100 - coalesce_num(sandtotal_r, 0) - coalesce_num(claytotal_r, 0))),
      SLMH = sanitize_char(hzname, "-99"),
      SLNI = NA_real_,
      SLHW = NA_real_,
      SLHB = NA_real_,
      SCEC = NA_real_,
      SADC = NA_real_
    ) %>%
    tidyr::fill(SDUL, SLLL, SSAT, SLSI, SLCL, SBDM, SLOC, .direction = "down") %>%
    dplyr::mutate(
      SLB = as.integer(round(hzdepb_r)),
      SBDM = coalesce_num(SBDM, dbthirdbar_r, dbovendry_r, 1.4),
      SALB = coalesce_num(as.numeric(component_row$albedodry_r), 0.13),
      SLDR = sldr_from_drainage(component_row$drainage),
      slope_r = suppressWarnings(as.numeric(component_row$slope_r)),
      hydgrp = as.character(component_row$hydgrp),
      SLRO = soil_ptf_curve_number(safe_first_non_na(slope_r, 0), ifelse(substr(hydgrp, 1, 1) %in% c("A","B","C","D"), match(substr(hydgrp, 1, 1), c("A","B","C","D")), NA_real_), SSKS, SLB)
    ) %>%
    dplyr::select(SLB, SLMH, SLLL, SDUL, SSAT, SRGF, SSKS, SBDM, SLOC,
           SLCL, SLSI, SLCF, SLNI, SLHW, SLHB, SCEC, SADC) %>%
    dplyr::filter(!is.na(SLB), SLB > 0) %>%
    dplyr::distinct(SLB, .keep_all = TRUE) %>%
    dplyr::arrange(SLB)

  if (nrow(hz) == 0) {
    soil_helper_log(log_file, "ERROR", "SSURGO_DOMINANT",
                    "Dominant component profile had 0 valid DSSAT layers after filtering",
                    point_id = point_id)
    return(NULL)
  }

  soil_helper_log(log_file, "INFO", "SSURGO_DOMINANT",
                  sprintf("Dominant profile retained %d DSSAT layer(s): %s",
                          nrow(hz), paste(hz$SLB, collapse = ",")),
                  point_id = point_id)

  slu1 <- soil_ptf_slu1(hz$SSAT, hz$SLLL, hz$SLB)
  if (!is.finite(slu1) || is.na(slu1)) slu1 <- 6

  list(
    profile_id = as.character(point_id),
    latitude = as.numeric(lat),
    longitude = as.numeric(lon),
    site = as.character(point_id),
    country = "USA",
    scs_family = "",
    scom = "SC",
    salb = coalesce_num(as.numeric(component_row$albedodry_r), 0.13),
    slu1 = round(slu1, 1),
    sldr = sldr_from_drainage(component_row$drainage),
    slro = soil_ptf_curve_number(safe_first_non_na(suppressWarnings(as.numeric(component_row$slope_r)), 0),
                                 ifelse(substr(as.character(component_row$hydgrp), 1, 1) %in% c("A","B","C","D"),
                                        match(substr(as.character(component_row$hydgrp), 1, 1), c("A","B","C","D")),
                                        soil_ptf_nrcs_hsg(hz$SSKS, hz$SLB)),
                                 hz$SSKS, hz$SLB),
    slnf = 1,
    slpf = 1,
    smhb = "IB001",
    smpx = "IB001",
    smke = "IB001",
    layers = hz,
    metadata = data.frame(
      ID = as.character(point_id),
      SOIL_ID = as.character(point_id),
      mukey = as.character(component_row$mukey),
      cokey = as.character(component_row$cokey),
      compname = as.character(component_row$compname),
      comppct_r = suppressWarnings(as.numeric(component_row$comppct_r)),
      latitude = as.numeric(lat),
      longitude = as.numeric(lon),
      stringsAsFactors = FALSE
    )
  )
}

# ---- DSSAT writing -----------------------------------------------------------
format_dssat_decimal <- function(x, digits = 3, width = 5) {
  if (is.na(x)) return(sprintf(paste0("%", width, "s"), "-99"))
  out <- sprintf(paste0("%", width, ".", digits, "f"), x)
  sub("^0", " ", out)
}

format_dssat_numeric <- function(x, width = 5, digits = 1) {
  if (is.na(x)) return(sprintf(paste0("%", width, "s"), "-99"))
  sprintf(paste0("%", width, ".", digits, "f"), x)
}

write_dssat_soil_file <- function(profile, output_dir) {
  filename <- file.path(output_dir, paste0(profile$profile_id, ".SOL"))
  con <- file(filename, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("*SOILS: USA SSURGO Soil Profiles", con)
  writeLines("! Generated from SSURGO database using full-profile logic", con)
  writeLines("", con)
  writeLines(sprintf("*%-10s SSURGO        %9.3f %9.3f",
                     substr(profile$profile_id, 1, 10), profile$latitude, profile$longitude), con)
  writeLines("@SITE        COUNTRY          LAT     LONG SCS FAMILY", con)
  writeLines(sprintf(" %-11s %-10s %9.3f %9.3f %s",
                     substr(profile$site, 1, 11), profile$country,
                     profile$latitude, profile$longitude,
                     substr(profile$scs_family, 1, 20)), con)
  writeLines("@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE", con)
  writeLines(sprintf(" %5s %5.2f %5.1f %5.2f %5.0f %5.0f %5.0f %5s %5s %5s",
                     substr(profile$scom, 1, 5), profile$salb, profile$slu1,
                     profile$sldr, profile$slro, profile$slnf, profile$slpf,
                     profile$smhb, profile$smpx, profile$smke), con)
  writeLines("@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC", con)

  for (i in seq_len(nrow(profile$layers))) {
    lyr <- profile$layers[i, ]
    line <- sprintf(
      "%5d %5s %5s %5s %5s %5.2f %5.2f %5.2f %5.2f %5.1f %5.1f %5.0f %5s %5s %5s %5s %5s",
      as.integer(lyr$SLB),
      substr(sanitize_char(lyr$SLMH, "-99"), 1, 5),
      format_dssat_decimal(lyr$SLLL, 3, 5),
      format_dssat_decimal(lyr$SDUL, 3, 5),
      format_dssat_decimal(lyr$SSAT, 3, 5),
      coalesce_num(lyr$SRGF, 1),
      coalesce_num(lyr$SSKS, -99),
      coalesce_num(lyr$SBDM, -99),
      coalesce_num(lyr$SLOC, -99),
      coalesce_num(lyr$SLCL, -99),
      coalesce_num(lyr$SLSI, -99),
      coalesce_num(lyr$SLCF, -99),
      ifelse(is.na(lyr$SLNI), "  -99", format_dssat_decimal(lyr$SLNI, 3, 5)),
      ifelse(is.na(lyr$SLHW), "  -99", format_dssat_decimal(lyr$SLHW, 3, 5)),
      ifelse(is.na(lyr$SLHB), "  -99", format_dssat_decimal(lyr$SLHB, 3, 5)),
      ifelse(is.na(lyr$SCEC), "  -99", format_dssat_decimal(lyr$SCEC, 3, 5)),
      ifelse(is.na(lyr$SADC), "  -99", format_dssat_decimal(lyr$SADC, 3, 5))
    )
    writeLines(line, con)
  }
  writeLines("", con)
}

calculate_soil_properties_fallback <- function(soil_properties, top_depth, bottom_depth) {
  ungrouped <- soil_properties %>%
    dplyr::mutate(
      adj_top = pmax(suppressWarnings(as.numeric(hzdept_r)), top_depth),
      adj_top = ifelse(is.na(adj_top), top_depth, adj_top),
      adj_bottom = pmin(suppressWarnings(as.numeric(hzdepb_r)), bottom_depth),
      adj_bottom = ifelse(is.na(adj_bottom), bottom_depth, adj_bottom),
      thickness = adj_bottom - adj_top,
      weighted_clay = suppressWarnings(as.numeric(claytotal_r)) * thickness * suppressWarnings(as.numeric(comppct_r)),
      weighted_sand = suppressWarnings(as.numeric(sandtotal_r)) * thickness * suppressWarnings(as.numeric(comppct_r)),
      weighted_om = suppressWarnings(as.numeric(om_r)) * thickness * suppressWarnings(as.numeric(comppct_r)),
      weighted_bd = suppressWarnings(as.numeric(dbthirdbar_r)) * thickness * suppressWarnings(as.numeric(comppct_r)),
      depth_range = paste(top_depth, "-", bottom_depth, "cm")
    ) %>%
    dplyr::filter(thickness > 0)

  if (nrow(ungrouped) == 0) return(NULL)

  grouped <- ungrouped %>%
    dplyr::group_by(depth_range) %>%
    dplyr::summarize(
      clay_pct = sum(weighted_clay, na.rm = TRUE) / sum(thickness * suppressWarnings(as.numeric(comppct_r)), na.rm = TRUE),
      sand_pct = sum(weighted_sand, na.rm = TRUE) / sum(thickness * suppressWarnings(as.numeric(comppct_r)), na.rm = TRUE),
      silt_pct = 100 - clay_pct - sand_pct,
      om_pct = sum(weighted_om, na.rm = TRUE) / sum(thickness * suppressWarnings(as.numeric(comppct_r)), na.rm = TRUE),
      bulk_density = sum(weighted_bd, na.rm = TRUE) / sum(thickness * suppressWarnings(as.numeric(comppct_r)), na.rm = TRUE),
      .groups = "drop"
    )
  grouped
}

build_simple_fallback_profile <- function(point_sf, point_id, lat, lon, comp_tbl = NULL, log_file = NULL) {
  soil_helper_log(log_file, "INFO", "SSURGO_FALLBACK", "Starting weighted-layer fallback profile build", point_id = point_id)
  spatial_hit <- robust_SDA_spatialQuery_alderman(point_sf, what = 'mukey')
  mukeys <- if (!is.null(spatial_hit) && 'mukey' %in% names(spatial_hit)) unique(na.omit(as.character(spatial_hit$mukey))) else character(0)
  if (length(mukeys) == 0) {
    soil_helper_log(log_file, "ERROR", "SSURGO_FALLBACK", "No mukey found from fallback spatial query", point_id = point_id)
    return(NULL)
  }
  soil_helper_log(log_file, "INFO", "SSURGO_FALLBACK",
                  sprintf("Fallback spatial query returned %d mukey(s): %s", length(mukeys), paste(head(mukeys, 5), collapse = ", ")),
                  point_id = point_id)

  q_bedrock <- paste0("SELECT mukey, brockdepmin FROM muaggatt WHERE mukey IN ", sql_in_from_values(mukeys))
  bedrock_data <- robust_SDA_query_alderman(q_bedrock)
  bedrock_depth <- 200
  if (!is.null(bedrock_data) && nrow(bedrock_data) > 0 && 'brockdepmin' %in% names(bedrock_data)) {
    bd_vals <- suppressWarnings(as.numeric(bedrock_data$brockdepmin))
    bd_vals <- bd_vals[is.finite(bd_vals)]
    if (length(bd_vals) > 0) bedrock_depth <- min(bd_vals)
  } else {
    soil_helper_log(log_file, "WARN", "SSURGO_FALLBACK",
                    "Bedrock query returned no usable rows; defaulting bedrock depth to 200 cm",
                    point_id = point_id)
  }
  if (!is.finite(bedrock_depth) || is.na(bedrock_depth) || bedrock_depth <= 0) bedrock_depth <- 200

  all_layers <- list(
    '0-5cm' = c(0, 5), '5-20cm' = c(5, 20), '20-35cm' = c(20, 35),
    '35-50cm' = c(35, 50), '50-65cm' = c(50, 65), '65-80cm' = c(65, 80),
    '80-95cm' = c(80, 95), '95-110cm' = c(95, 110), '110-125cm' = c(110, 125),
    '125-140cm' = c(125, 140), '140-155cm' = c(140, 155), '155-170cm' = c(155, 170),
    '170-185cm' = c(170, 185), '185-200cm' = c(185, 200)
  )
  valid_layers <- all_layers[sapply(all_layers, function(x) x[1] < bedrock_depth)]
  if (length(valid_layers) > 0) {
    last <- length(valid_layers)
    if (valid_layers[[last]][2] > bedrock_depth) valid_layers[[last]][2] <- bedrock_depth
    names(valid_layers)[last] <- paste0(valid_layers[[last]][1], '-', valid_layers[[last]][2], 'cm')
  } else {
    valid_layers <- list(c(0, bedrock_depth))
    names(valid_layers) <- paste0('0-', bedrock_depth, 'cm')
  }
  soil_helper_log(log_file, "INFO", "SSURGO_FALLBACK",
                  sprintf("Fallback using bedrock depth %.1f cm and %d target layer(s)", bedrock_depth, length(valid_layers)),
                  point_id = point_id)

  q_soil <- paste0(
    "SELECT component.mukey, component.cokey, component.comppct_r, ",
    "chorizon.hzdept_r, chorizon.hzdepb_r, chorizon.claytotal_r, ",
    "chorizon.sandtotal_r, chorizon.om_r, chorizon.dbthirdbar_r ",
    "FROM component INNER JOIN chorizon ON component.cokey = chorizon.cokey ",
    "WHERE component.mukey IN ", sql_in_from_values(mukeys)
  )
  props <- robust_SDA_query_alderman(q_soil)
  if (is.null(props)) {
    soil_helper_log(log_file, "ERROR", "SSURGO_FALLBACK", "Fallback property query returned NULL", point_id = point_id)
    return(NULL)
  }
  if (nrow(props) == 0) {
    soil_helper_log(log_file, "ERROR", "SSURGO_FALLBACK", "Fallback property query returned 0 rows", point_id = point_id)
    return(NULL)
  }
  props <- as.data.frame(props)
  soil_helper_log(log_file, "INFO", "SSURGO_FALLBACK",
                  sprintf("Fallback property query returned %d horizon row(s)", nrow(props)),
                  point_id = point_id)

  results_list <- lapply(names(valid_layers), function(layer_name) {
    d <- valid_layers[[layer_name]]
    cp <- calculate_soil_properties_fallback(props, d[1], d[2])
    if (!is.null(cp) && nrow(cp) > 0) cp$depth_range <- layer_name
    cp
  })
  layer_ok <- vapply(results_list, function(x) !is.null(x) && nrow(x) > 0, logical(1))
  results_list <- results_list[layer_ok]
  if (length(results_list) == 0) {
    soil_helper_log(log_file, "ERROR", "SSURGO_FALLBACK",
                    sprintf("Fallback produced 0 valid aggregated layers out of %d target layer(s)", length(valid_layers)),
                    point_id = point_id)
    return(NULL)
  }

  results_df <- dplyr::bind_rows(results_list)
  soil_helper_log(log_file, "INFO", "SSURGO_FALLBACK",
                  sprintf("Fallback aggregated %d valid layer row(s)", nrow(results_df)),
                  point_id = point_id)
  dom_comp <- if (!is.null(comp_tbl) && nrow(comp_tbl) > 0) choose_dominant_component(comp_tbl) else data.frame(cokey = NA_character_, mukey = mukeys[1], compname = 'SSURGO', comppct_r = NA_real_, hydgrp = NA_character_, slope_r = NA_real_, drainage = NA_character_, albedodry_r = NA_real_, stringsAsFactors = FALSE)

  layers <- results_df %>%
    dplyr::mutate(
      depth_num = as.numeric(sub('.*-', '', sub('cm', '', depth_range))),
      soc = om_pct / 1.724,
      SLLL = clip01(soil_ptf_saxton_slll(silt_pct, clay_pct, soc, coalesce_num(bulk_density, 1.4), 0)),
      SDUL = clip01(soil_ptf_saxton_sdul(silt_pct, clay_pct, soc, coalesce_num(bulk_density, 1.4), 0)),
      SSAT = clip01(soil_ptf_saxton_ssat(silt_pct, clay_pct, soc, coalesce_num(bulk_density, 1.4), 0)),
      SDUL = pmax(SDUL, SLLL + 0.005),
      SSAT = pmax(SSAT, SDUL + 0.01),
      SLB = as.integer(round(depth_num)),
      SLMH = '-99',
      SRGF = 1.0,
      SSKS = soil_ptf_saxton_ssks(SSAT, SDUL, SLLL, 0, coalesce_num(bulk_density, 1.4)),
      SBDM = coalesce_num(bulk_density, 1.4),
      SLOC = coalesce_num(soc, 0),
      SLCL = coalesce_num(clay_pct, -99),
      SLSI = coalesce_num(silt_pct, -99),
      SLCF = 0,
      SLNI = NA_real_,
      SLHW = NA_real_,
      SLHB = NA_real_,
      SCEC = NA_real_,
      SADC = NA_real_
    ) %>%
    dplyr::select(SLB, SLMH, SLLL, SDUL, SSAT, SRGF, SSKS, SBDM, SLOC,
           SLCL, SLSI, SLCF, SLNI, SLHW, SLHB, SCEC, SADC) %>%
    dplyr::arrange(SLB)

  if (nrow(layers) == 0) {
    soil_helper_log(log_file, "ERROR", "SSURGO_FALLBACK", "Fallback layers table is empty after DSSAT formatting", point_id = point_id)
    return(NULL)
  }

  soil_helper_log(log_file, "INFO", "SSURGO_FALLBACK",
                  sprintf("Fallback profile retained %d DSSAT layer(s): %s",
                          nrow(layers), paste(layers$SLB, collapse = ",")),
                  point_id = point_id)

  slu1 <- soil_ptf_slu1(layers$SSAT, layers$SLLL, layers$SLB)
  if (!is.finite(slu1) || is.na(slu1)) slu1 <- 6

  list(
    profile_id = as.character(point_id),
    latitude = as.numeric(lat),
    longitude = as.numeric(lon),
    site = as.character(point_id),
    country = 'USA',
    scs_family = '',
    scom = 'SC',
    salb = coalesce_num(suppressWarnings(as.numeric(dom_comp$albedodry_r)), 0.13),
    slu1 = round(slu1, 1),
    sldr = sldr_from_drainage(dom_comp$drainage %||% 'Well drained'),
    slro = soil_ptf_curve_number(
      safe_first_non_na(suppressWarnings(as.numeric(dom_comp$slope_r)), 0),
      ifelse(substr(as.character(dom_comp$hydgrp), 1, 1) %in% c('A', 'B', 'C', 'D'),
             match(substr(as.character(dom_comp$hydgrp), 1, 1), c('A', 'B', 'C', 'D')),
             2),
      layers$SSKS, layers$SLB
    ),
    slnf = 1,
    slpf = 1,
    smhb = 'IB001',
    smpx = 'IB001',
    smke = 'IB001',
    layers = layers,
    metadata = data.frame(
      ID = as.character(point_id),
      SOIL_ID = as.character(point_id),
      mukey = paste(mukeys, collapse = ';'),
      cokey = as.character(dom_comp$cokey[1] %||% NA_character_),
      compname = as.character(dom_comp$compname[1] %||% 'SSURGO'),
      comppct_r = suppressWarnings(as.numeric(dom_comp$comppct_r[1])),
      latitude = as.numeric(lat),
      longitude = as.numeric(lon),
      stringsAsFactors = FALSE
    )
  )
}

# ---- Point processor ---------------------------------------------------------
process_one_ssurgo_point <- function(point_data_row, id_col, lat_col, long_col, output_dir_individual, log_file = NULL) {
  point_id <- as.character(point_data_row[[id_col]])
  lat <- as.numeric(point_data_row[[lat_col]])
  lon <- as.numeric(point_data_row[[long_col]])
  point_sf <- sf::st_as_sf(point_data_row, coords = c(long_col, lat_col), crs = 4326)

  soil_helper_log(log_file, "INFO", "SSURGO_POINT", "Starting SSURGO point processing", point_id = point_id)
  profile <- NULL
  comp_tbl <- get_point_component_table(point_sf, point_id = point_id, log_file = log_file)
  if (!is.null(comp_tbl) && nrow(comp_tbl) > 0) {
    component_row <- choose_dominant_component(comp_tbl)
    soil_helper_log(log_file, "INFO", "SSURGO_POINT",
                    sprintf("Selected dominant component cokey=%s compname=%s comppct_r=%s",
                            as.character(component_row$cokey[1]),
                            as.character(component_row$compname[1]),
                            as.character(component_row$comppct_r[1])),
                    point_id = point_id)
    horizon_tbl <- get_component_horizons(component_row$cokey[1], point_id = point_id, log_file = log_file)
    if (!is.null(horizon_tbl) && nrow(horizon_tbl) > 0) {
      profile <- build_dssat_profile_from_component(component_row, horizon_tbl, point_id, lat, lon, log_file = log_file)
      if (!is.null(profile)) {
        soil_helper_log(log_file, "INFO", "SSURGO_POINT", "Built profile from dominant component/horizon data", point_id = point_id)
      } else {
        soil_helper_log(log_file, "WARN", "SSURGO_POINT", "Dominant component profile builder returned NULL", point_id = point_id)
      }
    } else {
      soil_helper_log(log_file, "WARN", "SSURGO_POINT", "No horizons available for dominant component; skipping to fallback", point_id = point_id)
    }
  } else {
    soil_helper_log(log_file, "WARN", "SSURGO_POINT", "No component table available; skipping to fallback", point_id = point_id)
  }

  if (is.null(profile)) {
    soil_helper_log(log_file, "WARN", "SSURGO_POINT", "Falling back to weighted-layer SSURGO profile logic", point_id = point_id)
    profile <- build_simple_fallback_profile(point_sf, point_id, lat, lon, comp_tbl, log_file = log_file)
  }
  if (is.null(profile)) {
    soil_helper_log(log_file, "ERROR", "SSURGO_POINT", "No soil profile could be generated for this point", point_id = point_id)
    return(NULL)
  }

  write_dssat_soil_file(profile, output_dir_individual)
  soil_helper_log(log_file, "INFO", "SSURGO_POINT", sprintf("Wrote soil profile to %s", file.path(output_dir_individual, paste0(point_id, ".SOL"))), point_id = point_id)
  profile$metadata
}

# ---- Main entry point expected by pipeline ----------------------------------
#' Process soils via SSURGO dominant-component/Saxton-Rawls PTF logic
#'
#' @export
process_soils_ssurgo_alderman <- function(grid_points, output_dir_csv, output_dir_individual, n_cores,
                                          id_col, lat_col, long_col, format_sql_func, log_file = NULL) {
  message("Starting SSURGO Processing (dominant component/measured tension fallbacks)...")
  soil_helper_log(log_file, "INFO", "SSURGO_MAIN", "Starting SSURGO Processing (dominant component/measured tension fallbacks)")

  grid_df <- grid_points %>% sf::st_drop_geometry() %>% as.data.frame()
  if (!lat_col %in% names(grid_df)) grid_df[[lat_col]] <- sf::st_coordinates(grid_points)[, 2]
  if (!long_col %in% names(grid_df)) grid_df[[long_col]] <- sf::st_coordinates(grid_points)[, 1]
  grid_df[[id_col]] <- as.character(grid_df[[id_col]])

  dir.create(dirname(output_dir_csv), recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir_individual, recursive = TRUE, showWarnings = FALSE)

  existing_files <- if (dir.exists(output_dir_individual)) tools::file_path_sans_ext(list.files(output_dir_individual, pattern = "\\.SOL$")) else character(0)
  missing_mask <- !(grid_df[[id_col]] %in% existing_files)
  points_to_process <- grid_df[missing_mask, , drop = FALSE]

  message(sprintf("Resume Check: Found %d existing profiles. Processing %d remaining points.",
                  nrow(grid_df) - nrow(points_to_process), nrow(points_to_process)))
  soil_helper_log(log_file, "INFO", "SSURGO_MAIN", sprintf("Resume check: %d existing profiles, %d remaining points.",
                  nrow(grid_df) - nrow(points_to_process), nrow(points_to_process)))

  if (nrow(points_to_process) == 0) {
    # ensure mapping CSV exists / is refreshed
    mapping_df <- grid_df %>%
      dplyr::transmute(
        !!id_col := .data[[id_col]],
        SOIL_ID = .data[[id_col]],
        latitude = .data[[lat_col]],
        longitude = .data[[long_col]]
      )
    readr::write_csv(mapping_df, output_dir_csv)
    message("All soil profiles already exist. SSURGO processing skipped.")
    soil_helper_log(log_file, "INFO", "SSURGO_MAIN", "All soil profiles already exist. SSURGO processing skipped.")
    return(TRUE)
  }

  process_point_wrapper <- function(row_df) {
    tryCatch(
      process_one_ssurgo_point(row_df, id_col, lat_col, long_col, output_dir_individual, log_file = log_file),
      error = function(e) {
        point_id <- as.character(row_df[[id_col]][1])
        message(sprintf("SSURGO point %s failed: %s", point_id, e$message))
        soil_helper_log(log_file, "ERROR", "SSURGO_POINT", conditionMessage(e), point_id = point_id)
        NULL
      }
    )
  }

  total_points <- nrow(points_to_process)
  chunk_size <- 1000
  num_chunks <- ceiling(total_points / chunk_size)
  message(sprintf("Processing %d points in %d chunks...", total_points, num_chunks))
  soil_helper_log(log_file, "INFO", "SSURGO_MAIN", sprintf("Processing %d points in %d chunks.", total_points, num_chunks))

  all_metadata <- list()
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(max(1, n_cores))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages({library(soilDB); library(sf); library(dplyr); library(tidyr); library(stringr); library(readr)}))
    parallel::clusterExport(cl,
                  varlist = c(
                    "robust_SDA_query_alderman", "robust_SDA_spatialQuery_alderman", "%||%",
                    "safe_first_non_na", "coalesce_num", "clip01", "sanitize_char", "sql_in_from_values",
                    "soil_ptf_saxton_slll", "soil_ptf_saxton_sdul", "soil_ptf_saxton_ssat",
                    "soil_ptf_saxton_ssks", "soil_ptf_slu1", "soil_ptf_nrcs_hsg",
                    "soil_ptf_curve_number", "sldr_from_drainage",
                    "get_point_component_table", "get_component_horizons", "choose_dominant_component",
                    "calc_sbdm", "is_coarse_texture", "calc_sdul_measured", "calc_ssat_measured",
                    "calculate_soil_properties_fallback", "build_simple_fallback_profile",
                    "build_dssat_profile_from_component", "format_dssat_decimal", "format_dssat_numeric",
                    "write_dssat_soil_file", "append_log_line", "soil_helper_log", "process_one_ssurgo_point", "process_point_wrapper",
                    "id_col", "lat_col", "long_col", "output_dir_individual", "log_file"
                  ),
                  envir = environment())
  } else {
    cl <- max(1, n_cores)
  }

  for (i in seq_len(num_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, total_points)
    message(sprintf("  > Chunk %d/%d (Points %d-%d)", i, num_chunks, start_idx, end_idx))
    soil_helper_log(log_file, "INFO", "SSURGO_CHUNK", sprintf("Chunk %d/%d (Points %d-%d)", i, num_chunks, start_idx, end_idx))
    chunk_data <- points_to_process[start_idx:end_idx, , drop = FALSE]
    chunk_list <- split(chunk_data, seq_len(nrow(chunk_data)))
    chunk_results <- pbapply::pblapply(chunk_list, process_point_wrapper, cl = cl)
    valid_results <- chunk_results[!vapply(chunk_results, is.null, logical(1))]
    if (length(valid_results) > 0) all_metadata[[length(all_metadata) + 1]] <- dplyr::bind_rows(valid_results)
    soil_helper_log(log_file, "INFO", "SSURGO_CHUNK", sprintf("Chunk %d produced %d valid profiles.", i, length(valid_results)))
    rm(chunk_data, chunk_list, chunk_results, valid_results)
    gc()
  }

  new_mapping <- if (length(all_metadata) > 0) dplyr::bind_rows(all_metadata) else data.frame()
  existing_files_after <- if (dir.exists(output_dir_individual)) tools::file_path_sans_ext(list.files(output_dir_individual, pattern = "\\.SOL$")) else character(0)
  if (nrow(new_mapping) == 0) {
    if (length(existing_files_after) > 0) {
      soil_helper_log(log_file, "WARN", "SSURGO_MAIN",
                      sprintf("No new soil profiles were created in this run, but %d existing individual .SOL files are already present. Continuing without abort.",
                              length(existing_files_after)))
      if (!file.exists(output_dir_csv)) {
        base_mapping_fallback <- grid_df %>%
          dplyr::transmute(
            !!id_col := .data[[id_col]],
            SOIL_ID = .data[[id_col]],
            latitude = .data[[lat_col]],
            longitude = .data[[long_col]]
          )
        readr::write_csv(base_mapping_fallback, output_dir_csv)
        soil_helper_log(log_file, "WARN", "SSURGO_MAIN",
                        sprintf("No existing mapping CSV found. Wrote fallback mapping CSV without component metadata: %s", output_dir_csv))
      } else {
        soil_helper_log(log_file, "INFO", "SSURGO_MAIN",
                        sprintf("Keeping existing soil mapping CSV unchanged: %s", output_dir_csv))
      }
      message("No new soil profiles were created in this run, but existing .SOL files were found. Continuing without abort.")
      return(TRUE)
    }
    soil_helper_log(log_file, "ERROR", "SSURGO_MAIN", "SSURGO_Alderman did not generate any individual .SOL files. Check SDA access and per-point query logic.")
    stop("SSURGO_Alderman did not generate any individual .SOL files. Check SDA access and per-point query logic.")
  }
  base_mapping <- grid_df %>%
    dplyr::transmute(
      !!id_col := .data[[id_col]],
      SOIL_ID = .data[[id_col]],
      latitude = .data[[lat_col]],
      longitude = .data[[long_col]]
    )

  if (nrow(new_mapping) > 0) {
    final_mapping <- base_mapping %>%
      dplyr::left_join(new_mapping, by = setNames(id_col, id_col), suffix = c("", ".new")) %>%
      dplyr::mutate(
        SOIL_ID = if ("SOIL_ID.new" %in% names(.)) dplyr::coalesce(.data[["SOIL_ID.new"]], SOIL_ID) else SOIL_ID,
        latitude = if ("latitude.new" %in% names(.)) dplyr::coalesce(.data[["latitude.new"]], latitude) else latitude,
        longitude = if ("longitude.new" %in% names(.)) dplyr::coalesce(.data[["longitude.new"]], longitude) else longitude
      ) %>%
      dplyr::select(dplyr::any_of(c(id_col, "SOIL_ID", "mukey", "cokey", "compname", "comppct_r", "latitude", "longitude"))) %>%
      dplyr::distinct(.data[[id_col]], .keep_all = TRUE)
  } else {
    final_mapping <- base_mapping
  }

  readr::write_csv(final_mapping, output_dir_csv)
  soil_helper_log(log_file, "INFO", "SSURGO_MAIN", sprintf("Wrote soil mapping CSV: %s", output_dir_csv))
  message("SSURGO Processing Complete.")
  soil_helper_log(log_file, "INFO", "SSURGO_MAIN", "SSURGO Processing Complete.")
  TRUE
}
